/* voneshot.c -- C port of voneshot.py, the one-shot certificate verifier of
 * github.com/AndrewVSutherland/OneShotPrimalityProofs, with a timing harness.  Same
 * checks in the same order, same x-only ladder formulas, same divide-and-conquer order
 * tree.  Field arithmetic is generic Montgomery multiplication on N-limb residues (GMP
 * mpn), nothing specific to any p, and a Miller-Rabin round is timed both with
 * mpz_powm and in that same arithmetic.  Differential test vs the Python verifier:
 * python3 ecpp/test_voneshot.py.
 *
 *   make -C ecpp voneshot
 *   ecpp/voneshot p A x0 m [q1 ... qk]      # print True/False, exit 0/1, like voneshot.py
 *   ecpp/voneshot --bench cert.txt ...      # time each cert file ("p A x0 m q1 ... qk")
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

#define MAXN 64                         /* p up to 4096 bits */
#define MAXQ 4096

static int N;                           /* limbs of p */
static mp_limb_t P[MAXN], PINV;         /* p and -p^{-1} mod 2^64 */
static mp_limb_t AM[MAXN], ONE[MAXN];   /* A and 1 in Montgomery form */
static mpz_t PZ;
static int64_t nmul;                    /* field multiplications */
static double ph[5];                    /* phase timestamps of the last verify() */

typedef struct { mp_limb_t X[MAXN], Z[MAXN]; } pt;

static double now(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* ---------------------------------------------------------------------------
 * Montgomery arithmetic mod p, R = 2^(64N).  REDC as in GMP's redc_1: the carry of
 * the i-th addmul lands in the (now zero) limb t[i] and is added in at the end.
 * ------------------------------------------------------------------------- */
static void redc(mp_limb_t *r, mp_limb_t *t)
{
	for (int i = 0; i < N; i++)
		t[i] = mpn_addmul_1(t + i, P, N, t[i] * PINV);
	if (mpn_add_n(r, t + N, t, N) || mpn_cmp(r, P, N) >= 0)
		mpn_sub_n(r, r, P, N);
}

static void fmul(mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b)
{
	mp_limb_t t[2 * MAXN];
	if (a == b)
		mpn_sqr(t, a, N);
	else
		mpn_mul_n(t, a, b, N);
	redc(r, t);
	nmul++;
}

static void fadd(mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b)
{
	if (mpn_add_n(r, a, b, N) || mpn_cmp(r, P, N) >= 0)
		mpn_sub_n(r, r, P, N);
}

static void fsub(mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b)
{
	if (mpn_sub_n(r, a, b, N))
		mpn_add_n(r, r, P, N);
}

static void to_mont(mp_limb_t *r, const mpz_t x)      /* x in [0, p) */
{
	mpz_t t;
	mpz_init(t);
	mpz_mul_2exp(t, x, 64 * N);
	mpz_mod(t, t, PZ);
	memset(r, 0, N * sizeof(mp_limb_t));
	mpz_export(r, NULL, -1, sizeof(mp_limb_t), 0, 0, t);
	mpz_clear(t);
}

static void mont_setup(const mpz_t p, const mpz_t A)
{
	mpz_t one;
	N = mpz_size(p);
	memset(P, 0, sizeof P);
	mpz_export(P, NULL, -1, sizeof(mp_limb_t), 0, 0, p);
	mp_limb_t inv = 1;                  /* Newton: p^{-1} mod 2^64 */
	for (int i = 0; i < 6; i++)
		inv *= 2 - P[0] * inv;
	PINV = -inv;
	mpz_set(PZ, p);
	to_mont(AM, A);
	mpz_init_set_ui(one, 1);
	to_mont(ONE, one);
	mpz_clear(one);
}

/* gcd(x, p) == 1; the Montgomery factor R is a unit, so the representative will do */
static int unit(const mp_limb_t *x)
{
	mpz_t g, xz;
	mpz_init(g);
	mpz_gcd(g, mpz_roinit_n(xz, x, N), PZ);
	int ok = mpz_cmp_ui(g, 1) == 0;
	mpz_clear(g);
	return ok;
}

/* ---------------------------------------------------------------------------
 * Montgomery x-only arithmetic, the formulas of voneshot.py (6 multiplications each).
 * ------------------------------------------------------------------------- */
static void xdbl(pt *r, const pt *q)    /* X' = (X^2-Z^2)^2 ; Z' = 4XZ(X^2+AXZ+Z^2) */
{
	mp_limb_t XX[MAXN], ZZ[MAXN], XZ[MAXN], t[MAXN], u[MAXN];
	fmul(XX, q->X, q->X);
	fmul(ZZ, q->Z, q->Z);
	fmul(XZ, q->X, q->Z);
	fsub(t, XX, ZZ);
	fmul(u, AM, XZ);
	fadd(u, u, XX);
	fadd(u, u, ZZ);
	fadd(XZ, XZ, XZ);
	fadd(XZ, XZ, XZ);
	fmul(r->X, t, t);
	fmul(r->Z, XZ, u);
}

/* a + b given d = a - b (d must not alias r) */
static void xadd(pt *r, const pt *a, const pt *b, const pt *d)
{
	mp_limb_t s[MAXN], e[MAXN], t1[MAXN], t2[MAXN];
	fsub(t1, a->X, a->Z);
	fadd(t2, b->X, b->Z);
	fmul(s, t1, t2);                    /* a */
	fadd(t1, a->X, a->Z);
	fsub(t2, b->X, b->Z);
	fmul(e, t1, t2);                    /* b */
	fadd(t1, s, e);
	fsub(t2, s, e);
	fmul(t1, t1, t1);
	fmul(t2, t2, t2);
	fmul(r->X, d->Z, t1);
	fmul(r->Z, d->X, t2);
}

static void ladder(pt *r, const mpz_t k, const pt *p)
{
	if (mpz_sgn(k) == 0) {
		memcpy(r->X, ONE, N * sizeof(mp_limb_t));
		memset(r->Z, 0, N * sizeof(mp_limb_t));
		return;
	}
	if (mpz_cmp_ui(k, 1) == 0) {
		*r = *p;
		return;
	}
	pt d = *p, R0 = *p, R1;
	xdbl(&R1, p);
	for (long i = (long)mpz_sizeinbase(k, 2) - 2; i >= 0; i--) {
		if (mpz_tstbit(k, i)) {
			xadd(&R0, &R0, &R1, &d);
			xdbl(&R1, &R1);
		} else {
			xadd(&R1, &R0, &R1, &d);
			xdbl(&R0, &R0);
		}
	}
	*r = R0;
}

/* (m/q)P != O for every q, given Q = (m/R)P, by the divide-and-conquer tree */
static int check_orders(const pt *q, mpz_t *pr, int t)
{
	if (t == 0)
		return 1;
	if (t == 1)
		return unit(q->Z);
	int mid = t / 2;
	mpz_t hL, hR;
	mpz_init_set_ui(hL, 1);
	mpz_init_set_ui(hR, 1);
	for (int i = 0; i < mid; i++)
		mpz_mul(hL, hL, pr[i]);
	for (int i = mid; i < t; i++)
		mpz_mul(hR, hR, pr[i]);
	pt L, R;
	ladder(&L, hL, q);
	ladder(&R, hR, q);
	mpz_clear(hL);
	mpz_clear(hR);
	return check_orders(&L, pr + mid, t - mid) && check_orders(&R, pr, mid);
}

/* ---------------------------------------------------------------------------
 * The verifier (voneshot.py's verify, step for step).
 * ------------------------------------------------------------------------- */
static uint32_t *primes;
static int64_t nprimes;

static void sieve(uint64_t lim)         /* odd-only byte sieve, primes <= lim */
{
	uint64_t h = lim / 2 + 1;
	uint8_t *c = calloc(h, 1);          /* c[i] <=> 2i+1 composite */
	primes = realloc(primes, (lim / 2 + 2) * sizeof(uint32_t));
	nprimes = 0;
	if (lim >= 2)
		primes[nprimes++] = 2;
	for (uint64_t i = 1; i < h; i++) {
		if (c[i])
			continue;
		uint64_t q = 2 * i + 1;
		if (q > lim)
			break;
		primes[nprimes++] = (uint32_t)q;
		for (uint64_t j = q * q / 2; j < h; j += q)
			c[j] = 1;
	}
	free(c);
}

static int verify(const mpz_t p, const mpz_t A, const mpz_t x0, const mpz_t m,
                  mpz_t *qs, int k)
{
	int ok = 0, nd = 0;
	mpz_t sp, L, t, rough, *div = NULL;
	mpz_inits(sp, L, t, rough, NULL);
	ph[0] = now();
	if (mpz_size(p) > MAXN) {
		fprintf(stderr, "voneshot: p exceeds %d bits (raise MAXN)\n", 64 * MAXN);
		exit(2);
	}
	if (mpz_cmp_ui(p, 3) <= 0 || mpz_even_p(p))
		goto done;
	if (mpz_sgn(A) < 0 || mpz_cmp(A, p) >= 0 || mpz_sgn(x0) < 0 || mpz_cmp(x0, p) >= 0)
		goto done;

	/* (i) size window */
	mpz_sqrt(sp, p);
	mpz_mul_2exp(t, sp, 2);
	mpz_sqrt(t, t);
	mpz_add(L, sp, t);
	mpz_add_ui(L, L, 1);
	if (mpz_cmp(m, L) <= 0)
		goto done;
	mpz_mul_2exp(t, p, 2);
	mpz_sqrt(t, t);
	mpz_add(t, t, p);
	mpz_add_ui(t, t, 1);
	if (mpz_cmp(m, t) > 0)
		goto done;

	/* (ii) nonsingular */
	mpz_mul(t, A, A);
	mpz_sub_ui(t, t, 4);
	mpz_mod(t, t, p);
	mpz_gcd(t, t, p);
	if (mpz_cmp_ui(t, 1) != 0)
		goto done;

	/* (iii) small primes of m by trial division up to n^2, then the q_i */
	uint64_t n = mpz_sizeinbase(p, 2), n2 = n * n, n4 = n2 * n2;
	sieve(n2);
	ph[1] = now();
	div = malloc((64 * mpz_size(m) + k + 1) * sizeof(mpz_t));
	mpz_set(rough, m);
	for (int64_t i = 0; i < nprimes; i++)
		if (mpz_divisible_ui_p(m, primes[i])) {
			mpz_init_set_ui(div[nd++], primes[i]);
			while (mpz_divisible_ui_p(rough, primes[i]))
				mpz_divexact_ui(rough, rough, primes[i]);
		}
	mpz_set_ui(t, n2);
	for (int i = 0; i < k; i++) {
		if (mpz_cmp(qs[i], t) <= 0 || mpz_cmp_ui(qs[i], n4) >= 0)
			goto done;
		mpz_set(t, qs[i]);
	}
	for (int i = 0; i < k; i++) {
		if (!mpz_divisible_p(rough, qs[i]))
			goto done;
		while (mpz_divisible_p(rough, qs[i]))
			mpz_divexact(rough, rough, qs[i]);
		mpz_init_set(div[nd++], qs[i]);
	}
	if (mpz_cmp_ui(rough, 1) != 0 || nd == 0)
		goto done;

	/* (iv) minimality m < L r */
	mpz_mul(t, L, div[0]);
	if (mpz_cmp(m, t) >= 0)
		goto done;
	ph[2] = now();

	/* (v) mP = O, reached as (X:0) with X a unit */
	mont_setup(p, A);
	pt P0, Q;
	to_mont(P0.X, x0);
	memcpy(P0.Z, ONE, N * sizeof(mp_limb_t));
	ladder(&Q, m, &P0);
	if (!mpn_zero_p(Q.Z, N) || !unit(Q.X))
		goto done;
	ph[3] = now();

	/* (vi) (m/q)P != O for every prime q | m */
	mpz_set_ui(t, 1);
	for (int i = 0; i < nd; i++)
		mpz_mul(t, t, div[i]);
	mpz_divexact(t, m, t);
	ladder(&Q, t, &P0);
	ok = check_orders(&Q, div, nd);
	ph[4] = now();
done:
	for (int i = 0; i < nd; i++)
		mpz_clear(div[i]);
	free(div);
	mpz_clears(sp, L, t, rough, NULL);
	return ok;
}

/* ---------------------------------------------------------------------------
 * Miller-Rabin baselines.
 * ------------------------------------------------------------------------- */
static int mr_gmp(const mpz_t p, unsigned long a)      /* one round, mpz_powm */
{
	mpz_t d, x, pm1;
	mpz_inits(d, x, pm1, NULL);
	mpz_sub_ui(pm1, p, 1);
	mp_bitcnt_t s = mpz_scan1(pm1, 0);
	mpz_tdiv_q_2exp(d, pm1, s);
	mpz_set_ui(x, a);
	mpz_powm(x, x, d, p);
	int ok = mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, pm1) == 0;
	for (mp_bitcnt_t i = 1; i < s && !ok; i++) {
		mpz_mul(x, x, x);
		mpz_mod(x, x, p);
		ok = mpz_cmp(x, pm1) == 0;
	}
	mpz_clears(d, x, pm1, NULL);
	return ok;
}

static int mr_mont(const mpz_t p, unsigned long a)     /* one round, 4-bit window, fmul */
{
	mp_limb_t T[16][MAXN], x[MAXN], pm1[MAXN];
	mpz_t d, az, e;
	mpz_inits(d, az, e, NULL);
	mpz_sub_ui(e, p, 1);
	mp_bitcnt_t s = mpz_scan1(e, 0);
	mpz_tdiv_q_2exp(d, e, s);
	to_mont(pm1, e);
	mpz_set_ui(az, a);
	memcpy(T[0], ONE, N * sizeof(mp_limb_t));
	to_mont(T[1], az);
	for (int i = 2; i < 16; i++)
		fmul(T[i], T[i - 1], T[1]);
	long nb = (long)mpz_sizeinbase(d, 2), top = (nb + 3) / 4 * 4;
	memcpy(x, ONE, N * sizeof(mp_limb_t));
	for (long i = top - 4; i >= 0; i -= 4) {
		if (i != top - 4)
			for (int j = 0; j < 4; j++)
				fmul(x, x, x);
		int w = (mpz_tstbit(d, i + 3) << 3) | (mpz_tstbit(d, i + 2) << 2) |
		        (mpz_tstbit(d, i + 1) << 1) | mpz_tstbit(d, i);
		if (w)
			fmul(x, x, T[w]);
	}
	int ok = mpn_cmp(x, ONE, N) == 0 || mpn_cmp(x, pm1, N) == 0;
	for (mp_bitcnt_t i = 1; i < s && !ok; i++) {
		fmul(x, x, x);
		ok = mpn_cmp(x, pm1, N) == 0;
	}
	mpz_clears(d, az, e, NULL);
	return ok;
}

/* ---------------------------------------------------------------------------
 * Driver.
 * ------------------------------------------------------------------------- */
static int parse(int argc, char **argv, mpz_t *v)
{
	for (int i = 0; i < argc; i++)
		if (mpz_init_set_str(v[i], argv[i], 0) != 0)
			return 0;
	return 1;
}

static int cmpd(const void *a, const void *b)
{
	double x = *(const double *)a, y = *(const double *)b;
	return (x > y) - (x < y);
}

static void bench(const char *file)
{
	static char *tok[MAXQ + 4];
	static mpz_t v[MAXQ + 4];
	enum { REPS = 2001 };
	static double tv[REPS], pa[REPS], pb[REPS], pc[REPS], pd[REPS], tg[REPS], tm[REPS];
	FILE *f = fopen(file, "r");
	if (!f) {
		perror(file);
		exit(2);
	}
	char buf[65536];
	size_t len = fread(buf, 1, sizeof buf - 1, f);
	buf[len] = 0;
	fclose(f);
	int c = 0;
	for (char *s = strtok(buf, " \t\r\n"); s && c < MAXQ + 4; s = strtok(NULL, " \t\r\n"))
		tok[c++] = s;
	if (c < 4 || !parse(c, tok, v)) {
		fprintf(stderr, "%s: bad certificate\n", file);
		exit(2);
	}
	nmul = 0;
	if (!verify(v[0], v[1], v[2], v[3], v + 4, c - 4)) {
		fprintf(stderr, "%s: certificate REJECTED\n", file);
		exit(1);
	}
	int64_t vm = nmul;
	nmul = 0;
	mr_mont(v[0], 2);
	int64_t mm = nmul;
	for (int r = 0; r < REPS; r++) {
		verify(v[0], v[1], v[2], v[3], v + 4, c - 4);
		tv[r] = ph[4] - ph[0];
		pa[r] = ph[1] - ph[0];
		pb[r] = ph[2] - ph[1];
		pc[r] = ph[3] - ph[2];
		pd[r] = ph[4] - ph[3];
		double t0 = now();
		mr_gmp(v[0], 2 + r % 64);
		double t1 = now();
		mr_mont(v[0], 2 + r % 64);
		tg[r] = t1 - t0;
		tm[r] = now() - t1;
	}
	double *arr[] = { tv, pa, pb, pc, pd, tg, tm };
	for (int i = 0; i < 7; i++)
		qsort(arr[i], REPS, sizeof(double), cmpd);
	int h = REPS / 2;
	printf("%s: %zu bits, %d listed q_i, m %zu bits\n", file, mpz_sizeinbase(v[0], 2), c - 4,
	       mpz_sizeinbase(v[3], 2));
	printf("  verify      median %8.1f us   (min %.1f)   %ld field mults\n",
	       1e6 * tv[h], 1e6 * tv[0], (long)vm);
	printf("    sieve to n^2 %6.1f   trial div + q_i %6.1f   ladder(m) %6.1f   order tree %6.1f us\n",
	       1e6 * pa[h], 1e6 * pb[h], 1e6 * pc[h], 1e6 * pd[h]);
	printf("  MR round    mpz_powm %6.2f us   same-arith %6.2f us (%ld mults)\n",
	       1e6 * tg[h], 1e6 * tm[h], (long)mm);
	printf("  ratio       verify/powm %.1f   verify/same-arith %.1f   mults %.1f   EC only/same-arith %.1f\n",
	       tv[h] / tg[h], tv[h] / tm[h], (double)vm / mm, (pc[h] + pd[h]) / tm[h]);
	for (int i = 0; i < c; i++)
		mpz_clear(v[i]);
}

int main(int argc, char **argv)
{
	mpz_init(PZ);
	if (argc >= 3 && !strcmp(argv[1], "--bench")) {
		for (int i = 2; i < argc; i++)
			bench(argv[i]);
		return 0;
	}
	static mpz_t v[MAXQ + 4];
	int c = argc - 1;
	if (c < 4 || c > MAXQ + 4) {
		fprintf(stderr, "usage: voneshot p A x0 m [q1 ... qk]  |  voneshot --bench cert.txt ...\n");
		return 2;
	}
	if (!parse(c, argv + 1, v)) {
		fprintf(stderr, "error: p, A, x0, m, q1, ... must be integers\n");
		return 2;
	}
	int ok = verify(v[0], v[1], v[2], v[3], v + 4, c - 4);
	puts(ok ? "True" : "False");
	return ok ? 0 : 1;
}
