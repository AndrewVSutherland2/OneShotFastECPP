/* voneshot.c -- C port of voneshot.py, the one-shot certificate verifier of
 * github.com/AndrewVSutherland/OneShotPrimalityProofs, with a timing harness.  Same
 * checks in the same order, same divide-and-conquer order tree, and x-only ladder
 * formulas computing exactly the same residues (xladder.h: generic N-limb Montgomery
 * arithmetic, nothing specific to any p).  A Miller-Rabin round is timed both with
 * mpz_powm and in that same arithmetic.  Differential test vs the Python verifier:
 * python3 ecpp/test_voneshot.py.
 *
 *   make -C ecpp voneshot
 *   ecpp/voneshot p A x0 m [q1 ... qk]      # print True/False, exit 0/1, like voneshot.py
 *   ecpp/voneshot --bench cert.txt ...      # time each cert file ("p A x0 m q1 ... qk")
 */
#include "xladder.h"

#define MAXQ 4096

static double ph[5];                    /* phase timestamps of the last verify() */

/* ---------------------------------------------------------------------------
 * The verifier (voneshot.py's verify, step for step).
 * ------------------------------------------------------------------------- */
static uint32_t *primes;
static int64_t nprimes;

static int verify(mpz_t p, mpz_t A, mpz_t x0, mpz_t m,
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
	free(primes);
	nprimes = sieve(n2, &primes);
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
 * Driver.
 * ------------------------------------------------------------------------- */
static int parse(int argc, char **argv, mpz_t *v)
{
	for (int i = 0; i < argc; i++)
		if (mpz_init_set_str(v[i], argv[i], 0) != 0)
			return 0;
	return 1;
}

static void bench(char *file)
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
	xl_init();
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
