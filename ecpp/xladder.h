/* xladder.h -- Montgomery-form field arithmetic and x-only ladder shared by the C
 * certificate verifiers (voneshot.c, vshort2.c).
 *
 * Residues mod an odd p are N-limb Montgomery representatives (R = 2^(64N), GMP mpn),
 * nothing specific to any p.  The x-only formulas compute exactly the residues of the
 * Python verifiers' xdbl/xadd (voneshot.py, vshort2.py), with fewer multiplications:
 *   xdbl: X' = (X+Z)^2 (X-Z)^2 = (X^2-Z^2)^2,
 *         Z' = C ((X-Z)^2 + a24 C) = 4XZ (X^2+AXZ+Z^2),  C = (X+Z)^2-(X-Z)^2 = 4XZ,
 *         a24 = (A+2)/4 mod p  (2S + 3M instead of 6 multiplications);
 *   xadd: as in Python, skipping the multiplication by Zd when Zd = 1 (5 or 6).
 * These are identities in Z/pZ for every odd p, so accept/reject behaviour is that of
 * the Python verifiers, composite p included.  Single-threaded, file-scope state.
 */
#ifndef XLADDER_H
#define XLADDER_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

#define MAXN 64                         /* p up to 4096 bits */

static int N;                           /* limbs of p */
static mp_limb_t P[MAXN], PINV;         /* p and -p^{-1} mod 2^64 */
static mp_limb_t A24[MAXN], ONE[MAXN];  /* (A+2)/4 and 1 in Montgomery form */
static mpz_t PZ;                        /* p (mpz_init'ed by xl_init) */
static int64_t nmul;                    /* field multiplications */

typedef struct { mp_limb_t X[MAXN], Z[MAXN]; } pt;

static inline void xl_init(void)
{
	mpz_init(PZ);
}

static inline double now(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* REDC as in GMP's redc_1: the carry of the i-th addmul lands in the (now zero)
 * limb t[i] and is added in at the end */
static inline void redc(mp_limb_t *r, mp_limb_t *t)
{
	for (int i = 0; i < N; i++)
		t[i] = mpn_addmul_1(t + i, P, N, t[i] * PINV);
	if (mpn_add_n(r, t + N, t, N) || mpn_cmp(r, P, N) >= 0)
		mpn_sub_n(r, r, P, N);
}

static inline void fmul(mp_limb_t *r, mp_limb_t *a, mp_limb_t *b)
{
	mp_limb_t t[2 * MAXN];
	if (a == b)
		mpn_sqr(t, a, N);
	else
		mpn_mul_n(t, a, b, N);
	redc(r, t);
	nmul++;
}

static inline void fadd(mp_limb_t *r, mp_limb_t *a, mp_limb_t *b)
{
	if (mpn_add_n(r, a, b, N) || mpn_cmp(r, P, N) >= 0)
		mpn_sub_n(r, r, P, N);
}

static inline void fsub(mp_limb_t *r, mp_limb_t *a, mp_limb_t *b)
{
	if (mpn_sub_n(r, a, b, N))
		mpn_add_n(r, r, P, N);
}

static inline void to_mont(mp_limb_t *r, mpz_t x)      /* x in [0, p) */
{
	mpz_t t;
	mpz_init(t);
	mpz_mul_2exp(t, x, 64 * N);
	mpz_mod(t, t, PZ);
	memset(r, 0, N * sizeof(mp_limb_t));
	mpz_export(r, NULL, -1, sizeof(mp_limb_t), 0, 0, t);
	mpz_clear(t);
}

/* set up arithmetic mod p (odd, at most MAXN limbs) for the curve coefficient A */
static inline void mont_setup(mpz_t p, mpz_t A)
{
	mpz_t t, f;
	N = mpz_size(p);
	memset(P, 0, sizeof P);
	mpz_export(P, NULL, -1, sizeof(mp_limb_t), 0, 0, p);
	mp_limb_t inv = 1;                  /* Newton: p^{-1} mod 2^64 */
	for (int i = 0; i < 6; i++)
		inv *= 2 - P[0] * inv;
	PINV = -inv;
	mpz_set(PZ, p);
	mpz_init_set_ui(t, 1);
	to_mont(ONE, t);
	mpz_init_set_ui(f, 4);
	mpz_invert(f, f, p);                /* p odd */
	mpz_add_ui(t, A, 2);
	mpz_mul(t, t, f);
	mpz_mod(t, t, p);
	to_mont(A24, t);
	mpz_clears(t, f, NULL);
}

/* gcd(x, p) == 1; the Montgomery factor R is a unit, so the representative will do */
static inline int unit(mp_limb_t *x)
{
	mpz_t g, xz;
	mpz_init(g);
	mpz_gcd(g, mpz_roinit_n(xz, x, N), PZ);
	int ok = mpz_cmp_ui(g, 1) == 0;
	mpz_clear(g);
	return ok;
}

static inline void xdbl(pt *r, pt *q)
{
	mp_limb_t s[MAXN], d[MAXN], ss[MAXN], dd[MAXN], c[MAXN];
	fadd(s, q->X, q->Z);
	fsub(d, q->X, q->Z);
	fmul(ss, s, s);
	fmul(dd, d, d);
	fsub(c, ss, dd);                    /* 4XZ */
	fmul(r->X, ss, dd);
	fmul(s, A24, c);
	fadd(s, s, dd);                     /* X^2 + AXZ + Z^2 */
	fmul(r->Z, c, s);
}

/* a + b given d = a - b (d must not alias r); zd1: d->Z is the Montgomery 1 */
static inline void xadd(pt *r, pt *a, pt *b, pt *d, int zd1)
{
	mp_limb_t s[MAXN], e[MAXN], t1[MAXN], t2[MAXN];
	fsub(t1, a->X, a->Z);
	fadd(t2, b->X, b->Z);
	fmul(s, t1, t2);
	fadd(t1, a->X, a->Z);
	fsub(t2, b->X, b->Z);
	fmul(e, t1, t2);
	fadd(t1, s, e);
	fsub(t2, s, e);
	if (zd1)
		fmul(r->X, t1, t1);
	else {
		fmul(t1, t1, t1);
		fmul(r->X, d->Z, t1);
	}
	fmul(t2, t2, t2);
	fmul(r->Z, d->X, t2);
}

/* (X:Z) of kP, maintaining (R0,R1) = (jP,(j+1)P); k = 0 gives (1:0) */
static inline void ladder(pt *r, mpz_t k, pt *p)
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
	int zd1 = mpn_cmp(p->Z, ONE, N) == 0;
	pt d = *p, R0 = *p, R1;
	xdbl(&R1, p);
	for (long i = (long)mpz_sizeinbase(k, 2) - 2; i >= 0; i--) {
		if (mpz_tstbit(k, i)) {
			xadd(&R0, &R0, &R1, &d, zd1);
			xdbl(&R1, &R1);
		} else {
			xadd(&R1, &R0, &R1, &d, zd1);
			xdbl(&R0, &R0);
		}
	}
	*r = R0;
}

/* (r/q)Q != O (mod every prime divisor of p) for every q in pr[0..t), given
 * Q = (o/r)P with r = prod pr: voneshot.py's divide-and-conquer tree (count split) */
static inline int check_orders(pt *q, mpz_t *pr, int t)
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

/* primes <= lim, odd-only byte sieve; returns the count, *out malloc'ed */
static inline int64_t sieve(uint64_t lim, uint32_t **out)
{
	uint64_t h = lim / 2 + 1;
	uint8_t *c = calloc(h, 1);          /* c[i] <=> 2i+1 composite */
	uint32_t *pr = malloc((lim / 2 + 2) * sizeof(uint32_t));
	int64_t np = 0;
	if (lim >= 2)
		pr[np++] = 2;
	for (uint64_t i = 1; i < h; i++) {
		if (c[i])
			continue;
		uint64_t q = 2 * i + 1;
		if (q > lim)
			break;
		pr[np++] = (uint32_t)q;
		for (uint64_t j = q * q / 2; j < h; j += q)
			c[j] = 1;
	}
	free(c);
	*out = pr;
	return np;
}

/* one Miller-Rabin round with mpz_powm */
static inline int mr_gmp(mpz_t p, unsigned long a)
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

/* one Miller-Rabin round in the arithmetic above (4-bit fixed window), after
 * mont_setup(p, .) */
static inline int mr_mont(mpz_t p, unsigned long a)
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

static inline int cmpd(const void *a, const void *b)
{
	double x = *(const double *)a, y = *(const double *)b;
	return (x > y) - (x < y);
}

#endif
