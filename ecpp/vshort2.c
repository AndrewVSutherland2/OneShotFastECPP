/* vshort2.c -- C port of vshort2.py, the verifier for the radical-capped short ECPP
 * format (August 2026, terminal-prime revision), with a timing harness.
 *
 * Two modes; ecpp/test_vshort2.py checks both against vshort2.py on every input.
 *   faithful (default): vshort2.py step for step -- per level a direct ladder for
 *     [o]P = O, then voneshot.py's count-split exactness tree over the primes of o.
 *     Since o = m p_next with m small, the count split multiplies by the big prime
 *     p_next in 1 + ceil(log2 t) ladders per level (t = number of primes of o).
 *   --fast: one big ladder per level.  With R = rad(m) p_next and Q = (o/R)P
 *     (normalized to Z = 1; a non-unit Z(Q) is rejected, as every leaf of the tree
 *     would be), V = p_next Q is the only ladder by p_next; then
 *       (o/p_next)P = rad(m) Q,   (o/q)P = (rad(m)/q) V  for q | m  (tree over the
 *       small primes),   [o]P = rad(m) V.
 *     The leaves are compositions of ladders exactly as in the faithful tree; the one
 *     change is that [o]P = O is reached as rad(m) V instead of by a direct ladder.
 *     A degenerate step modulo some l | p (a summand or the difference point O, or
 *     the difference the 2-torsion point (0,0)) produces (0:0) mod l, which every
 *     later step preserves and which fails "X a unit" just as the direct ladder's
 *     collapse does.  Terminal levels with p_next = 1 use the faithful path.
 *     The radicals g = gcd(P_B mod o, o) are computed without forming the primorial
 *     P_B: products of primes <= B are folded in, in chunks, modulo M = prod o_i (one
 *     pass, ~n bits for any valid chain; capped at 2n+64 bits, else per level mod o
 *     after rejecting o >= B L, which the window o < r L, r <= B, implies).
 *
 *   make -C ecpp vshort2
 *   ecpp/vshort2 [--fast] p0 A0 x0 o0 [A1 x1 o1 ...]   # True/False, exit 0/1
 *   ecpp/vshort2 --bench certs/short2/certs.csv [bits ...]   # time chains
 */
#include <math.h>
#include "xladder.h"

#define MAXSEQ 4096

static double t_prim, t_arith, t_ec;    /* phase accumulators */

/* the small primes of g = rad(m), ascending: trial division by 2 and the odd d with
 * d^2 <= g (g is squarefree with every prime <= B, so odd composites never divide
 * what is left), the cofactor > 1 being prime -- the list vshort2.py builds by
 * running over the primes <= B */
static int factor_rad(mpz_t *out, mpz_t g)
{
	int k = 0;
	mpz_t gg;
	mpz_init_set(gg, g);
	for (unsigned long d = 2; mpz_cmp_ui(gg, 1) > 0; d += (d == 2 ? 1 : 2)) {
		if (mpz_cmp_ui(gg, d * d) < 0)
			break;
		if (mpz_divisible_ui_p(gg, d)) {
			mpz_init_set_ui(out[k++], d);
			mpz_divexact_ui(gg, gg, d);
		}
	}
	if (mpz_cmp_ui(gg, 1) > 0)
		mpz_init_set(out[k++], gg);
	mpz_clear(gg);
	return k;
}

/* --fast: products of consecutive primes <= B, each < 2^64 */
static int64_t prime_words(uint64_t B, uint64_t **out)
{
	uint32_t *pr;
	int64_t np = sieve(B, &pr), nw = 0;
	uint64_t *w = malloc((np + 1) * sizeof(uint64_t)), acc = 1;
	for (int64_t i = 0; i < np; i++) {
		if (acc > UINT64_MAX / pr[i]) {
			w[nw++] = acc;
			acc = 1;
		}
		acc *= pr[i];
	}
	w[nw++] = acc;
	free(pr);
	*out = w;
	return nw;
}

/* --fast: P_B mod o without forming P_B -- word products gathered into chunks of
 * about o's size, folded in mod o (needs o < B L, so o is about n/2 bits) */
static void primorial_mod(mpz_t r, mpz_t o, uint64_t *w, int64_t nw)
{
	size_t lim = mpz_size(o);
	mpz_t c;
	mpz_init_set_ui(c, 1);
	mpz_set_ui(r, 1);
	for (int64_t i = 0; i < nw; i++) {
		mpz_mul_ui(c, c, w[i]);
		if (mpz_size(c) >= lim) {
			mpz_mul(r, r, c);
			mpz_mod(r, r, o);
			mpz_set_ui(c, 1);
		}
	}
	mpz_mul(r, r, c);
	mpz_mod(r, r, o);
	mpz_clear(c);
}

/* the exactness and [o]P = O checks of one level, after mont_setup(p, A) */
static int level_faithful(pt *P0, mpz_t o, mpz_t *div, int nd)
{
	pt Q;
	mpz_t t;
	ladder(&Q, o, P0);
	if (!mpn_zero_p(Q.Z, N) || !unit(Q.X))
		return 0;
	mpz_init_set_ui(t, 1);
	for (int i = 0; i < nd; i++)
		mpz_mul(t, t, div[i]);
	mpz_divexact(t, o, t);
	ladder(&Q, t, P0);
	mpz_clear(t);
	return check_orders(&Q, div, nd);
}

/* --fast: div[0..ns) the small primes, p_next > 1 */
static int level_fast(pt *P0, mpz_t o, mpz_t *div, int ns, mpz_t pn)
{
	pt Q, V, W;
	mpz_t rad, t, zi, xr, zr;           /* xr, zr: read-only views, never cleared */
	int ok = 0;
	mpz_inits(rad, t, zi, NULL);
	mpz_set_ui(rad, 1);
	for (int i = 0; i < ns; i++)
		mpz_mul(rad, rad, div[i]);
	mpz_mul(t, rad, pn);
	mpz_divexact(t, o, t);
	ladder(&Q, t, P0);                  /* Q = (o/R) P, R = rad(m) p_next */
	if (!unit(Q.Z))
		goto done;
	mpz_invert(zi, mpz_roinit_n(zr, Q.Z, N), PZ);      /* Montgomery reps: */
	mpz_mul(t, mpz_roinit_n(xr, Q.X, N), zi);          /* X R / (Z R) = X/Z */
	mpz_mod(t, t, PZ);
	to_mont(Q.X, t);
	memcpy(Q.Z, ONE, N * sizeof(mp_limb_t));
	ladder(&W, rad, &Q);                /* (o/p_next) P */
	if (!unit(W.Z))
		goto done;
	ladder(&V, pn, &Q);                 /* (o/rad(m)) P: the one big ladder */
	if (!check_orders(&V, div, ns))     /* (o/q) P, q | m */
		goto done;
	ladder(&W, rad, &V);                /* o P */
	ok = mpn_zero_p(W.Z, N) && unit(W.X);
done:
	mpz_clears(rad, t, zi, NULL);
	return ok;
}

static int verify(mpz_t *seq, int len, int fast)
{
	int ok = 0, nd = 0;
	mpz_t p, Pr, g, m, pn, oo, L, t, B2, div[160];
	double t0;
	if (len < 4 || (len - 1) % 3 != 0)
		return 0;
	for (int i = 0; i < len; i++)
		if (mpz_sgn(seq[i]) < 0)
			return 0;
	if (mpz_cmp_ui(seq[0], 5) < 0 || mpz_even_p(seq[0]))
		return 0;
	if (mpz_size(seq[0]) > MAXN) {
		fprintf(stderr, "vshort2: p0 exceeds %d bits (raise MAXN)\n", 64 * MAXN);
		exit(2);
	}
	uint64_t n = mpz_sizeinbase(seq[0], 2);
	double lg = log2((double)n);
	uint64_t B = (uint64_t)ceil((double)n * (double)n / lg);
	uint64_t radlim = (uint64_t)ceil((double)n / lg);
	int levels = (len - 1) / 3;
	if ((uint64_t)levels > n)
		return 0;
	for (int i = 3; i < len; i += 3)
		if (mpz_cmp_ui(seq[i], 2) < 0 || mpz_sizeinbase(seq[i], 2) > 2 * n + 2)
			return 0;

	mpz_inits(p, Pr, g, m, pn, oo, L, t, B2, NULL);
	uint64_t *words = NULL;
	int64_t nw = 0;
	int onepass = 0;
	t0 = now();
	if (fast) {
		nw = prime_words(B, &words);
		size_t mb = 0;
		for (int i = 3; i < len; i += 3)
			mb += mpz_sizeinbase(seq[i], 2);
		if (mb <= 2 * n + 64) {         /* every valid chain: sum bits(o_i) <~ n + O(k log n) */
			mpz_set_ui(t, 1);
			for (int i = 3; i < len; i += 3)
				mpz_mul(t, t, seq[i]);
			primorial_mod(Pr, t, words, nw);    /* P_B mod prod o_i */
			onepass = 1;
		}
	} else
		mpz_primorial_ui(Pr, B);
	t_prim += now() - t0;
	mpz_set_ui(B2, B);
	mpz_mul(B2, B2, B2);
	mpz_set(p, seq[0]);
	for (int i = 1; i < len; i += 3) {
		mpz_t *A = &seq[i], *x = &seq[i + 1], *o = &seq[i + 2];
		t0 = now();
		for (int j = 0; j < nd; j++)
			mpz_clear(div[j]);
		nd = 0;
		if (mpz_cmp_ui(p, 3) < 0 || mpz_even_p(p))
			goto done;
		if (mpz_cmp(*A, p) >= 0 || mpz_cmp(*x, p) >= 0)
			goto done;
		mpz_mul(t, *A, *A);
		mpz_sub_ui(t, t, 4);
		mpz_mod(t, t, p);
		mpz_gcd(t, t, p);
		if (mpz_cmp_ui(t, 1) != 0)
			goto done;
		if (mpz_cmp_ui(*o, 2) < 0)
			goto done;
		mpz_sqrt(L, p);                 /* L = isqrt(p) + 1 + isqrt(4 isqrt p) */
		mpz_mul_2exp(t, L, 2);
		mpz_sqrt(t, t);
		mpz_add(L, L, t);
		mpz_add_ui(L, L, 1);
		if (fast) {
			/* o < r L <= B L is implied by the window (r <= B), so rejecting
			 * here first changes no outcome and bounds the chunked residue */
			mpz_mul_ui(t, L, B);
			if (mpz_cmp(*o, t) >= 0)
				goto done;
			if (onepass)
				mpz_mod(g, Pr, *o);
			else
				primorial_mod(g, *o, words, nw);
		} else
			mpz_mod(g, Pr, *o);     /* rad(m) = gcd(P_B mod o, o) */
		mpz_gcd(g, g, *o);
		if (mpz_cmp_ui(g, 1) <= 0 || mpz_sizeinbase(g, 2) > radlim)
			goto done;
		int ns = factor_rad(div, g);
		nd = ns;
		mpz_set_ui(m, 1);               /* m = prod q^{v_q(o)} */
		for (int j = 0; j < ns; j++) {
			mpz_set(oo, *o);
			while (mpz_divisible_p(oo, div[j])) {
				mpz_divexact(oo, oo, div[j]);
				mpz_mul(m, m, div[j]);
			}
		}
		mpz_divexact(pn, *o, m);
		mpz_mul(t, div[0], L);
		if (mpz_cmp(L, *o) >= 0 || mpz_cmp(*o, t) >= 0)
			goto done;
		if (i + 3 == len) {             /* last level: p_next = 1 or < B^2 */
			if (mpz_cmp_ui(pn, 1) != 0 && mpz_cmp(pn, B2) >= 0)
				goto done;
		} else {                        /* B^2 < p_next, p_next^2 < p */
			mpz_mul(t, pn, pn);
			if (mpz_cmp(pn, B2) <= 0 || mpz_cmp(t, p) >= 0)
				goto done;
		}
		t_arith += now() - t0;
		t0 = now();
		mont_setup(p, *A);
		pt P0;
		to_mont(P0.X, *x);
		memcpy(P0.Z, ONE, N * sizeof(mp_limb_t));
		int lok;
		if (fast && mpz_cmp_ui(pn, 1) > 0)
			lok = level_fast(&P0, *o, div, ns, pn);
		else {
			if (mpz_cmp_ui(pn, 1) > 0)
				mpz_init_set(div[nd++], pn);
			lok = level_faithful(&P0, *o, div, nd);
		}
		t_ec += now() - t0;
		if (!lok)
			goto done;
		mpz_set(p, pn);
	}
	ok = 1;
done:
	for (int j = 0; j < nd; j++)
		mpz_clear(div[j]);
	free(words);
	mpz_clears(p, Pr, g, m, pn, oo, L, t, B2, NULL);
	return ok;
}

/* ------------------------------------------------------------------------- */
static void bench_chain(mpz_t *seq, int len)
{
	enum { MAXR = 2001 };
	static double tf[MAXR], tq[MAXR], tg[MAXR], tm[MAXR], fp[MAXR], fa[MAXR], fe[MAXR];
	int64_t mf, mq, mm;
	nmul = 0;
	if (!verify(seq, len, 0) || !verify(seq, len, 1)) {
		gmp_fprintf(stderr, "chain for %Zd REJECTED\n", seq[0]);
		exit(1);
	}
	nmul = 0;
	verify(seq, len, 0);
	mf = nmul;
	nmul = 0;
	t_prim = t_arith = t_ec = 0;
	verify(seq, len, 1);
	mq = nmul;
	double t1 = now();
	verify(seq, len, 0);
	int reps = (int)fmin(MAXR, fmax(21, 0.6 / fmax(now() - t1, 1e-7)));
	mpz_t a;
	mpz_init(a);
	mont_setup(seq[0], a);
	nmul = 0;
	mr_mont(seq[0], 2);
	mm = nmul;
	t_prim = t_arith = t_ec = 0;
	for (int r = 0; r < reps; r++) {
		double s = now();
		verify(seq, len, 0);
		double u = now();
		double p0 = t_prim, a0 = t_arith, e0 = t_ec;
		verify(seq, len, 1);
		double v = now();
		fp[r] = t_prim - p0;            /* fast-mode phases */
		fa[r] = t_arith - a0;
		fe[r] = t_ec - e0;
		mr_gmp(seq[0], 2 + r % 64);
		double w = now();
		mont_setup(seq[0], a);
		mr_mont(seq[0], 2 + r % 64);
		tf[r] = u - s;
		tq[r] = v - u;
		tg[r] = w - v;
		tm[r] = now() - w;
	}
	mpz_clear(a);
	double *arr[] = { tf, tq, tg, tm, fp, fa, fe };
	for (int i = 0; i < 7; i++)
		qsort(arr[i], reps, sizeof(double), cmpd);
	int h = reps / 2;
	/* phase split over both modes' runs (faithful: primorial; fast: sieve + words) */
	printf("%5zu bits, %2d levels:  faithful %8.3f ms (%6ld mults)   fast %8.3f ms (%6ld mults)\n",
	       mpz_sizeinbase(seq[0], 2), (len - 1) / 3, 1e3 * tf[h], (long)mf, 1e3 * tq[h], (long)mq);
	printf("      MR round: mpz_powm %7.3f ms, same-arith %7.3f ms (%ld mults)  "
	       "=> faithful %5.1f MR, fast %5.1f MR\n",
	       1e3 * tg[h], 1e3 * tm[h], (long)mm, tf[h] / tg[h], tq[h] / tg[h]);
	printf("      fast, in MR rounds: sieve + P_B mod prod(o_i) %.2f, level checks %.2f, curve arithmetic %.2f\n",
	       fp[h] / tg[h], fa[h] / tg[h], fe[h] / tg[h]);
}

static void bench(char *file, int argc, char **argv)
{
	FILE *f = fopen(file, "r");
	if (!f) {
		perror(file);
		exit(2);
	}
	static char line[1 << 20];
	static mpz_t seq[MAXSEQ];
	while (fgets(line, sizeof line, f)) {
		int len = 0;
		for (char *s = strtok(line, ", \t\r\n"); s && len < MAXSEQ; s = strtok(NULL, ", \t\r\n"))
			mpz_init_set_str(seq[len++], s, 10);
		if (len < 4)
			continue;
		size_t bits = mpz_sizeinbase(seq[0], 2);
		int want = argc == 0;
		for (int i = 0; i < argc; i++)
			want |= (size_t)atol(argv[i]) == bits;
		if (want)
			bench_chain(seq, len);
		for (int i = 0; i < len; i++)
			mpz_clear(seq[i]);
	}
	fclose(f);
}

int main(int argc, char **argv)
{
	xl_init();
	if (argc >= 3 && !strcmp(argv[1], "--bench")) {
		bench(argv[2], argc - 3, argv + 3);
		return 0;
	}
	int fast = argc >= 2 && !strcmp(argv[1], "--fast");
	static mpz_t seq[MAXSEQ];
	int c = argc - 1 - fast;
	if (c < 1 || c > MAXSEQ) {
		fprintf(stderr, "usage: vshort2 [--fast] p0 A0 x0 o0 [A1 x1 o1 ...]  |  vshort2 --bench certs.csv [bits ...]\n");
		return 2;
	}
	for (int i = 0; i < c; i++)
		if (mpz_init_set_str(seq[i], argv[1 + fast + i], 10) != 0) {
			fprintf(stderr, "usage: vshort2 [--fast] p0 A0 x0 o0 [A1 x1 o1 ...]\n");
			return 2;
		}
	int ok = verify(seq, c, fast);
	puts(ok ? "True" : "False");
	return ok ? 0 : 1;
}
