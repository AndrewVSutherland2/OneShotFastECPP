#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <gmp.h>
#include "cornacchia.h"

/*
    Factor-base discriminant scan (fastECPP substeps 1+2, residue-factor variant).

      dscan  p=<decimal> | pbits=<n> [seed=<n>]   B=<bound>   [dump] [verify] [nocorn]

    Enumerate every negative fundamental discriminant D with |D| < B all of whose
    prime-discriminant factors q-star are quadratic residues mod p, i.e. the
    Kronecker symbol (q-star | p) = 1, and report the traces t for which
    4p = t^2 - v^2 D = t^2 + |D| v^2 is solvable.

    Each D is a product of distinct prime discriminants from a factor base
    Q = { q-star : (q-star | p) = 1 }, where q-star = q (q==1 mod4), -q (q==3 mod4),
    or one of -4,+8,-8.  We precompute sqrt(q-star) mod p for each (the only modular
    square roots), and a DFS over products maintains sqrt(D mod p) = product of the
    sqrt(q-star) incrementally, so Cornacchia runs with a supplied root and never
    recomputes a per-D Tonelli.
*/

#define MAXDEPTH 48

static double wall (void)
{ struct timespec ts; clock_gettime (CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }

// Legendre symbol (a/q), q an odd prime not dividing a, via the binary algorithm.
static int ui_legendre (unsigned long a, unsigned long q)
{
    int k = 1;
    a %= q;
    while ( a ) {
        while ( !(a & 1) ) { a >>= 1; if ( (q & 7) == 3 || (q & 7) == 5 ) k = -k; }
        unsigned long t = a; a = q; q = t;                  // swap
        if ( (a & 3) == 3 && (q & 3) == 3 ) k = -k;
        a %= q;
    }
    return q == 1 ? k : 0;
}

// ---- factor base ----
typedef struct {
    long  *sval;        // signed prime discriminant q*
    int   *sgn;         // sign of q* (+1/-1)
    unsigned long *aval;// |q*|
    mpz_t *sq;          // sqrt(q* mod p)
    size_t n, cap;
} fbase;

// reserve room for n factor-base entries
static void fb_reserve (fbase *fb, size_t n)
{
    fb->cap = n;
    fb->sval = malloc (n * sizeof(*fb->sval));
    fb->sgn  = malloc (n * sizeof(*fb->sgn));
    fb->aval = malloc (n * sizeof(*fb->aval));
    fb->sq   = malloc (n * sizeof(*fb->sq));
}

// append a prime discriminant q*; its square root is filled in later, in parallel
static void fb_push (fbase *fb, long sval)
{
    size_t i = fb->n++;
    fb->sval[i] = sval;
    fb->sgn[i]  = sval < 0 ? -1 : 1;
    fb->aval[i] = sval < 0 ? (unsigned long)(-sval) : (unsigned long) sval;
    mpz_init (fb->sq[i]);
}

// ---- DFS state ----
static fbase oddfb;                 // odd prime discriminants with (q*/p)=1, ascending |q*|
static long  tp_val[4]; static mpz_t tp_sq[4]; static int n_tp;  // 2-part options incl. "1"
static unsigned long Bbound;
static unsigned long Bminb;         // emit only d >= Bminb (incremental rescans)
static int do_dump, do_verify, no_corn, do_dumpscan;

// Per-thread scan state.  The DFS partitions cleanly by the smallest prime factor
// of D (each subtree only extends by larger primes), so the top-level loop over the
// starting prime is data-parallel: each thread runs whole subtrees with its own
// context, sqrt-product stack, scratch, counters, and collected output.
typedef struct {
    cornacchia_ctx cc;
    mpz_t xstk[MAXDEPTH], x0, t, v, chk, dv2;   // xstk[k] = sqrt(odd-part mod p) at depth k
    unsigned long n_scan, n_solv;
    char *obuf; size_t olen, ocap;              // dump/dumpscan output, merged at the end
} scanstate;

static void obuf_add (scanstate *S, const char *s, int n)
{
    if ( S->olen + n + 1 > S->ocap ) {
        if ( ! S->ocap ) S->ocap = 8192;
        while ( S->olen + n + 1 > S->ocap ) S->ocap *= 2;
        S->obuf = realloc (S->obuf, S->ocap);
    }
    memcpy (S->obuf + S->olen, s, n);  S->olen += n;
}

static void try_D (scanstate *S, unsigned long absprod, int sign, const mpz_t xprod)
{
    for ( int j = 0 ; j < n_tp ; j++ ) {
        unsigned long tabs = tp_val[j] < 0 ? (unsigned long)(-tp_val[j]) : (unsigned long) tp_val[j];
        unsigned long dabs = absprod * tabs;
        if ( dabs >= Bbound || dabs < 3 ) continue;
        int dsign = sign * (tp_val[j] < 0 ? -1 : 1);
        if ( dsign >= 0 ) continue;                          // need D = -dabs < 0
        S->n_scan++;
        if ( dabs < Bminb ) continue;                        // below the incremental lower bound
        if ( do_dumpscan ) { char b[24]; int n = snprintf (b, sizeof b, "%lu\n", dabs); obuf_add (S, b, n); }
        if ( no_corn ) continue;
        mpz_mul (S->x0, xprod, tp_sq[j]);  mpz_mod (S->x0, S->x0, S->cc.p);   // sqrt(D mod p)
        if ( cornacchia_solve_x0 (&S->cc, dabs, S->x0, S->t, S->v) ) {
            // self-check 4p = t^2 + dabs*v^2
            mpz_mul (S->chk, S->t, S->t);  mpz_mul (S->dv2, S->v, S->v);  mpz_mul_ui (S->dv2, S->dv2, dabs);
            mpz_add (S->chk, S->chk, S->dv2);
            if ( mpz_cmp (S->chk, S->cc.p4) != 0 ) { gmp_fprintf (stderr, "SELF-CHECK FAIL D=-%lu\n", dabs); exit (2); }
            if ( do_verify ) {                               // cross-check vs full Tonelli solve
                mpz_t t2, v2; mpz_init (t2); mpz_init (v2);
                if ( ! cornacchia_solve (&S->cc, dabs, t2, v2) || mpz_cmp (t2, S->t) != 0 )
                    { gmp_fprintf (stderr, "VERIFY MISMATCH D=-%lu\n", dabs); exit (3); }
                mpz_clear (t2); mpz_clear (v2);
            }
            S->n_solv++;
            if ( do_dump ) { char *b; int n = gmp_asprintf (&b, "%lu %Zd %Zd\n", dabs, S->t, S->v); obuf_add (S, b, n); free (b); }
        }
    }
}

// serial DFS over a subtree: process the D's at this node then recurse into children.
static void dfs_serial (scanstate *S, int depth, size_t start, unsigned long absprod, int sign, const mpz_t xprod)
{
    try_D (S, absprod, sign, xprod);
    for ( size_t i = start ; i < oddfb.n ; i++ ) {
        if ( absprod > (Bbound - 1) / oddfb.aval[i] ) break;
        mpz_mul (S->xstk[depth + 1], xprod, oddfb.sq[i]);  mpz_mod (S->xstk[depth + 1], S->xstk[depth + 1], S->cc.p);
        dfs_serial (S, depth + 1, i + 1, absprod * oddfb.aval[i], sign * oddfb.sgn[i], S->xstk[depth + 1]);
    }
}

// A "seed" is an independent unit of parallel work: children [ilo,ihi) of a node
// (absprod,sign,xprod), each of whose subtrees is small (<~ SEED_D discriminants),
// batched so the seed is ~SEED_D discriminants total.  Processing a seed = a serial
// DFS of each of those children -- fully independent, no shared state, no suspension.
typedef struct { unsigned long absprod; int sign; size_t ilo, ihi; mpz_t xprod; } seed_t;
static seed_t *seeds; static size_t nseeds, seeds_cap;
static unsigned long SEED_D = 2000;

static void emit_seed (unsigned long absprod, int sign, const mpz_t xprod, size_t ilo, size_t ihi)
{
    if ( nseeds == seeds_cap ) { seeds_cap = seeds_cap ? 2 * seeds_cap : 4096; seeds = realloc (seeds, seeds_cap * sizeof(*seeds)); }
    seed_t *s = &seeds[nseeds++];
    s->absprod = absprod; s->sign = sign; s->ilo = ilo; s->ihi = ihi;
    mpz_init (s->xprod); mpz_set (s->xprod, xprod);
}

// Serial pass: process the D's at "upper" nodes (few) and emit seeds covering all
// the small subtrees (the bulk).  Recurses only into big subtrees (est. size >= SEED_D).
static void gen_seeds (scanstate *S0, int depth, size_t start, unsigned long absprod, int sign, const mpz_t xprod)
{
    try_D (S0, absprod, sign, xprod);
    size_t i = start;
    while ( i < oddfb.n ) {                                   // big children: recurse to break them up
        if ( absprod > (Bbound - 1) / oddfb.aval[i] ) return;
        unsigned long ca = absprod * oddfb.aval[i];
        if ( Bbound / ca < SEED_D ) break;                   // this and all later children are small
        mpz_mul (S0->xstk[depth + 1], xprod, oddfb.sq[i]);  mpz_mod (S0->xstk[depth + 1], S0->xstk[depth + 1], S0->cc.p);
        gen_seeds (S0, depth + 1, i + 1, ca, sign * oddfb.sgn[i], S0->xstk[depth + 1]);
        i++;
    }
    while ( i < oddfb.n && absprod <= (Bbound - 1) / oddfb.aval[i] ) {   // small children: batch into seeds
        size_t j = i; unsigned long batch = 0;
        while ( j < oddfb.n && absprod <= (Bbound - 1) / oddfb.aval[j] && batch < SEED_D )
            { batch += Bbound / (absprod * oddfb.aval[j]); j++; }
        emit_seed (absprod, sign, xprod, i, j);
        i = j;
    }
}

int main (int argc, char *argv[])
{
    mpz_t p;  mpz_init (p);  mpz_set_ui (p, 0);
    unsigned long pbits = 0, seed = 1;
    int n_threads = 0;
    Bbound = 1000000;
    for ( int i = 1 ; i < argc ; i++ ) {
        if ( ! strncmp (argv[i], "p=", 2) ) mpz_set_str (p, argv[i] + 2, 10);
        else if ( ! strncmp (argv[i], "pbits=", 6) ) pbits = strtoul (argv[i] + 6, 0, 10);
        else if ( ! strncmp (argv[i], "B=", 2) ) Bbound = strtoul (argv[i] + 2, 0, 10);
        else if ( ! strncmp (argv[i], "Bmin=", 5) ) Bminb = strtoul (argv[i] + 5, 0, 10);
        else if ( ! strncmp (argv[i], "seed=", 5) ) seed = strtoul (argv[i] + 5, 0, 10);
        else if ( ! strncmp (argv[i], "threads=", 8) ) n_threads = atoi (argv[i] + 8);
        else if ( ! strncmp (argv[i], "sd=", 3) ) SEED_D = strtoul (argv[i] + 3, 0, 10);
        else if ( ! strcmp (argv[i], "dump") ) do_dump = 1;
        else if ( ! strcmp (argv[i], "verify") ) do_verify = 1;
        else if ( ! strcmp (argv[i], "nocorn") ) no_corn = 1;
        else if ( ! strcmp (argv[i], "dumpscan") ) { do_dumpscan = 1; no_corn = 1; }
    }
    if ( ! mpz_sgn (p) ) {
        if ( ! pbits ) pbits = 256;
        gmp_randstate_t rs; gmp_randinit_default (rs); gmp_randseed_ui (rs, seed);
        mpz_urandomb (p, rs, pbits);  mpz_setbit (p, pbits - 1);  mpz_nextprime (p, p);
        gmp_randclear (rs);
    }
    if ( ! mpz_probab_prime_p (p, 25) ) { fprintf (stderr, "p not prime\n"); return 1; }

    cornacchia_ctx ctx;  cornacchia_init (&ctx, p);        // for the 2-part sqrts and p mod 8

    // ---- build factor base ----
    int pmod8 = ctx.pmod8;
    if ( n_threads > 0 ) omp_set_num_threads (n_threads);
    int nth = omp_get_max_threads ();

    // 2-part prime discriminants -4,8,-8 (serial; only three)
    n_tp = 0;  tp_val[n_tp] = 1; mpz_init_set_ui (tp_sq[n_tp], 1); n_tp++;   // the "no 2-part" option
    long tps[3] = { -4, 8, -8 };
    for ( int k = 0 ; k < 3 ; k++ ) {
        int qr = 0; long e = tps[k];
        if ( e == -4 ) qr = (pmod8 == 1 || pmod8 == 5);
        else if ( e == 8 ) qr = (pmod8 == 1 || pmod8 == 7);
        else /* -8 */ qr = (pmod8 == 1 || pmod8 == 3);
        if ( qr ) { tp_val[n_tp] = e; mpz_init (tp_sq[n_tp]); mpz_t a; mpz_init_set_si (a, e); mpz_mod (a, a, p);
                    cornacchia_sqrtmodp (&ctx, tp_sq[n_tp], a); mpz_clear (a); n_tp++; }
    }

    // sieve odd primes < B
    double t_sieve = wall ();
    unsigned long N = Bbound;
    size_t nb = (N + 63) / 64;
    uint64_t *bits = calloc (nb, sizeof(uint64_t));          // bit q set => q composite
    for ( unsigned long q = 3 ; q * q < N ; q += 2 )
        if ( ! (bits[q >> 6] & (1ULL << (q & 63))) )
            for ( unsigned long m = q * q ; m < N ; m += 2 * q ) bits[m >> 6] |= 1ULL << (m & 63);
    size_t pcap = (size_t)(N / (log ((double) N) - 1.1)) + 1024;
    uint32_t *primes = malloc (pcap * sizeof(uint32_t));
    size_t np = 0;
    for ( unsigned long q = 3 ; q < N ; q += 2 )
        if ( ! (bits[q >> 6] & (1ULL << (q & 63))) ) primes[np++] = (uint32_t) q;
    free (bits);
    t_sieve = wall () - t_sieve;

    // PASS 1 (parallel, Legendre only): which primes are QR mod p
    double t_det = wall ();
    uint8_t *qr = malloc (np);
    #pragma omp parallel for schedule(static)
    for ( size_t i = 0 ; i < np ; i++ ) {
        unsigned long q = primes[i];
        unsigned long pmq = mpz_fdiv_ui (p, q);              // read-only p -> thread-safe
        qr[i] = (uint8_t)(pmq && ui_legendre (pmq, q) == 1);
    }
    t_det = wall () - t_det;

    // collect QR prime discriminants, ascending
    size_t n_qr = 0; for ( size_t i = 0 ; i < np ; i++ ) n_qr += qr[i];
    fb_reserve (&oddfb, n_qr);
    for ( size_t i = 0 ; i < np ; i++ )
        if ( qr[i] ) { unsigned long q = primes[i]; fb_push (&oddfb, (q & 3) == 1 ? (long) q : -(long) q); }

    // PASS 2 (parallel, Tonelli): square roots of the factor-base elements
    double t_sqrt = wall ();
    cornacchia_ctx *tctx = malloc (nth * sizeof(*tctx));
    for ( int k = 0 ; k < nth ; k++ ) cornacchia_init (&tctx[k], p);
    #pragma omp parallel
    {
        cornacchia_ctx *lc = &tctx[omp_get_thread_num ()];
        mpz_t a; mpz_init (a);
        #pragma omp for schedule(dynamic, 256)
        for ( size_t i = 0 ; i < oddfb.n ; i++ ) {
            mpz_set_si (a, oddfb.sval[i]);  mpz_mod (a, a, p);
            cornacchia_sqrtmodp (lc, oddfb.sq[i], a);
        }
        mpz_clear (a);
    }
    t_sqrt = wall () - t_sqrt;
    for ( int k = 0 ; k < nth ; k++ ) cornacchia_clear (&tctx[k]);
    free (tctx); free (primes); free (qr);
    double tfb = t_sieve + t_det + t_sqrt;

    // ---- DFS scan (parallel over the starting prime, i.e. smallest factor of D) ----
    double tsc = wall ();
    scanstate *states = malloc (nth * sizeof(*states));
    for ( int k = 0 ; k < nth ; k++ ) {
        scanstate *S = &states[k];  memset (S, 0, sizeof(*S));
        cornacchia_init (&S->cc, p);
        for ( int d2 = 0 ; d2 < MAXDEPTH ; d2++ ) mpz_init (S->xstk[d2]);
        mpz_set_ui (S->xstk[0], 1);
        mpz_init (S->x0); mpz_init (S->t); mpz_init (S->v); mpz_init (S->chk); mpz_init (S->dv2);
    }
    // serial pass: emit balanced seeds (and process the few "upper" nodes' D's on states[0])
    double t_gen = wall ();
    mpz_t ONE; mpz_init_set_ui (ONE, 1);
    gen_seeds (&states[0], 0, 0, 1, 1, ONE);
    mpz_clear (ONE);
    t_gen = wall () - t_gen;
    // parallel pass: each seed is an independent chunk of ~SEED_D discriminants
    double t_par = wall ();
    #pragma omp parallel for schedule(dynamic, 4)
    for ( size_t s = 0 ; s < nseeds ; s++ ) {
        scanstate *S = &states[omp_get_thread_num ()];
        seed_t *sd = &seeds[s];
        for ( size_t k = sd->ilo ; k < sd->ihi ; k++ ) {
            mpz_mul (S->xstk[1], sd->xprod, oddfb.sq[k]);  mpz_mod (S->xstk[1], S->xstk[1], S->cc.p);
            dfs_serial (S, 1, k + 1, sd->absprod * oddfb.aval[k], sd->sign * oddfb.sgn[k], S->xstk[1]);
        }
    }
    t_par = wall () - t_par;
    for ( size_t s = 0 ; s < nseeds ; s++ ) mpz_clear (seeds[s].xprod);
    free (seeds);
    unsigned long n_scan = 0, n_solv = 0, n_desc = 0;
    for ( int k = 0 ; k < nth ; k++ ) { n_scan += states[k].n_scan; n_solv += states[k].n_solv; n_desc += states[k].cc.n_descent; }
    tsc = wall () - tsc;

    for ( int k = 0 ; k < nth ; k++ ) if ( states[k].olen ) fwrite (states[k].obuf, 1, states[k].olen, stdout);
    for ( int k = 0 ; k < nth ; k++ ) {
        scanstate *S = &states[k];
        cornacchia_clear (&S->cc);
        for ( int d2 = 0 ; d2 < MAXDEPTH ; d2++ ) mpz_clear (S->xstk[d2]);
        mpz_clear (S->x0); mpz_clear (S->t); mpz_clear (S->v); mpz_clear (S->chk); mpz_clear (S->dv2);
        free (S->obuf);
    }
    free (states);

    if ( ! do_dump && ! do_dumpscan ) {
        gmp_fprintf (stderr, "p (%lu bits, %lu mod 4) = %Zd\n", (unsigned long) mpz_sizeinbase (p, 2), mpz_fdiv_ui (p, 4), p);
        fprintf (stderr, "B = %lu   factor base: %zu odd prime discriminants + %d two-part\n", Bbound, oddfb.n, n_tp - 1);
        fprintf (stderr, "factor-base build: %.3f s  [sieve %.3f + determine %.3f (Legendre) + %zu roots %.3f (Tonelli)]\n",
                 tfb, t_sieve, t_det, oddfb.n, t_sqrt);
        fprintf (stderr, "scan: %.3f s   D scanned: %lu (%.3f us/D, %.0f D/s, %d threads)\n",
                 tsc, n_scan, n_scan ? 1e6 * tsc / n_scan : 0, n_scan ? n_scan / tsc : 0, nth);
        fprintf (stderr, "  seeds: %zu (SEED_D=%lu)   gen_seeds(serial): %.3f s   parallel: %.3f s\n",
                 nseeds, SEED_D, t_gen, t_par);
        fprintf (stderr, "solvable: %lu (1 per %.1f D)   descents: %lu\n",
                 n_solv, n_solv ? (double) n_scan / n_solv : 0, n_desc);
    }
    return 0;
}
