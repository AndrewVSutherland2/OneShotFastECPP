#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <gmp.h>
#include "smooth.h"

/*
    Driver / test harness for the batched n^4-smoothness engine (component b).

    Modes:
      smoothtest pbuild  y=<Y> [threads=N] [save=<file>]
          Build P = prod_{q<=Y} q (optionally cache it to <file>).

      smoothtest parts   (y=<Y> | load=<file>) [threads=N]
          Read decimal integers from stdin (whitespace-separated); print
          "N  smoothpart  cofactor" for each.  For unit testing vs PARI.

      smoothtest gate    (p=<dec> | pbits=<n> [seed=<s>]) [y=<Y> | load=<file>]
                         [threads=N]
          Read solvable discriminant records "d t v" from stdin (as produced by
          `dscan ... dump`), form the two candidate curve orders N = p+1 -/+ t,
          batch-test their n^4-smooth parts, and report every candidate whose
          smooth part exceeds L (a valid one-shot m exists).  Each winner is
          self-checked (m|N, L<m<Hass, m<L*r, m is n^4-smooth).
*/

static double wall (void)
{ struct timespec ts; clock_gettime (CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }

// regenerate the same p as dscan for given pbits/seed
static void gen_p (mpz_t p, unsigned long pbits, unsigned long seed)
{
    gmp_randstate_t rs;  gmp_randinit_default (rs);  gmp_randseed_ui (rs, seed);
    mpz_urandomb (p, rs, pbits);  mpz_setbit (p, pbits - 1);  mpz_nextprime (p, p);
    gmp_randclear (rs);
}

static int get_base (smooth_base *sb, uint64_t y, const char *loadpath, int nth, int verbose)
{
    double t0 = wall ();
    if ( loadpath ) {
        if ( ! smooth_base_load (sb, loadpath) ) { fprintf (stderr, "load %s failed\n", loadpath); return 0; }
        if ( verbose ) fprintf (stderr, "loaded P: y=%lu, %lu primes, %zu MiB, %.2fs\n",
                                (unsigned long)sb->y, (unsigned long)sb->nprimes,
                                (mpz_sizeinbase(sb->P,2)/8)>>20, wall()-t0);
    } else {
        smooth_base_build (sb, y, nth);
        if ( verbose ) fprintf (stderr, "built P: y=%lu, %lu primes, %.0f Mbit (%zu MiB), %.2fs\n",
                                (unsigned long)sb->y, (unsigned long)sb->nprimes,
                                mpz_sizeinbase(sb->P,2)/1e6, (mpz_sizeinbase(sb->P,2)/8)>>20, wall()-t0);
    }
    if ( ! smooth_base_selfcheck (sb) ) {
        fprintf (stderr, "ERROR: P failed self-check (corrupt cache or build?)\n");
        return 0;
    }
    return 1;
}

// ---- winner self-check (independent of build_m) ----
static int check_winner (const mpz_t p, const mpz_t N, const mpz_t m,
                         const mpz_t L, const mpz_t Hass, uint64_t n2, uint64_t n4,
                         const uint64_t *qs, int nq)
{
    (void) p;
    if ( ! mpz_divisible_p (N, m) ) { fprintf (stderr, "  CHECK m does not divide N\n"); return 0; }
    if ( mpz_cmp (m, L) <= 0 )      { fprintf (stderr, "  CHECK m<=L\n"); return 0; }
    if ( mpz_cmp (m, Hass) > 0 )    { fprintf (stderr, "  CHECK m>Hasse\n"); return 0; }
    // factor m; confirm n^4-smooth and least prime r gives m < L*r
    uint64_t pr[64];  int ex[64];
    int k = factor_smooth (m, pr, ex, 64);
    if ( k == 0 ) return 0;
    if ( pr[k-1] > n4 ) { fprintf (stderr, "  CHECK m not n^4-smooth (max prime %lu)\n", (unsigned long)pr[k-1]); return 0; }
    mpz_t Lr;  mpz_init (Lr);  mpz_mul_ui (Lr, L, pr[0]);
    int lt = mpz_cmp (m, Lr) < 0;  mpz_clear (Lr);
    if ( ! lt ) { fprintf (stderr, "  CHECK m>=L*r\n"); return 0; }
    // q_i must be exactly the distinct primes of m in (n2,n4)
    int want = 0;
    for ( int j = 0 ; j < k ; j++ ) if ( pr[j] > n2 && pr[j] < n4 ) {
        if ( want >= nq || qs[want] != pr[j] ) { fprintf (stderr, "  CHECK q_i mismatch\n"); return 0; }
        want++;
    }
    if ( want != nq ) { fprintf (stderr, "  CHECK extra q_i\n"); return 0; }
    return 1;
}

static int mode_parts (int argc, char **argv)
{
    uint64_t y = 0;  const char *load = NULL;  int nth = 0;
    for ( int i = 2 ; i < argc ; i++ ) {
        if ( ! strncmp (argv[i], "y=", 2) ) y = strtoull (argv[i]+2, 0, 10);
        else if ( ! strncmp (argv[i], "load=", 5) ) load = argv[i]+5;
        else if ( ! strncmp (argv[i], "threads=", 8) ) nth = atoi (argv[i]+8);
    }
    if ( ! y && ! load ) { fprintf (stderr, "parts: need y= or load=\n"); return 2; }
    smooth_base sb;
    if ( ! get_base (&sb, y, load, nth, 1) ) return 1;

    // read all integers
    size_t cap = 1024, n = 0;
    mpz_t *N = malloc (cap * sizeof(mpz_t));
    for (;;) {
        if ( n == cap ) { cap *= 2; N = realloc (N, cap*sizeof(mpz_t)); }
        mpz_init (N[n]);
        if ( mpz_inp_str (N[n], stdin, 10) == 0 ) { mpz_clear (N[n]); break; }
        n++;
    }
    mpz_t *S = malloc (n * sizeof(mpz_t));
    for ( size_t i = 0 ; i < n ; i++ ) mpz_init (S[i]);
    double t0 = wall ();
    smooth_parts (&sb, (const mpz_t*)N, n, S, nth);
    double dt = wall () - t0;
    mpz_t cof;  mpz_init (cof);
    for ( size_t i = 0 ; i < n ; i++ ) {
        mpz_divexact (cof, N[i], S[i]);
        gmp_printf ("%Zd %Zd %Zd\n", N[i], S[i], cof);
    }
    fprintf (stderr, "%zu numbers, smooth_parts %.3fs (%.1f us/num)\n", n, dt, n?1e6*dt/n:0);
    mpz_clear (cof);
    for ( size_t i = 0 ; i < n ; i++ ) { mpz_clear (N[i]); mpz_clear (S[i]); }
    free (N); free (S);
    smooth_base_clear (&sb);
    return 0;
}

static int mode_gate (int argc, char **argv)
{
    mpz_t p;  mpz_init (p);
    unsigned long pbits = 0, seed = 1;  uint64_t y = 0;  const char *load = NULL;  int nth = 0;
    for ( int i = 2 ; i < argc ; i++ ) {
        if ( ! strncmp (argv[i], "p=", 2) ) mpz_set_str (p, argv[i]+2, 10);
        else if ( ! strncmp (argv[i], "pbits=", 6) ) pbits = strtoul (argv[i]+6, 0, 10);
        else if ( ! strncmp (argv[i], "seed=", 5) ) seed = strtoul (argv[i]+5, 0, 10);
        else if ( ! strncmp (argv[i], "y=", 2) ) y = strtoull (argv[i]+2, 0, 10);
        else if ( ! strncmp (argv[i], "load=", 5) ) load = argv[i]+5;
        else if ( ! strncmp (argv[i], "threads=", 8) ) nth = atoi (argv[i]+8);
    }
    if ( ! mpz_sgn (p) ) { if ( ! pbits ) pbits = 256;  gen_p (p, pbits, seed); }
    if ( ! mpz_probab_prime_p (p, 25) ) { fprintf (stderr, "p not prime\n"); return 1; }

    mpz_t L, Hass;  mpz_init (L);  mpz_init (Hass);
    unsigned long n;  uint64_t n2, n4;
    cert_bounds (p, L, Hass, &n, &n2, &n4);
    if ( ! y ) y = n4;
    gmp_fprintf (stderr, "p = %Zd\n", p);
    fprintf (stderr, "n=%lu  n^2=%lu  n^4=%lu  (smoothness bound y=%lu)\n",
             n, (unsigned long)n2, (unsigned long)n4, (unsigned long)y);
    gmp_fprintf (stderr, "L = %Zd  (need smooth part > L)\n", L);

    smooth_base sb;
    if ( ! get_base (&sb, y, load, nth, 1) ) return 1;
    if ( sb.y != y ) fprintf (stderr, "WARNING: P built for y=%lu, wanted %lu\n",
                              (unsigned long)sb.y, (unsigned long)y);

    // read "d t v"; make candidates N = p+1 -/+ t
    size_t cap = 4096, nc = 0;
    mpz_t *N = malloc (cap * sizeof(mpz_t));
    long   *sgn = malloc (cap * sizeof(long));             // -1 => p+1-t, +1 => p+1+t
    unsigned long *dd = malloc (cap * sizeof(unsigned long));
    mpz_t *tt = malloc (cap * sizeof(mpz_t));
    mpz_t t, v, tmp;  mpz_init (t);  mpz_init (v);  mpz_init (tmp);
    unsigned long d;
    while ( gmp_fscanf (stdin, "%lu %Zd %Zd", &d, t, v) == 3 ) {
        for ( int s = 0 ; s < 2 ; s++ ) {
            if ( nc == cap ) {
                cap *= 2;
                N = realloc (N, cap*sizeof(mpz_t));  sgn = realloc (sgn, cap*sizeof(long));
                dd = realloc (dd, cap*sizeof(unsigned long));  tt = realloc (tt, cap*sizeof(mpz_t));
            }
            mpz_init (N[nc]);  mpz_add_ui (tmp, p, 1);
            if ( s == 0 ) mpz_sub (N[nc], tmp, t); else mpz_add (N[nc], tmp, t);
            sgn[nc] = s ? +1 : -1;  dd[nc] = d;  mpz_init_set (tt[nc], t);
            nc++;
        }
    }
    fprintf (stderr, "candidates: %zu  (from solvable D on stdin)\n", nc);

    mpz_t *S = malloc (nc * sizeof(mpz_t));
    for ( size_t i = 0 ; i < nc ; i++ ) mpz_init (S[i]);
    double t0 = wall ();
    smooth_parts (&sb, (const mpz_t*)N, nc, S, nth);
    double dt = wall () - t0;

    // gate + winners
    mpz_t m;  mpz_init (m);
    uint64_t qs[64];  int nq;
    size_t nwin = 0;
    for ( size_t i = 0 ; i < nc ; i++ ) {
        if ( mpz_cmp (S[i], L) <= 0 ) continue;
        // winner: build m and self-check
        if ( ! build_m (m, qs, &nq, S[i], L, n2, n4) ) {
            gmp_fprintf (stderr, "  build_m failed (edge) for D=-%lu sgn=%+ld\n", dd[i], sgn[i]);
            continue;
        }
        int ok = check_winner (p, N[i], m, L, Hass, n2, n4, qs, nq);
        nwin++;
        gmp_printf ("WIN D=-%lu t=%Zd order=%Zd  m=%Zd  smoothpart=%Zd  q=[",
                    dd[i], tt[i], N[i], m, S[i]);
        for ( int j = 0 ; j < nq ; j++ ) printf ("%s%lu", j?",":"", (unsigned long)qs[j]);
        printf ("]  %s\n", ok ? "OK" : "SELFCHECK-FAIL");
        if ( ! ok ) { fprintf (stderr, "SELF-CHECK FAILED, aborting\n"); return 3; }
    }
    fprintf (stderr, "smooth_parts: %.3fs (%.1f us/candidate, %d threads)\n",
             dt, nc?1e6*dt/nc:0, nth>0?nth:omp_get_max_threads());
    fprintf (stderr, "winners: %zu / %zu candidates  (yield %.3g)\n",
             nwin, nc, nc?(double)nwin/nc:0);

    mpz_clear (m);  mpz_clear (t);  mpz_clear (v);  mpz_clear (tmp);
    for ( size_t i = 0 ; i < nc ; i++ ) { mpz_clear (N[i]); mpz_clear (S[i]); mpz_clear (tt[i]); }
    free (N); free (S); free (sgn); free (dd); free (tt);
    mpz_clear (L);  mpz_clear (Hass);  mpz_clear (p);
    smooth_base_clear (&sb);
    return 0;
}

static int mode_pbuild (int argc, char **argv)
{
    uint64_t y = 0;  int nth = 0;  const char *save = NULL;
    for ( int i = 2 ; i < argc ; i++ ) {
        if ( ! strncmp (argv[i], "y=", 2) ) y = strtoull (argv[i]+2, 0, 10);
        else if ( ! strncmp (argv[i], "threads=", 8) ) nth = atoi (argv[i]+8);
        else if ( ! strncmp (argv[i], "save=", 5) ) save = argv[i]+5;
    }
    if ( ! y ) { fprintf (stderr, "pbuild: need y=\n"); return 2; }
    smooth_base sb;
    double t0 = wall ();
    smooth_base_build (&sb, y, nth);
    fprintf (stderr, "built P: y=%lu, %lu primes, %.0f Mbit (%zu MiB), %.2fs\n",
             (unsigned long)y, (unsigned long)sb.nprimes, mpz_sizeinbase(sb.P,2)/1e6,
             (mpz_sizeinbase(sb.P,2)/8)>>20, wall()-t0);
    if ( save ) {
        if ( ! smooth_base_save (&sb, save) ) { fprintf (stderr, "save failed\n"); return 1; }
        fprintf (stderr, "saved to %s\n", save);
    }
    smooth_base_clear (&sb);
    return 0;
}

int main (int argc, char **argv)
{
    if ( argc < 2 ) {
        fprintf (stderr,
          "usage:\n"
          "  smoothtest pbuild y=<Y> [threads=N] [save=<file>]\n"
          "  smoothtest parts  (y=<Y>|load=<file>) [threads=N]   < integers\n"
          "  smoothtest gate   (p=<dec>|pbits=<n> [seed=<s>]) [y=<Y>|load=<file>] [threads=N]  < 'd t v' lines\n");
        return 2;
    }
    if ( ! strcmp (argv[1], "pbuild") ) return mode_pbuild (argc, argv);
    if ( ! strcmp (argv[1], "parts") )  return mode_parts (argc, argv);
    if ( ! strcmp (argv[1], "gate") )   return mode_gate (argc, argv);
    fprintf (stderr, "unknown mode '%s'\n", argv[1]);
    return 2;
}
