#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include "cornacchia.h"

/*
    Driver for the baseline Cornacchia CM-discriminant search.

      cmsearch p=<decimal> | pbits=<n> [seed=<n>]  maxd=<N>  [dump] [verbose]

    Iterates fundamental discriminants D = -d with d <= maxd, solving
    4p = t^2 + d v^2.  Default: benchmark (throughput + counts, self-checking
    every solution).  `dump` prints "d solvable [t v]" for cross-checking against
    PARI.  `verbose` lists the solutions found.
*/

static double wall (void)
{
    struct timespec ts;  clock_gettime (CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

// Generate fundamental |D| values d in (0, maxd]: d == 3 mod 4 squarefree, or
// d == 0 mod 4 with d/4 == 1,2 mod 4 squarefree.  Returns count, fills *out (malloc'd).
static size_t gen_fundamental (unsigned long maxd, unsigned long **out)
{
    char *sf = malloc (maxd + 1);
    memset (sf, 1, maxd + 1);
    for ( unsigned long i = 2 ; i * i <= maxd ; i++ )
        for ( unsigned long j = i * i ; j <= maxd ; j += i * i ) sf[j] = 0;
    unsigned long *v = malloc ((maxd / 2 + 2) * sizeof(*v));
    size_t n = 0;
    for ( unsigned long d = 3 ; d <= maxd ; d++ ) {
        int ok = 0;
        if ( (d & 3) == 3 ) ok = sf[d];
        else if ( (d & 3) == 0 ) { unsigned long m = d >> 2; ok = ((m & 3) == 1 || (m & 3) == 2) && sf[m]; }
        if ( ok ) v[n++] = d;
    }
    free (sf);
    *out = v;
    return n;
}

int main (int argc, char *argv[])
{
    mpz_t p, t, v, chk;
    unsigned long maxd = 1000000, pbits = 0, seed = 1;
    int dump = 0, verbose = 0;
    mpz_init (p);  mpz_set_ui (p, 0);

    for ( int i = 1 ; i < argc ; i++ ) {
        if ( ! strncmp (argv[i], "p=", 2) ) mpz_set_str (p, argv[i] + 2, 10);
        else if ( ! strncmp (argv[i], "pbits=", 6) ) pbits = strtoul (argv[i] + 6, 0, 10);
        else if ( ! strncmp (argv[i], "maxd=", 5) ) maxd = strtoul (argv[i] + 5, 0, 10);
        else if ( ! strncmp (argv[i], "seed=", 5) ) seed = strtoul (argv[i] + 5, 0, 10);
        else if ( ! strcmp (argv[i], "dump") ) dump = 1;
        else if ( ! strcmp (argv[i], "verbose") ) verbose = 1;
    }
    if ( ! mpz_sgn (p) ) {
        if ( ! pbits ) pbits = 256;
        gmp_randstate_t rs;  gmp_randinit_default (rs);  gmp_randseed_ui (rs, seed);
        mpz_urandomb (p, rs, pbits);  mpz_setbit (p, pbits - 1);  mpz_nextprime (p, p);
        gmp_randclear (rs);
    }
    if ( ! mpz_probab_prime_p (p, 25) ) { fprintf (stderr, "p is not prime\n"); return 1; }

    unsigned long *ds;
    size_t nd = gen_fundamental (maxd, &ds);

    cornacchia_ctx ctx;
    double t0 = wall ();
    cornacchia_init (&ctx, p);
    double t_init = wall () - t0;

    mpz_init (t);  mpz_init (v);  mpz_init (chk);
    size_t n_solv = 0;
    double ts = wall ();
    for ( size_t i = 0 ; i < nd ; i++ ) {
        unsigned long d = ds[i];
        if ( cornacchia_solve (&ctx, d, t, v) ) {
            // self-check: 4p == t^2 + d v^2
            mpz_mul (chk, t, t);                         // t^2
            mpz_mul (ctx.t1, v, v);  mpz_mul_ui (ctx.t1, ctx.t1, d);  // d v^2
            mpz_add (chk, chk, ctx.t1);
            if ( mpz_cmp (chk, ctx.p4) != 0 ) {
                gmp_fprintf (stderr, "SELF-CHECK FAILED d=%lu t=%Zd v=%Zd\n", d, t, v);
                return 2;
            }
            n_solv++;
            if ( dump ) gmp_printf ("%lu 1 %Zd %Zd\n", d, t, v);
            if ( verbose ) gmp_printf ("  D=-%lu  t=%Zd v=%Zd  ->  curve orders p+1-+t\n", d, t, v);
        } else {
            if ( dump ) printf ("%lu 0\n", d);
        }
    }
    double dt = wall () - ts;

    if ( ! dump ) {
        gmp_fprintf (stderr, "p (%lu bits) = %Zd\n", (unsigned long) mpz_sizeinbase (p, 2), p);
        fprintf (stderr, "fundamental D with |D| <= %lu : %zu\n", maxd, nd);
        fprintf (stderr, "init: %.3f ms\n", 1e3 * t_init);
        fprintf (stderr, "search: %.3f s  (%.3f us / D, %.0f D/s)\n",
                 dt, 1e6 * dt / nd, nd / dt);
        fprintf (stderr, "  legendre-pass (sqrt done): %lu (%.1f%%)\n",
                 ctx.n_sqrt, 100.0 * ctx.n_sqrt / nd);
        fprintf (stderr, "  solvable: %zu (1 per %.0f D)\n", n_solv, (double) nd / (n_solv ? n_solv : 1));
    }

    mpz_clears (p, t, v, chk, NULL);
    cornacchia_clear (&ctx);
    free (ds);
    return 0;
}
