#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include "zp_poly.h"
#include "fproot.h"

/*
    Head-to-head benchmark + cross-validation: fproot's fp_find_root (Montgomery +
    Kronecker squaring + Barrett reduction) vs zp_poly's zp_poly_find_split_root
    (Harvey-Sutherland, tuned mpz engine) on the SAME completely-split monic h.

      zpbench pbits mind maxd [step-factor]

    For each degree d: build h = prod (x - r_i) from d distinct random roots, run
    both root-finders on copies, check each returned root is in {r_i}, and report
    wall time.  The fproot timing includes the [0,p) -> Montgomery load and the
    root conversion back, since that is what real usage pays.
*/

static double wall (void)
{ struct timespec ts; clock_gettime (CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }

int main (int argc, char **argv)
{
    if ( argc < 4 ) { fprintf (stderr, "usage: zpbench pbits mind maxd [mulfac]\n"); return 2; }
    int pbits = atoi (argv[1]), mind = atoi (argv[2]), maxd = atoi (argv[3]);
    double fac = argc > 4 ? atof (argv[4]) : 2.0;

    gmp_randstate_t rs;  gmp_randinit_default (rs);  gmp_randseed_ui (rs, 20260701);
    mpz_t P;  mpz_init (P);
    mpz_urandomb (P, rs, pbits);  mpz_setbit (P, pbits-1);  mpz_nextprime (P, P);
    gmp_printf ("p = %Zd (%d bits)\n", P, pbits);

    zp_mod_t zmod;  zp_mod_init (zmod, P);
    fp_ctx C;  fp_init (&C, P);

    printf ("%6s %14s %14s %9s\n", "deg", "fproot(s)", "zp_poly(s)", "ratio");
    for ( int d = mind ; d <= maxd ; d = (int)(d * fac) ) {
        // distinct random roots
        mpz_t *roots = malloc (d * sizeof(mpz_t));
        for ( int i = 0 ; i < d ; i++ ) { mpz_init (roots[i]); mpz_urandomm (roots[i], rs, P); }
        // h = prod (x - r_i) via zp_poly's product tree (fast, and independent of both)
        zp_t *h = zp_poly_alloc (d, zmod->p);
        zp_poly_from_roots (h, roots, d, zmod->p);

        // --- fproot ---
        double t0 = wall ();
        fp_poly hf;  fpoly_init (&C, &hf, d + 2);
        for ( int i = 0 ; i <= d ; i++ ) fp_set_mpz (&C, hf.c + (size_t)i*C.s, h[i]);
        hf.deg = d;
        mp_limb_t rm[64];
        fp_find_root (&C, &hf, rm, 0xFEED + d);
        mpz_t r1;  mpz_init (r1);  fp_get_mpz (&C, r1, rm);
        double tf = wall () - t0;
        fpoly_clear (&hf);

        // --- zp_poly ---
        zp_t *hc = zp_poly_alloc (d, zmod->p);
        for ( int i = 0 ; i <= d ; i++ ) mpz_set (hc[i], h[i]);
        mpz_t r2;  mpz_init (r2);
        t0 = wall ();
        zp_poly_find_split_root (r2, hc, d, zmod);
        double tz = wall () - t0;
        zp_poly_free (hc, d);

        // --- validate: both roots in {r_i} ---
        int ok1 = 0, ok2 = 0;
        for ( int i = 0 ; i < d ; i++ ) {
            if ( ! mpz_cmp (r1, roots[i]) ) ok1 = 1;
            if ( ! mpz_cmp (r2, roots[i]) ) ok2 = 1;
        }
        printf ("%6d %11.4f %s %11.4f %s %8.2fx\n", d, tf, ok1?"OK":"BAD", tz, ok2?"OK":"BAD", tf/tz);
        fflush (stdout);

        zp_poly_free (h, d);
        for ( int i = 0 ; i < d ; i++ ) mpz_clear (roots[i]);
        free (roots);  mpz_clear (r1);  mpz_clear (r2);
    }
    fp_clear (&C);  zp_mod_clear (zmod);  mpz_clear (P);  gmp_randclear (rs);
    return 0;
}
