#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include "fproot.h"

/*
    End-to-end component (c): compute the class polynomial H_D mod p (via the
    validated classpoly binary), then find a root j0 of H_D over F_p with the
    equal-degree-splitting root finder.  For a discriminant D with 4p = t^2+|D|v^2
    solvable (as produced by dscan), H_D splits completely over F_p, so j0 is a
    valid j-invariant of a curve with CM by D -- the seed for building the curve.

      hdroot D=<disc<-4> (p=<dec> | pbits=<n> [seed=<s>]) [inv=0]

    Requires the classpoly environment (CLASSPOLY_PHI_DIR etc); source setenv.sh,
    or the vars must already be set.  Cross-checks the root by Horner evaluation.
*/

static double wall (void)
{ struct timespec ts; clock_gettime (CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }

static void gen_p (mpz_t p, unsigned long pbits, unsigned long seed)
{
    gmp_randstate_t rs;  gmp_randinit_default (rs);  gmp_randseed_ui (rs, seed);
    mpz_urandomb (p, rs, pbits);  mpz_setbit (p, pbits - 1);  mpz_nextprime (p, p);
    gmp_randclear (rs);
}

int main (int argc, char **argv)
{
    long D = 0;  int inv = 0;
    unsigned long pbits = 0, seed = 1;
    mpz_t p;  mpz_init (p);
    for ( int i = 1 ; i < argc ; i++ ) {
        if ( ! strncmp (argv[i], "D=", 2) ) D = strtol (argv[i]+2, 0, 10);
        else if ( ! strncmp (argv[i], "p=", 2) ) mpz_set_str (p, argv[i]+2, 10);
        else if ( ! strncmp (argv[i], "pbits=", 6) ) pbits = strtoul (argv[i]+6, 0, 10);
        else if ( ! strncmp (argv[i], "seed=", 5) ) seed = strtoul (argv[i]+5, 0, 10);
        else if ( ! strncmp (argv[i], "inv=", 4) ) inv = atoi (argv[i]+4);
    }
    if ( D >= -4 ) { fprintf (stderr, "need D < -4\n"); return 2; }
    if ( ! mpz_sgn (p) ) { if ( ! pbits ) pbits = 256;  gen_p (p, pbits, seed); }
    if ( ! mpz_probab_prime_p (p, 25) ) { fprintf (stderr, "p not prime\n"); return 1; }

    // 1. compute H_D mod p via classpoly binary -> temp file
    char pdec[4096];  gmp_snprintf (pdec, sizeof pdec, "%Zd", p);
    char outf[256];   snprintf (outf, sizeof outf, "/tmp/hdroot_%ld_%d.txt", -D, inv);
    char cmd[8192];   snprintf (cmd, sizeof cmd, "classpoly %ld %d %s %s -1", D, inv, pdec, outf);
    double t0 = wall ();
    if ( system (cmd) != 0 ) { fprintf (stderr, "classpoly failed (env set? source setenv.sh)\n"); return 1; }
    double t_hd = wall () - t0;

    // 2. parse coefficients (lines "c*X^e + "); low-degree first, monic top term
    FILE *f = fopen (outf, "r");
    if ( ! f ) { fprintf (stderr, "cannot read %s\n", outf); return 1; }
    int cap = 64, deg = -1;
    mpz_t *hc = malloc (cap * sizeof(mpz_t));
    char line[8192];
    while ( fgets (line, sizeof line, f) ) {
        char *x = strstr (line, "*X^");
        if ( ! x ) continue;
        int e = atoi (x + 3);
        if ( e >= cap ) { int nc = cap; while ( e >= nc ) nc *= 2; hc = realloc (hc, nc*sizeof(mpz_t)); cap = nc; }
        while ( deg < e ) { deg++; mpz_init (hc[deg]); }     // init up to e
        *x = 0;  mpz_set_str (hc[e], line, 10);
    }
    fclose (f);
    if ( deg < 1 ) { fprintf (stderr, "degenerate H_D (deg %d)\n", deg); return 1; }

    // 3. load into fp_poly (Montgomery) and find a root
    fp_ctx C;  fp_init (&C, p);
    fp_poly h;  fpoly_init (&C, &h, deg + 2);
    for ( int i = 0 ; i <= deg ; i++ ) { mpz_mod (hc[i], hc[i], p); fp_set_mpz (&C, h.c + (size_t)i*C.s, hc[i]); }
    h.deg = deg;

    mp_limb_t rootm[64];
    t0 = wall ();
    fp_find_root (&C, &h, rootm, 0xC0FFEE + seed);
    double t_root = wall () - t0;

    // 4. verify: Horner eval of H_D at the root == 0 mod p
    mpz_t j, ev;  mpz_init (j);  mpz_init (ev);
    fp_get_mpz (&C, j, rootm);
    mpz_set (ev, hc[deg]);
    for ( int i = deg - 1 ; i >= 0 ; i-- ) { mpz_mul (ev, ev, j); mpz_add (ev, ev, hc[i]); mpz_mod (ev, ev, p); }
    int ok = (mpz_sgn (ev) == 0);

    gmp_printf ("D=%ld  p=%lu bits  h(D)=deg=%d  inv=%d\n", D, (unsigned long) mpz_sizeinbase(p,2), deg, inv);
    printf ("H_D mod p (classpoly): %.3f s\n", t_hd);
    printf ("root-find (EDS+Barrett): %.4f s\n", t_root);
    gmp_printf ("root j0 = %Zd\n", j);
    printf ("verify H_D(j0) == 0 mod p: %s\n", ok ? "OK" : "FAIL");

    for ( int i = 0 ; i <= deg ; i++ ) mpz_clear (hc[i]);
    free (hc);  mpz_clear (j);  mpz_clear (ev);  fpoly_clear (&h);  fp_clear (&C);  mpz_clear (p);
    return ok ? 0 : 3;
}
