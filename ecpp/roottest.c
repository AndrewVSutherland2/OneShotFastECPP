#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include "fproot.h"

/*
    Validate fp_find_root: build h = prod (x - r_i) over F_p from known distinct
    roots r_i, find a root, and check (a) it is one of the r_i, (b) it really is a
    root of h (independent Horner evaluation mod p), (c) [subset] PARI's root set
    of h mod p equals {r_i}.  Also benchmarks.

      roottest [pbits] [deg] [ntrials]        # random test at one size
      roottest pari  [pbits] [deg]            # + PARI cross-check
      roottest bench                          # timing sweep
*/

static double wall (void)
{ struct timespec ts; clock_gettime (CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }

// build monic h = prod(x - roots[i]) mod p, coefficients into hc[0..d] (mpz).
// Multiply h (degree cur) by (x - r): new_j = old_{j-1} - r*old_j, high->low.
static void build_from_roots (mpz_t *hc, const mpz_t *roots, int d, const mpz_t p)
{
    mpz_set_ui (hc[0], 1);  int cur = 0;                    // h = 1
    mpz_t t;  mpz_init (t);
    for ( int i = 0 ; i < d ; i++ ) {
        mpz_set_ui (hc[cur + 1], 0);                        // extend to degree cur+1
        for ( int j = cur + 1 ; j >= 0 ; j-- ) {
            mpz_mul (t, roots[i], hc[j]);                   // r * old_j
            if ( j >= 1 ) mpz_sub (hc[j], hc[j-1], t);      // old_{j-1} - r*old_j
            else mpz_neg (hc[j], t);                        // j=0: -r*old_0
            mpz_mod (hc[j], hc[j], p);
        }
        cur++;
    }
    mpz_clear (t);
}

// Horner: eval sum hc[i] x^i at x mod p
static void poly_eval (mpz_t out, const mpz_t *hc, int d, const mpz_t x, const mpz_t p)
{
    mpz_set (out, hc[d]);
    for ( int i = d - 1 ; i >= 0 ; i-- ) { mpz_mul (out, out, x); mpz_add (out, out, hc[i]); mpz_mod (out, out, p); }
}

static int run (unsigned long pbits, int d, int ntrials, int do_pari, uint64_t seed0)
{
    gmp_randstate_t rs;  gmp_randinit_default (rs);  gmp_randseed_ui (rs, 424242 + seed0);
    mpz_t p;  mpz_init (p);  mpz_urandomb (p, rs, pbits);  mpz_setbit (p, pbits-1);  mpz_nextprime (p, p);

    fp_ctx C;  fp_init (&C, p);
    mpz_t *roots = malloc (d * sizeof(mpz_t));
    mpz_t *hc = malloc ((d+1) * sizeof(mpz_t));
    for ( int i = 0 ; i < d ; i++ ) mpz_init (roots[i]);
    for ( int i = 0 ; i <= d ; i++ ) mpz_init (hc[i]);

    int fails = 0;  double tsum = 0;
    for ( int tr = 0 ; tr < ntrials ; tr++ ) {
        // distinct random roots
        for ( int i = 0 ; i < d ; i++ ) {
            int dup;
            do { dup = 0; mpz_urandomm (roots[i], rs, p);
                 for ( int k = 0 ; k < i ; k++ ) if ( !mpz_cmp (roots[i], roots[k]) ) { dup = 1; break; } }
            while ( dup );
        }
        build_from_roots (hc, roots, d, p);

        fp_poly h;  fpoly_init (&C, &h, d + 2);
        for ( int i = 0 ; i <= d ; i++ ) fp_set_mpz (&C, h.c + (size_t)i*C.s, hc[i]);
        h.deg = d;

        mp_limb_t rootm[64];
        double t0 = wall ();
        fp_find_root (&C, &h, rootm, 0x1234567 + tr*2654435761u + seed0);
        tsum += wall () - t0;
        fpoly_clear (&h);

        mpz_t rz, ev;  mpz_init (rz);  mpz_init (ev);
        fp_get_mpz (&C, rz, rootm);
        // (a) membership
        int member = 0;  for ( int i = 0 ; i < d ; i++ ) if ( !mpz_cmp (rz, roots[i]) ) { member = 1; break; }
        // (b) independent root check
        poly_eval (ev, hc, d, rz, p);
        int isroot = (mpz_sgn (ev) == 0);
        if ( !member || !isroot ) {
            fails++;
            gmp_printf ("  FAIL trial %d: root=%Zd member=%d isroot=%d\n", tr, rz, member, isroot);
        }
        mpz_clear (rz);  mpz_clear (ev);
    }

    // (c) PARI cross-check on the last h
    int pari_ok = 1;
    if ( do_pari ) {
        char *cmd = malloc (64 + (size_t)(d+1)*(pbits/3+8));
        int off = sprintf (cmd, "p=");
        off += gmp_sprintf (cmd+off, "%Zd", p);
        off += sprintf (cmd+off, ";h=Mod(1,p)*(");
        for ( int i = d ; i >= 0 ; i-- ) { char *b; gmp_asprintf (&b, "%s%Zd*x^%d", i==d?"":"+", hc[i], i); off += sprintf (cmd+off, "%s", b); free (b); }
        off += sprintf (cmd+off, ");v=vecsort(apply(z->lift(z),polrootsmod(h)));print(#v);for(i=1,#v,print(v[i]))");
        FILE *f = fopen ("/tmp/rt_pari.gp", "w");  fprintf (f, "%s\n", cmd);  fclose (f);
        FILE *out = popen ("gp -q < /tmp/rt_pari.gp 2>/dev/null", "r");
        int npar = -1;  if ( out ) {
            if ( fscanf (out, "%d", &npar) == 1 ) {
                mpz_t pr;  mpz_init (pr);  int matched = 0;
                for ( int i = 0 ; i < npar ; i++ ) { if ( gmp_fscanf (out, "%Zd", pr) == 1 ) {
                    for ( int k = 0 ; k < d ; k++ ) if ( !mpz_cmp (pr, roots[k]) ) { matched++; break; } } }
                pari_ok = (npar == d && matched == d);
                mpz_clear (pr);
            } else pari_ok = 0;
            pclose (out);
        }
        printf ("  PARI: rootset size %d (expected %d), %s\n", npar, d, pari_ok ? "matches {r_i}" : "MISMATCH");
        free (cmd);
    }

    printf ("p=%lub d=%d trials=%d: %s   (avg %.4f s/root)\n",
            pbits, d, ntrials, (fails==0 && pari_ok) ? "OK" : "FAIL", tsum/ntrials);

    for ( int i = 0 ; i < d ; i++ ) mpz_clear (roots[i]);
    for ( int i = 0 ; i <= d ; i++ ) mpz_clear (hc[i]);
    free (roots);  free (hc);  fp_clear (&C);  mpz_clear (p);  gmp_randclear (rs);
    return (fails==0 && pari_ok) ? 0 : 1;
}

int main (int argc, char **argv)
{
    if ( argc > 1 && !strcmp (argv[1], "bench") ) {
        for ( int pb = 256 ; pb <= 384 ; pb += 128 )
            for ( int d = 10 ; d <= 640 ; d *= 2 ) run (pb, d, 1, 0, d);
        return 0;
    }
    if ( argc > 1 && !strcmp (argv[1], "pari") ) {
        unsigned long pb = argc > 2 ? strtoul (argv[2], 0, 10) : 256;
        int d = argc > 3 ? atoi (argv[3]) : 40;
        return run (pb, d, 1, 1, 0);
    }
    unsigned long pb = argc > 1 ? strtoul (argv[1], 0, 10) : 256;
    int d = argc > 2 ? atoi (argv[2]) : 40;
    int nt = argc > 3 ? atoi (argv[3]) : 20;
    return run (pb, d, nt, 0, 0);
}
