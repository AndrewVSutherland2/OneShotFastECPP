#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cminv.h"

/* cmjtest p inv f0 : print the j-invariant(s) cm_j_from_inv derives from the
   invariant value f0 (code inv) over F_p.  Cross-checked against classpoly's
   word-size invtoj by test_cminv.py.  Needs CLASSPOLY_PHI_DIR. */
int main (int argc, char **argv)
{
    if ( argc < 4 ) { fprintf (stderr, "usage: cmjtest p inv f0\n"); return 2; }
    mpz_t p, f0;  mpz_init_set_str (p, argv[1], 10);  int inv = atoi (argv[2]);
    mpz_init_set_str (f0, argv[3], 10);
    const char *phidir = getenv ("CLASSPOLY_PHI_DIR");
    mpz_t J[8];  for ( int i = 0 ; i < 8 ; i++ ) mpz_init (J[i]);
    int nj = cm_j_from_inv (J, 8, f0, p, inv, phidir);
    if ( nj < 0 ) { fprintf (stderr, "cm_j_from_inv failed (inv=%d)\n", inv); return 1; }
    for ( int i = 0 ; i < nj ; i++ ) gmp_printf ("%Zd ", J[i]);
    printf ("\n");
    for ( int i = 0 ; i < 8 ; i++ ) mpz_clear (J[i]);
    mpz_clear (p);  mpz_clear (f0);
    return 0;
}
