#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "fproot.h"
#include "cornacchia.h"
#include "curve.h"

/* asmtest p j0 N t m : run mont_assemble and print "A x0" (or FAIL). */
int main (int argc, char **argv)
{
    if ( argc < 6 ) { fprintf (stderr, "usage: asmtest p j0 N t m\n"); return 2; }
    mpz_t p, j0, N, t, m, A, x0;
    mpz_init_set_str (p, argv[1], 10);  mpz_init_set_str (j0, argv[2], 10);
    mpz_init_set_str (N, argv[3], 10);  mpz_init_set_str (t, argv[4], 10);
    mpz_init_set_str (m, argv[5], 10);  mpz_init (A);  mpz_init (x0);
    fp_ctx C;  fp_init (&C, p);
    cornacchia_ctx cc;  cornacchia_init (&cc, p);
    int ok = mont_assemble (&C, &cc, j0, N, t, m, A, x0, 0xBEEF15);
    if ( ok ) gmp_printf ("%Zd %Zd\n", A, x0);  else printf ("FAIL\n");
    fp_clear (&C);  cornacchia_clear (&cc);
    mpz_clears (p, j0, N, t, m, A, x0, NULL);
    return ok ? 0 : 1;
}
