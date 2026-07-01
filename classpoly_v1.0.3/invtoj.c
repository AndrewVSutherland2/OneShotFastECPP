#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "ff_poly.h"
#include "class_inv.h"
#include "mpzutil.h"
#include "cstd.h"

/*
    invtoj -- map class-invariant values to j-invariants over F_p.

    Usage:  invtoj p inv

    Sets up the prime field F_p (p must be a word-size prime), then reads
    invariant values x (integers in [0,p)) from stdin, one per line.  For each
    x it uses classpoly's f->j modular polynomial Phi_inv(x, J) to compute every
    j-invariant J in F_p with Phi_inv(x, J) = 0, and prints

        x j1 j2 ...

    (a line with just x if there are no roots in F_p).  This is the inverse of
    the j->invariant map classpoly uses internally; PARI/GP cannot do it for the
    invariants it does not support (Atkin A_41/47/59/71, the Ramanujan t, etc.),
    so it is the engine for Test 2 -- verifying that the roots of an invariant
    class polynomial mod p map onto the roots of the Hilbert class polynomial.
*/

int main (int argc, char *argv[])
{
    unsigned long p, xv;
    int inv, i, k;
    char line[512];
    ff_t x, jl[PHI_MAX_JDEG + 1];

    if (argc < 3) {
        printf("usage: invtoj p inv   (reads invariant values x mod p from stdin)\n");
        return 0;
    }
    p = strtoul(argv[1], 0, 10);
    inv = atoi(argv[2]);

    mpz_util_init();
    ff_setup_ui(p);
    if (inv && !inv_good_invariant(inv)) { err_printf("invalid invariant %d\n", inv); return 1; }
    inv_load_phi_fj(inv);              // load Phi_inv(f,j) from $CLASSPOLY_PHI_DIR/phi_<inv>_j.txt

    while (fgets(line, sizeof(line), stdin)) {
        if (line[0] == '\n' || line[0] == 0) continue;
        xv = strtoul(line, 0, 10);
        _ff_set_ui(x, xv);
        // every j in F_p with Phi_inv(x,J)=0 (power-aware, non-aborting)
        k = ff_all_j_from_inv(jl, &x, inv);
        printf("%lu", xv);
        for (i = 0; i < k; i++) printf(" %lu", _ff_get_ui(jl[i]));
        printf("\n");
    }
    return 0;
}
