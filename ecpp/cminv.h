#ifndef _CMINV_H_
#define _CMINV_H_

#include <gmp.h>

/*
    Convert a class-invariant value X (with classpoly invariant code `inv`) to the
    j-invariant(s) of the corresponding elliptic curve over F_p (p a large prime
    in an mpz_t).  This is the big-F_p analog of classpoly's word-size invtoj /
    ff_all_j_from_inv, i.e. the mpz_j_from_inv the header declares.

    The hardwired families (gamma2, Weber f/f2/f3/f4/f8, t/t2/t6, u/u2/u8) use the
    exact pure-mpz formulas from A. Sutherland's class_inv_mpz.c (classpoly, GPL).
    The modular-polynomial families (Atkin A_N, single/double eta) apply the same
    mpz_inv_power and then find a root of Phi_inv(X', J) over F_p -- classpoly does
    this with the zp_poly library; here we use fproot (see invj.{c,h}), which is
    equivalent and self-contained.

    Writes up to maxj candidate j-invariants into J[] (in [0,p)); returns the count
    (1 for the hardwired invariants, 1-2 for the modular-polynomial ones -- the CM
    trace test in cm_method selects the correct one).  Returns -1 on error.
    phidir = CLASSPOLY_PHI_DIR (for the Phi_inv modular polynomials).

    cm_j_from_inv is the production entry point: it calls classpoly's authoritative
    mpz_j_from_inv (class_inv_mpz.c, linked against the zp_poly library) for every
    invariant that function supports, falling back to the fproot-based path only for
    the generic single-eta range (400..499) that class_inv_mpz.c does not handle.
    cm_j_from_inv_ref is the self-contained fproot-based implementation (Sutherland's
    formulas re-ported + Phi_inv root-find), kept as a cross-validation reference.
*/
int cm_j_from_inv     (mpz_t *J, int maxj, mpz_t X, mpz_t P, int inv, const char *phidir);
int cm_j_from_inv_ref (mpz_t *J, int maxj, mpz_t X, mpz_t P, int inv, const char *phidir);

#endif
