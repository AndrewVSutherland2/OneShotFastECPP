#ifndef _CURVE_H_
#define _CURVE_H_

#include <stdint.h>
#include <gmp.h>
#include "fproot.h"
#include "cornacchia.h"

/*
    Certificate assembly: from a CM j-invariant j0 (curve order N = p+1 -/+ t) and
    a smooth m | N, produce a Montgomery coefficient A and an x-coordinate x0 such
    that the point (x0) has order exactly m on the Montgomery curve E_A (or its
    twist) -- i.e. the (A, x0, m, ...) part of a one-shot ECPP certificate that
    voneshot.py accepts.  All arithmetic is x-only Montgomery over F_p (fp_ctx).

      j0 -> A:   solve 256(A^2-3)^3 = j0 (A^2-4) for A^2 = u, A = sqrt(u).
      order:     determine #E_A (p+1-t or p+1+t) with the ladder; pick E_A or its
                 twist so the working curve has order N.
      x0:        random x with the right quadratic character; x0 = [N/m] x; check
                 [m] x0 = O and [m/l] x0 != O for each prime l | m.

    Returns 1 with A_out, x0_out in [0,p) on success; 0 if no valid (A, x0) exists
    (e.g. E_A has no Montgomery model, or m does not divide the group exponent).
*/
int mont_assemble (const fp_ctx *C, cornacchia_ctx *cc, const mpz_t j0, const mpz_t N,
                   const mpz_t t, const mpz_t m, mpz_t A_out, mpz_t x0_out, uint64_t seed);

#endif
