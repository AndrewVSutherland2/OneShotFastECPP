#ifndef _FPROOT_H_
#define _FPROOT_H_

#include <stdint.h>
#include <gmp.h>

/*
    Root-finding of a class polynomial H_D over F_p, p a few-hundred-bit prime.

    H_D splits completely over F_p (p splits in the ring class field for the CM
    discriminant we selected), so there is no distinct-degree step: every
    irreducible factor is linear.  We go straight to equal-degree splitting for
    degree 1 (Cantor-Zassenhaus): pick random a, form b = (x+a)^((p-1)/2) mod h,
    split h by g = gcd(h, b-1), keep the smaller-degree side, recurse to degree 1.

    Inner loop (per the design):
      * (x+a)*f mod h with h = x^d + h'(x) monic: a shift, a scalar multiply and
        two additions -- NO polynomial multiplication.  Cheap.
      * f^2 mod h: the squaring dominates.  We form f^2 by Kronecker substitution
        (pack coefficients into one big integer, square it with GMP, unpack), then
        reduce coefficients mod p, then reduce the degree-(2d-2) result mod h.
      * All scalar mod-p arithmetic is in Montgomery form (multi-limb REDC).

    Field elements live in Montgomery form throughout; polynomials are flat arrays
    of s-limb coefficients.  The Montgomery limb count s is padded so that a
    Kronecker product coefficient (a sum of up to d products, < d*p^2) still lies
    below p*R and REDC returns it directly.
*/

// ---- Montgomery F_p context ----
typedef struct {
    int n;                 // bit length of p
    int s;                 // limb count of coefficients / R = 2^(64 s)  (padded)
    mp_limb_t *p;          // modulus, s limbs
    mp_limb_t  pinv;       // -p^{-1} mod 2^64
    mp_limb_t *R1;         // R mod p    (Montgomery form of 1), s limbs
    mp_limb_t *R2;         // R^2 mod p  (to convert into Montgomery form), s limbs
    mpz_t      pz;         // p as mpz (for setup / debugging)
} fp_ctx;

void fp_init  (fp_ctx *C, const mpz_t p);   // pads s so REDC covers Kronecker coeffs (< 2^34 * p^2)
void fp_clear (fp_ctx *C);

// s-limb Montgomery ops (operands and results are s limbs, in Montgomery form)
void fp_mul (const fp_ctx *C, mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b);
void fp_add (const fp_ctx *C, mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b);
void fp_sub (const fp_ctx *C, mp_limb_t *r, const mp_limb_t *a, const mp_limb_t *b);
void fp_set_mpz  (const fp_ctx *C, mp_limb_t *r, const mpz_t a);   // a in [0,p) -> Montgomery
void fp_get_mpz  (const fp_ctx *C, mpz_t r, const mp_limb_t *a);   // Montgomery -> [0,p)
void fp_set_ui   (const fp_ctx *C, mp_limb_t *r, unsigned long a);
int  fp_is_zero  (const fp_ctx *C, const mp_limb_t *a);
int  fp_equal    (const fp_ctx *C, const mp_limb_t *a, const mp_limb_t *b);
void fp_inv      (const fp_ctx *C, mp_limb_t *r, const mp_limb_t *a);   // Fermat: a^(p-2)

// ---- polynomial over F_p (flat array of s-limb Montgomery coefficients) ----
// A poly is (deg+1) coefficients c[0..deg], c[i] at limbs (buf + i*s).  deg = -1
// is the zero polynomial.  All routines assume enough capacity in the buffer.
typedef struct { mp_limb_t *c; int deg; int cap; } fp_poly;   // cap = #coeff slots

void fpoly_init  (const fp_ctx *C, fp_poly *f, int cap);
void fpoly_clear (fp_poly *f);

// Find one root of h (which must split completely over F_p).  h is monic of
// degree d>=1 (Montgomery-form coefficients), given in place and destroyed.
// Returns the root in Montgomery form (root, s limbs).  seed drives the RNG.
void fp_find_root (const fp_ctx *C, fp_poly *h, mp_limb_t *root, uint64_t seed);

// Find ALL distinct roots of h in F_p (h need NOT split completely).  Writes up
// to maxr roots (Montgomery form) into roots[i*s]; returns the count.  Works by
// g = gcd(h, x^p - x) then peeling linear factors.  h is not modified.
int  fp_find_all_roots (const fp_ctx *C, const fp_poly *h, mp_limb_t *roots, int maxr, uint64_t seed);

// low-level pieces reused by the certificate assembly:
void fp_pow (const fp_ctx *C, mp_limb_t *r, const mp_limb_t *a, const mpz_t e);  // Montgomery a^e

#endif
