#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "curve.h"
#include "smooth.h"          // factor_smooth

#define ML 64                // matches fproot MAXLIMB (supports p up to ~1900 bits)

/* ---- Montgomery x-only arithmetic on E_A: y^2 = x^3 + A x^2 + x (X:Z coords) ---- *
 * Formulas depend only on A and are valid on E_A and its twist (Kummer line).      */

static void xdbl (const fp_ctx *C, mp_limb_t *X2, mp_limb_t *Z2,
                  const mp_limb_t *X, const mp_limb_t *Z, const mp_limb_t *A, const mp_limb_t *four)
{
    int s = C->s;  mp_limb_t XX[ML], ZZ[ML], XZ[ML], t[ML], u[ML];
    fp_mul (C, XX, X, X);  fp_mul (C, ZZ, Z, Z);  fp_mul (C, XZ, X, Z);
    fp_sub (C, t, XX, ZZ);  fp_mul (C, t, t, t);                       // (X^2-Z^2)^2
    fp_mul (C, u, A, XZ);  fp_add (C, u, u, XX);  fp_add (C, u, u, ZZ);// X^2 + A XZ + Z^2
    fp_mul (C, u, u, XZ);  fp_mul (C, u, u, four);                     // 4 XZ (...)
    mpn_copyi (X2, t, s);  mpn_copyi (Z2, u, s);
}

static void xadd (const fp_ctx *C, mp_limb_t *X3, mp_limb_t *Z3,
                  const mp_limb_t *X1, const mp_limb_t *Z1, const mp_limb_t *X2, const mp_limb_t *Z2,
                  const mp_limb_t *Xd, const mp_limb_t *Zd)
{
    int s = C->s;  mp_limb_t a[ML], b[ML], t[ML], u[ML];
    fp_sub (C, a, X1, Z1);  fp_add (C, t, X2, Z2);  fp_mul (C, a, a, t);   // (X1-Z1)(X2+Z2)
    fp_add (C, b, X1, Z1);  fp_sub (C, t, X2, Z2);  fp_mul (C, b, b, t);   // (X1+Z1)(X2-Z2)
    fp_add (C, t, a, b);  fp_mul (C, t, t, t);  fp_mul (C, t, t, Zd);      // Zd (a+b)^2
    fp_sub (C, u, a, b);  fp_mul (C, u, u, u);  fp_mul (C, u, u, Xd);      // Xd (a-b)^2
    mpn_copyi (X3, t, s);  mpn_copyi (Z3, u, s);
}

// (Xk:Zk) = [k](XP:ZP)   (Montgomery ladder; difference point is P)
static void ladder (const fp_ctx *C, mp_limb_t *Xk, mp_limb_t *Zk, const mpz_t k,
                    const mp_limb_t *XP, const mp_limb_t *ZP, const mp_limb_t *A, const mp_limb_t *four)
{
    int s = C->s;
    if ( mpz_sgn (k) == 0 ) { mpn_copyi (Xk, C->R1, s);  memset (Zk, 0, (size_t)s*sizeof(mp_limb_t)); return; }
    int nb = (int) mpz_sizeinbase (k, 2);
    if ( nb == 1 ) { mpn_copyi (Xk, XP, s);  mpn_copyi (Zk, ZP, s); return; }
    mp_limb_t X0[ML], Z0[ML], X1[ML], Z1[ML], tX[ML], tZ[ML];
    mpn_copyi (X0, XP, s);  mpn_copyi (Z0, ZP, s);          // 1P
    xdbl (C, X1, Z1, XP, ZP, A, four);                      // 2P
    for ( int i = nb - 2 ; i >= 0 ; i-- ) {
        if ( ! mpz_tstbit (k, i) ) {
            xadd (C, tX, tZ, X0, Z0, X1, Z1, XP, ZP);  mpn_copyi (X1, tX, s);  mpn_copyi (Z1, tZ, s);
            xdbl (C, tX, tZ, X0, Z0, A, four);          mpn_copyi (X0, tX, s);  mpn_copyi (Z0, tZ, s);
        } else {
            xadd (C, tX, tZ, X0, Z0, X1, Z1, XP, ZP);  mpn_copyi (X0, tX, s);  mpn_copyi (Z0, tZ, s);
            xdbl (C, tX, tZ, X1, Z1, A, four);          mpn_copyi (X1, tX, s);  mpn_copyi (Z1, tZ, s);
        }
    }
    mpn_copyi (Xk, X0, s);  mpn_copyi (Zk, Z0, s);
}

// [k](x:1) == O ?  (point at infinity has Z-coordinate 0)
static int ladder_is_O (const fp_ctx *C, const mpz_t k, const mp_limb_t *Xp, const mp_limb_t *A, const mp_limb_t *four)
{
    mp_limb_t Xk[ML], Zk[ML];
    ladder (C, Xk, Zk, k, Xp, C->R1, A, four);
    return fp_is_zero (C, Zk);
}

int mont_assemble (const fp_ctx *C, cornacchia_ctx *cc, const mpz_t j0, const mpz_t N,
                   const mpz_t t, const mpz_t m, mpz_t A_out, mpz_t x0_out, uint64_t seed)
{
    int s = C->s;
    mp_limb_t four[ML];  fp_set_ui (C, four, 4);
    mpz_t p;  mpz_init_set (p, C->pz);                      // local copy of the modulus

    // 1. A from j0: roots u of  u^3 - 9u^2 + (6912-j0)/256 u + (4j0-6912)/256 = 0,  A = sqrt(u).
    mpz_t inv256, c0, c1, tmp, Az, uz;  mpz_inits (inv256, c0, c1, tmp, Az, uz, NULL);
    mpz_set_ui (inv256, 256);  mpz_invert (inv256, inv256, p);
    mpz_ui_sub (c1, 6912, j0);  mpz_mul (c1, c1, inv256);  mpz_mod (c1, c1, p);   // (6912-j0)/256
    mpz_mul_ui (c0, j0, 4);  mpz_sub_ui (c0, c0, 6912);  mpz_mul (c0, c0, inv256); mpz_mod (c0, c0, p);
    fp_poly cub;  fpoly_init (C, &cub, 4);
    fp_set_mpz (C, cub.c + 0*s, c0);  fp_set_mpz (C, cub.c + 1*s, c1);
    mpz_set_si (tmp, -9);  mpz_mod (tmp, tmp, p);  fp_set_mpz (C, cub.c + 2*s, tmp);
    mpn_copyi (cub.c + 3*s, C->R1, s);  cub.deg = 3;
    mp_limb_t uroots[3*ML];  int nu = fp_find_all_roots (C, &cub, uroots, 3, seed);
    fpoly_clear (&cub);
    int haveA = 0;
    mp_limb_t Am[ML];
    for ( int i = 0 ; i < nu ; i++ ) {
        fp_get_mpz (C, uz, uroots + (size_t)i*s);
        if ( mpz_cmp_ui (uz, 4) == 0 ) continue;                     // A^2=4 -> A=+-2, singular
        if ( ! cornacchia_sqrtmodp (cc, Az, uz) ) continue;         // u must be a QR to get A=sqrt(u)
        mpz_sub_ui (tmp, p, 2);
        if ( mpz_cmp_ui (Az, 2) == 0 || mpz_cmp (Az, tmp) == 0 ) continue;
        fp_set_mpz (C, Am, Az);  haveA = 1;  break;
    }
    if ( ! haveA ) { mpz_clears (inv256, c0, c1, tmp, Az, uz, NULL); mpz_clear (p); return 0; }   // no Montgomery model

    // 2. #E_A: pick a point on E_A (f(x) a QR) and see which of p+1-/+t annihilates it.
    mpz_t np1mt, np1pt, x, fx;  mpz_inits (np1mt, np1pt, x, fx, NULL);
    mpz_add_ui (np1mt, p, 1);  mpz_sub (np1mt, np1mt, t);
    mpz_add_ui (np1pt, p, 1);  mpz_add (np1pt, np1pt, t);
    gmp_randstate_t rs;  gmp_randinit_default (rs);  gmp_randseed_ui (rs, seed ? seed : 1);
    int cardEA_is_mt = -1;                                            // 1 => #E_A=p+1-t, 0 => p+1+t
    mp_limb_t Xp[ML];
    for ( int tr = 0 ; tr < 60 && cardEA_is_mt < 0 ; tr++ ) {
        mpz_urandomm (x, rs, p);
        mpz_mul (fx, x, x);  mpz_add (fx, fx, Az);  mpz_mul (fx, fx, x);  mpz_add (fx, fx, x);  mpz_mod (fx, fx, p);
        if ( mpz_legendre (fx, p) != 1 ) continue;                  // want a point on E_A
        fp_set_mpz (C, Xp, x);
        int a = ladder_is_O (C, np1mt, Xp, Am, four);
        int b = ladder_is_O (C, np1pt, Xp, Am, four);
        if ( a && ! b ) cardEA_is_mt = 1;  else if ( b && ! a ) cardEA_is_mt = 0;   // else ambiguous: retry
    }
    if ( cardEA_is_mt < 0 ) { gmp_randclear (rs); mpz_clears (inv256,c0,c1,tmp,Az,uz,np1mt,np1pt,x,fx,NULL); mpz_clear (p); return 0; }

    // Which quadratic character puts a point on the order-N curve (E_A or its twist)?
    int cardEA_matches_N = ( mpz_cmp (cardEA_is_mt ? np1mt : np1pt, N) == 0 );
    int target_leg = cardEA_matches_N ? 1 : -1;                      // f(x) must have this Legendre symbol

    // 3. point of order exactly m on the order-N curve
    uint64_t mp_[64];  int mex[64];  int nmp = factor_smooth (m, mp_, mex, 64);
    mpz_t Nm, ml;  mpz_inits (Nm, ml, NULL);
    mpz_divexact (Nm, N, m);                                         // cofactor N/m
    int ok = 0;
    mp_limb_t X0[ML], Z0[ML], Zk[ML], Zinv[ML], x0m[ML];
    for ( int tr = 0 ; tr < 400 && ! ok ; tr++ ) {
        mpz_urandomm (x, rs, p);
        mpz_mul (fx, x, x);  mpz_add (fx, fx, Az);  mpz_mul (fx, fx, x);  mpz_add (fx, fx, x);  mpz_mod (fx, fx, p);
        if ( mpz_legendre (fx, p) != target_leg ) continue;        // put x on the order-N curve
        fp_set_mpz (C, Xp, x);
        ladder (C, X0, Z0, Nm, Xp, C->R1, Am, four);                // [N/m] (x:1)
        if ( fp_is_zero (C, Z0) ) continue;                         // landed on O
        // affine x-coord of the candidate point; then order test
        mp_limb_t Xk[ML];
        ladder (C, Xk, Zk, m, X0, Z0, Am, four);                    // [m] P
        if ( ! fp_is_zero (C, Zk) ) continue;                       // ord does not divide m
        int good = 1;
        for ( int i = 0 ; i < nmp && good ; i++ ) {
            mpz_divexact_ui (ml, m, mp_[i]);                        // m/l
            ladder (C, Xk, Zk, ml, X0, Z0, Am, four);
            if ( fp_is_zero (C, Zk) ) good = 0;                     // [m/l]P = O -> order < m
        }
        if ( ! good ) continue;
        fp_inv (C, Zinv, Z0);  fp_mul (C, x0m, X0, Zinv);           // x0 = X0/Z0 (affine)
        fp_get_mpz (C, x0_out, x0m);
        mpz_set (A_out, Az);                                        // A in [0,p)
        ok = 1;
    }

    gmp_randclear (rs);
    mpz_clears (inv256, c0, c1, tmp, Az, uz, np1mt, np1pt, x, fx, Nm, ml, NULL);  mpz_clear (p);
    return ok;
}
