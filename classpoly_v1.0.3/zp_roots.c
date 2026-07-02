/*

zp_roots.c: functions for computing rth roots and solving degree 2 and 3 equations in Fp=Z/pZ and Fp^2

Copyright (C) 2011-2012, Andrew V. Sutherland

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gmp.h>
#include "zp_poly_util.h"
#include "zp_poly_mul.h"
#include "zp_poly_main.h"
#include "zp2.h"
#include "zp_roots.h"

// computes b = a^(1/4) mod p.  Uses just one exponentiation when p = 3 mod 4.  Returns 0 if a is not a 4th-power residue, 1 ow.
int zp_fourth_root (zp_t b, zp_t a, zp_mod_t p)
{
    mpz_t e;
    zp_t c, d;
    int sts;
    
    // if p=1 mod 4 just iterate square roots - note that the choice of root doesn't matter since -1 is a QR
    if ( p->p1mod4 ) { if ( ! zp_sqrt (b, a, p) ) return 0; return zp_sqrt (b, b, p); }

    /*
        For p=3 mod 4, we may compute y=x^{1/4} as x^{(p+1)^2/16} since (p+1)^2/16 = 2^2/16 = 1/4 mod p-1.
        This yields y^4 = x^{(p+1)^2/4} = x^{(p-1)/2+1}^2 = x^{((p-1)/2)^2 + (p-1) + 1} = +/- x, 
        If y^4 = -x then x has no fourth root, since only one of x or -x has a square root for p=3 mod 4.
    */  
    mpz_init2 (e, 2*mpz_bits(p->p));  zp_init (c, p);  zp_init (d, p);
    mpz_add_ui (e, p->p, 1);  mpz_div_2exp (e, e, 2);  mpz_mul (e, e, e); mpz_mod (e, e, p->pm1);
    zp_set (c, a);
    zp_exp (b, a, e, p);
    zp_sqr (d, b, p);  zp_sqr (d, d, p);
    sts = zp_equal (c, d);
    zp_clear (c); zp_clear (d);
    mpz_clear (e);
    return sts;
}

// computes b = a^(1/6) mod p.  Uses just one exponentiation when p = 11 mod 12.   Returns 0 if a is not a 6th-power residue, 1 ow.
int zp_sixth_root (zp_t b, zp_t a, zp_mod_t p)
{
    // if p = 1 mod 4 or p = 1 mod 3 just iterate a sqrt and a cbrt
    if ( p->p1mod4 || p->p1mod3 ) { if ( ! zp_sqrt (b, a, p) ) return 0; return zp_cbrt (b, b, p); }
    puts ("fast sixth root not yet implmented");
    if ( ! zp_sqrt (b, a, p) ) return 0;
    return zp_cbrt (b, b, p); 
}

// computes b = a^(1/8) mod p.  Uses just one exponentiation when p = 3 mod 4.   Returns 0 if a is not a 8th-power residue, 1 ow.
int zp_eighth_root (zp_t b, zp_t a, zp_mod_t p)
{
    mpz_t e;
    zp_t c, d;
    int sts;
    
    // if p=1 mod 4 just iterate square roots - note that the choice of root doesn't matter since -1 is a QR
    if ( p->p1mod4 ) { if ( ! zp_sqrt (b, a, p) ) return 0; if ( ! zp_sqrt (b, b, p) ) return 0; return zp_sqrt (b, b, p); }

    /*
        For p=3 mod 4, we may compute y=x^{1/8} as x^{(p+1)^3/64} since (p+1)^3/64 = 2^3/64 = 1/8 mod p-1.
        This yields y^8 = x^{(p+1)^3/8} = x^{(p-1)/2+1}^3 = x^{((p-1)/2)^3 + 3*((p-1)/2)^2 + 3*(p-1)/2 + 1} = +/- x, 
        If y^8 = -x then x has no eighth root, since only one of x or -x has a square root for p=3 mod 4.
    */
    mpz_init2 (e, 2*mpz_bits(p->p));  zp_init (c, p);  zp_init (d, p);
    mpz_add_ui (d, p->p, 1);  mpz_div_2exp (d, d, 2);  mpz_mul (e, d, d); mpz_mod (e, e, p->pm1); mpz_mul (e, e, d); mpz_mod (d, d, p->pm1);
    zp_set (c, a);
    zp_exp (b, a, e, p);
    zp_sqr (d, b, p);  zp_sqr (d, d, p); zp_sqr (d, d, p);
    sts = zp_equal (c, d);
    zp_clear (c); zp_clear (d);
    mpz_clear (e);
    return sts;
    
}

// computes r = sqrt(d) mod p, returns 0 if none exists, 1 ow
// uses algorithm of S. Muller "On Probable Prime Testing and the Computation of Square Roots mod n," ANTS IV.
// however our goal is to find a sqrt, not test primality, so we only report compositeness if it prevents us from finding a sqrt
// this algorithm is typically 2 to 2.5 times slower for p=5 mod 8, but is much faster when the 2-Sylow is big and in the trivial case p=3mod4
int zp_sqrt_muller (zp_t c, zp_t a, zp_mod_t p)
{
    mpz_t n, h1, h0, w, x;
    int t, b, sts;

    if ( zp_is_zero(a) ) { zp_set_zero(c); return 1; }
    assert ( mpz_sgn(a) > 0 && mpz_cmp(a,p->p) < 0);        // just to be sure
    if ( mpz_jacobi (a, p->p) < 0) return 0;
    // p = 3 mod 4 case is a simple exponentiation o = a^((p+1)/4)
    if ( p->k == 1 ) {
        // it's ok to use p->w here, we aren't making any calls to z2 functions
        mpz_add_ui (p->w[0], p->p, 1);  mpz_div_2exp (p->w[0], p->w[0], 2);
        zp_exp (p->w[1], a, p->w[0], p);
        zp_sqr (p->w[0], p->w[1], p);
        if ( ! zp_equal (p->w[0], a) ) { gmp_printf ("Non-prime %Zd detected in zp_sqrt_muller\n", p); return -1; }
        zp_set (c, p->w[1]);
        return 1;
    }
    mpz_init_modp (n, p->p); mpz_init_modp_w (h0, p->p); mpz_init_modp_w (h1, p->p); mpz_init_modp_w (w, p->p);  mpz_init_modp_w (x, p->p);
    mpz_mul_2exp (x, a, 2);
    for ( b = 0 ; ; b++ ) {
        mpz_sub_ui (w, x, b*b);  mpz_neg(w, w);         // compute D = P^2-4Q where P=b is small and Q=a is the number whose sqrt we want
        if ( mpz_jacobi (w, p->p) == -1 ) break;
    }
    assert ( mpz_cmp_ui(p->p,b) > 0 );
    mpz_add_ui(n, p->p, 1); mpz_div_2exp (n, n, 1);     // n = (p->p+1)/2
    mpz_set_ui (h1,0); mpz_set_ui (h0,1);
    // compute x^{(p->p+1)/2} mod (x^2-Px+Q = x^2-bx+a)
    for ( t = mpz_sizeinbase (n, 2) - 1 ; t >= 0 ; t-- ) {
        // square h1*x+h0 mod (x^2-bx+a)
        mpz_mul_mod (w, h1, h1, p->p);
        mpz_mul (h1, h1, h0);  mpz_add (h1, h1, h1);  mpz_mul_ui (x, w, b);  mpz_add (h1, h1, x);  mpz_mod (h1, h1, p->p);
        mpz_mul (h0, h0, h0);  mpz_mul (w, w, a);  mpz_sub (h0, h0, w);  mpz_mod (h0, h0, p->p);
        if ( mpz_tstbit(n,t) ) {
            // comput x*(h1*x+h0) mod (x^2-bx+a)
            mpz_mul (w, h1, a);  mpz_mul_ui (x, h1, b);  mpz_add (h1, h0, x);  mpz_mod (h1, h1, p->p);
            mpz_neg (h0, w);  mpz_mod (h0, h0, p->p);
        }
    }
    if ( mpz_sgn (h1) ) {
        sts = -1;
    } else {
        mpz_mul_mod (x, h0, h0, p->p);
        sts = ( mpz_cmp (x, a) == 0 ? 1 : -1 );
    }
    if ( sts == 1 ) mpz_set (c, h0);
    mpz_clear (n); mpz_clear (h0); mpz_clear (h1); mpz_clear (w);  mpz_clear (x);
    if ( sts < 0 ) gmp_printf ("Non-prime %Zd detected in mpz_sqrt_mod2\n", p->p);
    return sts;
}


// c = sqrt(a) mod p.  returns 1 if sqrt is found, 0 if no sqrt exists, -1 if p is found to be composite
// note that we might find a sqrt even when p is composite, and in this case we just return the square root, we are not obliged to prove that p is prime
// Algorithm 1.5.1 in [CANT] (or Algorithm 11.23 in [handbook of HECC]) - standard Tonelli-Shanks algorithm, modified to recognize composite p
// c = a is OK
int zp_sqrt (zp_t c, zp_t a, zp_mod_t p)
{
    zp_t b, t, x, y;
    int  i, r, m, sts;

    if ( zp_is_zero (a) ) { zp_set_zero (c);  return 1; }
    assert (mpz_sgn(a) > 0 && mpz_cmp(a, p->p) < 0);
    // revert to alternative algorithm if the 2-Sylow is huge (or if p==3 mod 4, in which case life is easy)
    if ( p->k == 1 || p->k*p->k > zp_bits(p) ) return zp_sqrt_muller (c, a, p);
    if ( mpz_legendre (a, p->p) < 0 ) return 0;                     // for large p this is well worth doing, and we need to do this if p is not necessarily prime

    zp_init (b, p); zp_init (t, p); zp_init (x, p); zp_init (y, p); 
    zp_exp (x, a, p->o, p);                                     // x = a^((m-1)/2) where p-1=2^k*m  (p->o precomputed by zp_init_mod)
    zp_mul (t, a, x, p);  zp_mul (b, t, x, p);  zp_set (x, t);              // b = ax^2, x = ax, note that x^2 = ab
    zp_set (y, p->g2);                                          // generator p->g2 for the 2-sylow precomputed by zp_init_mod
    r = p->k;                                                   // 2-adic valuation of p-1
    while ( ! zp_is_one (b) ) {
        for ( m = 1, zp_sqr (t, b, p) ; m <= r && ! zp_is_one(t) ; m++ ) zp_sqr (t, t, p);
        if ( m == r ) { sts = 0; goto done; }                       // this can happen only if a is not a residue (which is impossible since we check the legendre symbol above, but we'll leave it in case this changes)
        if ( m > r ) { sts = -1; goto done; }                       // this is impossible if p is actually prime
        zp_set (t, y); for ( i = 0 ; i < r-m-1 ; i++ ) zp_sqr (t, t, p);        // t = y^(r-m-1)
        zp_sqr (y, t, p);  r = m;
        zp_mul (x, x, t, p);  zp_mul (b, b, y, p);                      // x = tx, b = t^2b so we still have x^2=ab
    }
    if ( ! zp_is_one(b) ) { i = -1; goto done; }    // impossible if p is prime, we checked that a is a residue above
    // verify result (only necessary in the case that p isn't actually prime)
    zp_sqr (y, x, p);  if ( zp_equal (y, a) ) { zp_set (c, x);  sts = 1; } else { zp_set_zero(c); sts = -1; }
done:
    zp_clear (b); zp_clear (t); zp_clear (x); zp_clear (y); 
    if ( sts == -1 ) gmp_printf ("zp_sqrt detected non-prime p=%Zd\n", p->p);
    return sts;
}

// randomized algorithm to find a generator g for the r-Sylow subgroup of (Z/pZ)*
// also computes a primitive rth root of unity z, if one exists in (Z/pZ)* (ow z is set to zero)
// returns the r-adic valuation of order of the r-Sylow subgroup
int zp_sylow_gen (zp_t g, zp_t z, int r, zp_mod_t p)
{
    mpz_t n, x;
    int k;

    if ( r == 2 ) { zp_set (g, p->g2);  zp_set_one(z); zp_neg(z, z, p); return p->k; }
    mpz_init_modp (n, p->p);  mpz_init_set_ui (x,r);
    k = mpz_remove (n, p->pm1, x);                      // make p-1 prime to r by removing all factors of  r
    if ( ! k ) { zp_set_one(g); zp_set_zero(z); goto done; }        // if n is not divisible by r than the r-Sylow subgroup is trivial
    mpz_ui_pow_ui (x, r, k-1);                          // x = r^(k-1)
    do {
        do { zp_random (z, p); } while ( zp_is_zero(z) );
        zp_exp (g, z, n, p);  zp_exp (z, g, x, p);
    } while ( zp_is_one (z) );
done:
    mpz_clear (n);  mpz_clear (x);
    return k;
}


// attempts to compute an rth root of a, where r is a small prime (e.g. 3), given a generator g for the r-Sylow subgroup of (Z/pZ)*
// for p = 1 mod r, we use AMM with a naive dlog but organize the exponentiations carefully so that the cost should be roughly equal to a single exponentiation
int zp_rth_root (zp_t b, zp_t a, int r, zp_t g, zp_mod_t p)
{
    zp_t c, f, w;
    mpz_t e, n, t, x, y;
    int i, ei, k, sts;

    if ( zp_is_zero (a) ) { zp_set_zero(b); return 1; }
    
    mpz_init (n);
    mpz_sub_ui (n, p->p, 1);                                // n = p-1 is the order of (Z/pZ)*
    if ( ! mpz_divisible_ui_p (n, r) ) {
        // if the r-Sylow is trivial, then the rth root of a is just a^t where rx =1 mod n
        mpz_init (x);  mpz_set_ui(x,r);
        mpz_invert (x, x, n);
        zp_exp (b, a, x, p);
        mpz_clear (x);  mpz_clear (n);
        return 1;
    }
    
    zp_init (c, p); zp_init  (f, p); zp_init  (w, p);
    mpz_init (e); mpz_init(t); mpz_init (x); mpz_init (y);          

    // in order to reduce the amount of exponentiating, our plan is to first compute b=a^((ut+1)/r), then use this to quickly get c = b^r/a = a^(ut)
    // we will then compute e=dlog(g,a^(-ut)) so that g^(e/r)*a^((ut+1)/r) = a^(1/r)
    mpz_set_ui(e,r);
    k = mpz_remove (t, n, e);                           // n = (r^s)t where (r,t)=1
    mpz_invert (x, t, e);  mpz_sub (x, e, x);                   // u = -1/t mod r stored in x
    mpz_mul (t, t, x);  mpz_add_ui (t, t, 1);
    mpz_divexact_ui (t, t, r);
    zp_invert (w, a, p);                                // invert a here in case b and a coincide
    zp_exp (b, a, t, p);                                // b = a^((ut+1)/r)
    zp_exp_ui (c, b, r, p);
    zp_mul (c, c, w, p);                                // c = a^(ut)
    
    mpz_ui_pow_ui (x, r, k-1); mpz_set_one(y);              // x=r^(k-1), y =r^0
    zp_set (f, g);  mpz_set_zero (e);                       // f=g^(r^0), c = g^e*a^(ut)
    // compute e = dlog(g,a^-(ut)) so that g^e*a^(ut) = 1
    // use naive digit by digit approach -- we expect k=1 on average
    for ( i = 0 ; i < k ; i++ ) {                           // i ranges over the digits of e
        if ( i ) {
            mpz_divexact_ui (x, x, r);                      // x = r^(k-1-i)
            mpz_mul_ui (y, y, r);                       // y = r^i
            zp_exp_ui (f, f, r, p);                         // f = g^(r^i)
        }
        // compute ith digit of e
        for ( ei = 0 ; ei < r ; ei++ ) {
            zp_exp (w, c, x, p);  if ( zp_is_one (w) ) break;
            mpz_add (e, e, y);                          // e  = e + e_ir^i
            zp_mul (c, c, f, p);                        // c = g^e*a^(ut)
        }
        assert ( ei < r );
    }
    if ( ! mpz_divisible_ui_p (e, r) ) { sts = 0; goto done; }
    mpz_divexact_ui (e, e, r);
    zp_exp (w, g, e, p);                                // w = g^(e/r) = a^(-ut/r)
    zp_mul (b, b, w, p);                                // b = a^((ut+1)/r) * a^(-ut/r) = a^(1/r)
    sts = 1;
done:
    zp_clear (c); zp_clear (f); zp_clear (w);
    mpz_clear (x); mpz_clear (y); mpz_clear (e); mpz_clear (t); mpz_clear (n);
    
    return sts;
}


/*
    sqrt algorithm over F_p^2, reduces problem to two sqrts in F_p (and an inversion)

    We are given an element of the form a1z+a0 where z^2=s (a non-residue in F_p)
    If a1 is zero, we just compute sqrt(a0) in F_p, and if it doesn't exist, sqrt(a0)=sqrt(a0/s)*sqrt(s) = sqrt(a0/s)*z.

    We now assume a1!=0.  If (b1z+b0)^2 = (a1z+a0) then

       b0^2+sb1^2 = a0  and  2b0b1=a1

    and we know that both b0 and b1 are non-zero.  We then obtain

       b0^4 - a0b0^2 + a1^2s/4 = 0 =>  b0^2 = (a0 +/- sqrt(N(a)))/2
 
    and similarly

      sb1^4 - a0b1^2 + a1^2/4 = 0 => b1^2 = (a0 -/+ sqrt(N(a)))/(2s)

    Note that b0^2*b1^2 = a1^2/4 is a QR, but 1/s is not a QR, so exactly one of (a0+sqrt(N(a)))/2 and (a0-sqrt(N(a)))/2 is a QR
    
    We could speed up this code using an invsqrt operation (as in ff.c)
*/
int zp2_sqrt (zp2_t b, zp2_t a, zp_mod_t p)
{
    zp_t t0, t1;
    
    if ( zp2_is_in_zp(a) ) {
        if ( zp_sqrt (b[0], a[0], p) ) {                // b0 = sqrt(a0) in Fp
            zp_set_zero (b[1]);
        } else {
            zp_mul (b[0], p->si, a[0], p);
            if ( ! zp_sqrt (b[1], b[0], p) )            // sqrt(a0) = sqrt(a0/s)*z
                { gmp_printf ("Fatal error in zp2_sqrt, s=%Zd, p=%Zd, a=%Zd*z+%Zdn",p->s,p,a[1],a[0]); abort(); }
            zp_set_zero (b[0]);
        }
        return 1;
    }
    zp_init (t0, p); zp_init (t1, p);
    zp2_norm(t0, a, p);
    if ( ! zp_sqrt (t1, t0, p) )
        { zp_clear (t0); zp_clear (t1); return 0; } // a is a QR in F_p^2 iff N(a) is a QR in F_p
    zp_add (t0, t1, a[0], p);
    zp_div2 (t1, t0, p);                        // t1 = (a0 + sqrt(N(a))) / 2
    if ( ! zp_sqrt (b[0], t1, p) ) {                    // compute b0 as sqrt(a0 +/- sqrt(N(a))
        zp_sub (t1, t1, a[0], p);
        zp_neg (t1, t1, p);
        if ( ! zp_sqrt (b[0], t1, p) )              // this must work.  Note that mpz_sqrt_mod applies a Jacobi test first, so we really only incur the full cost of one square root
            { gmp_printf ("Fatal error in zp2_sqrt, s=%Zd, p=%Zd, a=%Zd*z+%Zdn",p->s,p,a[1],a[0]); abort(); }
    }
    zp_add (t1, b[0], b[0], p);
    zp_invert (t1, t1, p);
    zp_mul (b[1], t1, a[1], p);
    zp_clear (t0); zp_clear (t1);
    return 1;
}


// randomized algorithm to find a generator for the r-Sylow subgroup of Fp2*
int zp2_sylow_gen (zp2_t g, zp2_t z, int r, zp_mod_t p)
{
    mpz_t n, x;
    int k;

    assert ( r&1 ); // there is no need to support r==2 given zp2_sqrt above
    
    // if r-Sylow is contained in Fp*, take advantage of this
    if ( mpz_congruent_ui_p (p->p, 1, r) ) { mpz_set_zero (g[1]); mpz_set_zero(z[1]); return zp_sylow_gen (g[0], z[0], r, p); }

    // Ok, we now know r does not divide p-1, so (p^2-1)/r = (p-1)(p+1)/r
    // the r-Sylow subgroup lies inside the subgroup of elements with norm 1 (the image of the [p-1]-power map)
    // so it suffices to generate random elements of norm 1 and then exponentiate by (p+1)/r
    
    mpz_init2 (n, 2*zp_bits(p));  mpz_init (x);
    mpz_add_ui (n, p->p, 1);
    mpz_set_ui(x,r);
    k = mpz_remove (n, n, x);                                   // make n prime to r by removing all factors of r (which is assumed to be prime)
    if ( ! k ) { zp2_set_one (g); zp2_set_zero(z); goto done; }
    mpz_ui_pow_ui (x, r, k-1);                                  // x = r^(k-1)
    do {
        do { zp2_random (z, p); } while ( zp2_is_zero(z) );             // get  a random element z of Fp2*
        zp2_norm(g[0], z, p); zp_invert (g[0], g[0], p);                // compute g[0] = 1/N(z)
        zp2_sqr (z, z, p); zp2_mul_zp (z, z, g[0], p);                  // replace  z by z^2/N(z) with norm 1
        zp2_exp (g, z, n, p);                                   // g  is now in the r-Sylow
        zp2_exp (z, g, x, p);                                       // test whether g has maximal order r^k
    } while ( zp2_is_one(z) );
done:
    mpz_clear (n);  mpz_clear (x);
    return k;
}


// computes an rth root b of a in Fp[x]/(x^2-d) given a generator g of the r-Sylow subgroup (r must be prime)
// returns 0 if no rth root of a exists, 1 ow.  a and b may coincide
// for p = 1 mod r, we use AMM with a naive dlog but order the exponentiations carefully so that the total cost is about equal to a single exponentiation
// this is essentially an exact copy of zp_rth_root
int zp2_rth_root (zp2_t b, zp2_t a, int r, zp2_t g, zp_mod_t p)
{
    zp2_t c, f, w, z;
    mpz_t e, n, t, x, y;
    int i, ei, k, sts;

    if ( zp2_is_zero(a) ) { zp2_set_zero(b); return 1; }
    
    // if a is actually in Fp, take advantage of this
    if ( zp2_is_in_zp(a) ) {
        zp_set_zero(b[1]);
        // either the r-sylow of Fp2* lies in Fp*, or the r-Sylow of Fp* is trivial -- either way its safe to specify g[0] as a generator (it will be ignored in the trivial case)
        return zp_rth_root (b[0], a[0], r, g[0], p);
    }
    mpz_init2 (n, 2*mpz_bits(p->p));    
    mpz_mul (n, p->p, p->p);  mpz_sub_ui (n, n, 1);             // n = p^2-1 is the order of Fp2*

    if ( ! mpz_divisible_ui_p (n, r) ) {
        // if the r-Sylow is trivial, then the rth root of a is just a^t where rx =1 mod n
        mpz_init2 (x, 2*mpz_bits(p->p));  mpz_set_ui(x, r);
        mpz_invert (x, x, n);
        zp2_exp (b, a, x, p);
        mpz_clear (x);  mpz_clear (n);
        return 1;
    }
    zp2_init (c, p); zp2_init (f, p); zp2_init (w, p); zp2_init (z, p);
    mpz_init (e); mpz_init(t); mpz_init (x); mpz_init (y);

    // in order to reduce the amount of exponentiating, our plan is to first compute b=a^((ut+1)/r), where x = -1/t mod r, then use this to quickly get c = b^r/a = a^(ut)
    // we will then compute e=dlog(g,a^(-ut)) so that g^(e/r)*a^((ut+1)/r) = a^(1/r)
    mpz_set_ui(e,r);
    k = mpz_remove (t, n, e);                           // n = (r^s)t where (r,t)=1
    mpz_invert (x, t, e);  mpz_sub (x, e, x);                   // u = -1/t mod r stored in x
    mpz_mul (t, t, x);  mpz_add_ui (t, t, 1);
    mpz_divexact_ui (t, t, r);
    zp2_invert (w, a, p);                                   // invert a here (since b and a may coincide!)
    zp2_exp (b, a, t, p);                                   // b = a^((ut+1)/r)
    zp2_exp_ui (c, b, r, p);
    zp2_mul (c, c, w, p);                               // c = a^(ut)

    mpz_ui_pow_ui (x, r, k-1); mpz_set_one(y);              // x=r^(k-1), y =r^0
    zp2_set (f, g);  mpz_set_zero (e);                      // f=g^(r^0), c = g^e*a^(ut)
    // compute e = dlog(g,a^-t) so that g^e*a^t = 1
    // use naive digit by digit approach -- we expect k=1 on average
//gmp_printf ("r=%d, k=%d, x=%Zd, y=%Zd, f:=%Zd*z+%Zd, c=%Zd*z+%Zd\n", r, k, x, y, f[1], f[0], c[1], c[0]);
    for ( i = 0 ; i < k ; i++ ) {                           // i ranges over the digits of e
        if ( i ) {
            mpz_divexact_ui (x, x, r);                      // x = r^(k-1-i)
            mpz_mul_ui (y, y, r);                       // y = r^i
            zp2_exp_ui (f, f, r, p);                        // f = g^(r^i)
        }
        // compute ith digit of e
        for ( ei = 0 ; ei < r ; ei++ ) {
            zp2_exp (z, c, x, p);  if ( zp2_is_one(z) ) break;
            mpz_add (e, e, y);                          // e  = e + e_ir^i
            zp2_mul (c, c, f, p);                           // c = g^e*a^t
        }
        assert ( ei < r );
    }
    if ( ! mpz_divisible_ui_p (e, r) ) { sts = 0; goto done; }
    mpz_divexact_ui (e, e, r);
    zp2_exp (w, g, e, p);                               // w = g^(e/r) = a^(-ut/r)
    zp2_mul (b, b, w, p);                               // b = a^((ut+1)/r) * a^(-ut/r) = a^(1/r)
    sts = 1;
done:
    mpz_clear (e);  mpz_clear (n);mpz_clear (t); mpz_clear (x); mpz_clear (y); 
    zp2_clear (c); zp2_clear (f); zp2_clear (w); zp2_clear (z); 
    return sts;
}

static inline int zp2_cbrt (zp2_t b, zp2_t a, zp_mod_t p)
{
    if ( zp2_is_zero (p->g3) ) zp2_sylow_gen (p->g3, p->z3, 3, p);
    return zp2_rth_root (b, a, 3, p->g3, p);
}


// finds roots of a monic quadratic, returns number of roots found (0 or 2)
// f and r may overlap
int zp_poly_quadratic_roots (zp_t r[2], zp_t f[2], zp_mod_t p)
{
    zp_t w;
    
    mpz_init_modp_w (w, p->p);
    // w = f[1]^2 - 4*f[0]
    mpz_mul_2exp (w, f[0], 2);  mpz_submul (w, f[1], f[1]);  mpz_neg (w, w);  mpz_mod (w, w, p->p);
    if ( ! zp_sqrt (w, w, p) ) { mpz_clear (w); return 0; }
    zp_sub (r[0], w, f[1], p);
    zp_div2 (r[0], r[0], p);
    zp_sub (r[1], r[0], w, p);
    mpz_clear (w);
    return 2;
}


// finds roots of a monic quadratic, returns number of roots found (0 or 2)
// f and r may overlap
int zp2_poly_quadratic_roots (zp2_t r[2], zp2_t f[2], zp_mod_t p)
{
    zp2_t w;
    
    zp2_init (w, p);
    // w = f[1]^2 - 4*f[0]
    zp2_add (w, f[0], f[0], p); zp2_add (w, w, w, p); 
    zp2_sqr (r[0], f[1], p); zp2_sub (w, r[0], w, p);
    if ( ! zp2_sqrt (w, w, p) ) { zp2_clear (w); return 0; }
    zp2_sub (r[0], w, f[1], p);
    zp2_div2 (r[0], r[0], p);
    zp2_sub (r[1], r[0], w, p);
    zp2_clear (w);
    return 2;
}


// finds roots of x^3 + ax + b
// handles cases where r[0] or r[1] overlaps a or b
int zp_poly_depressed_cubic_roots (zp_t r[3], zp_t a, zp_t b, zp_mod_t p)
{
    zp2_t d, g, u, v;
    zp_t t;
    int i, j, k, m;

//gmp_printf ("solving x^3+%Zd*x+%Zd = 0 mod %Zd\n", a, b, p);
    
    // handle x^3+ax case separately
    if ( zp_is_zero (b) ) {
        if ( zp_is_zero(a) ) { for ( i = 0 ; i < 3 ; i++ ) zp_set_zero (r[i]);  return 3; }
        zp_set_zero (r[0]);
        zp_neg(r[1], a, p);
        if ( zp_sqrt (r[1], r[1], p) ) { zp_neg (r[2], r[1], p); return 3; } else return 1;
    }
    
    // handle x^3+b case separately
    if ( zp_is_zero (a) ) {
        if ( zp_is_zero(b) ) { for ( i = 0 ; i < 3 ; i++ ) zp_set_zero (r[i]);  return 3; }
        if ( ! zp_cbrt (r[0], b, p) ) return 0;
        zp_neg (r[0], r[0], p);
        if ( ! p->p1mod3 ) return 1;    // no cube root of unity when p=2 mod 3, so only one cube root
        zp_mul (r[1], r[0], p->z3[0], p);  zp_mul (r[2], r[1], p->z3[0], p);
        return 3;
    }
    zp_init (t, p); zp2_init (d, p);
    // compute t = b^2/4 + a^3/27 = disc(x^3+ax+b)/(-3)  (use d[0] and d[1] as temps)
    zp_div2 (d[0], b, p);  zp_sqr (d[0], d[0], p);  zp_div3 (d[1], a, p);  zp_sqr (t, d[1], p);  zp_mul (d[1], d[1], t, p);  zp_add (t, d[0], d[1], p);
    // check for double root -- occurs iff t == 0
    if ( zp_is_zero (t) ) {
        // double root must be r = -3b/2a (we know that a is nonzero)
        zp_invert (t, a, p);  zp_div2 (t, t, p);
        zp_mul_ui (t, t, 3, p);  zp_mul (t, t, b, p);
        zp_neg (t, t, p);
        // x^3 + ax + b = (x-r)(x^2+rx+2a/3)
        zp_div3 (r[1], a, p);  zp_add (r[1], r[1], r[1], p);  zp_set (r[2], t);  zp_set (r[0], t);
        zp_clear (t);  zp2_clear (d);
        if ( zp_poly_quadratic_roots (r+1, r+1, p) ) return 3; else return 1;
    }
    // We have distinct roots and we know a and b are nonzero, so we apply Cardona's method.
    // for simplicity we work in Fp2, even though this isn't necessary in every case
    // most of the zp2 functions are optimized to notice when their inputs lie in Fp, so this doesn't slow things down much
    
    zp2_init (u, p);  zp2_init (v, p);  zp2_init (g, p); 
    
    zp_set (u[0], t);  zp_set_zero (u[1]);
    if ( ! zp2_sqrt (d, u, p) ) { gmp_printf ("Fatal error, square root of %Zd mod %Zd not found in quadratic extension!\n", d[0], p); abort(); }
    // compute maximum number k of possible roots via Stickelberger parity theorem
    // we know that the disc and t have the same quad char if p=1 mod 3 and opposite quad char if p = 2 mod 3
    if ( zp2_is_in_zp (d) ) {
        k = ( p->p1mod3 ? 3 : 1 );
    } else {
        k = ( p->p1mod3 ? 1 : 3 );
    }
    zp_div2 (t, b, p);
    zp2_set (u, d);  zp_sub (u[0], u[0], t, p);
    if ( ! zp2_cbrt (u, u, p) ) { k = 0;  goto done; }
    zp2_mul_ui (v, u, 3, p);  zp2_invert (v, v, p);  zp2_mul_zp (v, v, a, p);  zp2_neg (v, v, p);       // v = -a/3u
    // test z3^i u+z3^j*v for i,j=0,1,2 to find roots in Fp (TODO: understand why this seems to be necessary and why it works!)
    // relocate a and b in case they overlap with r
    m = 0; zp_set (d[0], b);  zp_set (d[1], a);
    for ( i = 0 ; i < 3 ; i++ ) {
        for ( j = 0 ; j < 3 ; j++ ) {
            zp2_add (g, u, v, p);
            if ( zp2_is_in_zp (g) ) {
                zp_sqr (t, g[0], p);
                zp_add (t, t, d[1], p);
                zp_mul (t, t, g[0], p);
                zp_add (t, t, d[0], p);
                if ( zp_is_zero (t) ) { zp_set (r[m++], g[0]); if ( m==k ) goto done; }
            }
            zp2_mul (v, v, p->z3, p);
        }
        zp2_mul (u, u, p->z3, p);
    }
    if ( k == 1 || m ) { gmp_printf ("Fatal error in zp_poly_cubic_roots, only found %d of %d roots of x^3+%Zd*x+%Zd mod %Zd\n", m, k, a, b, p); abort(); }
    k = 0;
done:
    zp2_clear (d);  zp2_clear (g); zp2_clear (u);  zp2_clear (v); mpz_clear (t);
    return k;
}

// given f(x)=x^3+f[2]x^2+f[1]x+f[0] computes a, b, and s such that f(x-s) = x^3+ax+b
// a and b cannot overlap f!
void zp_poly_depress_monic_cubic (zp_t a, zp_t b, zp_t s, zp_t f[3], zp_mod_t p)
{
    zp_div3 (s, f[2], p);                   // store shift in s
    mpz_mul_ui (a, s, 3);
    mpz_mul_2exp (b, f[2], 1);
    mpz_sub (a, a, b);
    zp_mul (a, a, s, p);
    zp_add (a, a, f[1], p);             // a = f[1] + s(3s-2f[2])
    mpz_sub (b, f[2], s);
    zp_mul (b, b, s, p);
    mpz_sub (b, f[1], b);
    zp_mul (b, b, s, p);
    zp_sub (b, f[0], b, p);             // b = f[0] - s(f[1]-s(f[2]-s))
}

// finds roots of a monic cubic, returns number of roots found (0, 1, or 3).  r and f may coincide
int zp_poly_cubic_roots (zp_t r[3], zp_t f[3], zp_mod_t p)
{
    zp_t s, a, b;
    int i, n;

    if ( zp_is_zero (f[2]) ) return  zp_poly_depressed_cubic_roots (r, f[1], f[0], p);
        
    zp_init (s, p); zp_init (a, p); zp_init (b, p);
    zp_poly_depress_monic_cubic (a, b, s, f, p);
    n = zp_poly_depressed_cubic_roots (r, a, b, p);
    for ( i = 0 ; i < n ; i++ ) zp_sub (r[i], r[i], s, p);
    zp_clear (s); zp_clear (a); zp_clear (b);
    return n;
}



// finds roots of x^3 + ax + b
// handles cases where r[0] or r[1] overlaps a or b
int zp2_poly_depressed_cubic_roots (zp2_t r[3], zp2_t a, zp2_t b, zp_mod_t p)
{
    zp2_t d, g, t, u, v, w;
    int i, k;

//gmp_printf ("solving x^3+(%Zd*z+%Zd)*x+(%Zd*z+%Zd) = 0 mod %Zd (z^2 = %ld)\n", a[1], a[0], b[1], b[0], p->p, p->s);
    
    
    // handle x^3+ax case separately
    if ( zp2_is_zero (b) ) {
        if ( zp2_is_zero(a) ) { for ( i = 0 ; i < 3 ; i++ ) zp2_set_zero (r[i]);  return 3; }
        zp2_set_zero (r[0]);
        zp2_neg(r[1], a, p);
        if ( zp2_sqrt (r[1], r[1], p) ) { zp2_neg (r[2], r[1], p); return 3; } else return 1;
    }

    // handle x^3+b case separately
    if ( zp2_is_zero (a) ) {
        if ( zp2_is_zero(b) ) { for ( i = 0 ; i < 3 ; i++ ) zp2_set_zero (r[i]);  return 3; }
        if ( ! zp2_cbrt (r[0], b, p) ) return 0;
        zp2_neg (r[0], r[0], p);
        zp2_mul (r[1], r[0], p->z3, p);
        zp2_mul (r[2], r[1], p->z3, p);
        return 3;
    }
    zp2_init (t, p); zp2_init (d, p);
    // compute t = b^2/4 + a^3/27 = disc(x^3+ax+b)/(-3)
    zp2_div3 (d, a, p);  zp2_sqr (t, d, p);  zp2_mul (d, d, t, p); zp2_div2 (t, b, p);  zp2_sqr (t, t, p); zp2_add (d, d, t, p);
    
    // check for double root -- occurs iff t == 0
    if ( zp2_is_zero (d) ) {
        // double root must be r = -3b/2a (we know that a is nonzero)
        zp2_invert (t, a, p);  zp2_div2 (t, t, p); zp2_add (d, t, t, p); zp2_add (t, t, d, p);  zp2_mul (t, t, b, p);  zp2_neg (t, t, p);
        // x^3 + ax + b = (x-r)(x^2+rx+2a/3)
        zp2_div3 (r[1], a, p);  zp2_add (r[1], r[1], r[1], p);  zp2_set (r[2], t);  zp2_set (r[0], t);
        zp2_clear (t);  zp2_clear (d);
        if ( zp2_poly_quadratic_roots (r+1, r+1, p) ) return 3; else return 1;
    }
    
    // We have distinct roots and we know a and b are nonzero, so we apply Cardona's method.
    
    zp2_init (u, p);  zp2_init (v, p);  zp2_init (g, p); zp2_init (w, p); 
    
    if ( ! zp2_sqrt (d, d, p) ) { zp2_set_zero (r[0]); k=1; goto done; puts ("single root case in zp2_poly_depressed_cubic_roots not implemented yet."); abort(); } // TODO: compute the root using gcd(X^p-X,f) by calling zp_poly_linear_factor_product

    // Over Fp^2, -3 is always a QR, so d is a QR iff disc(x^3+ax+b) is a QR iff f has 0 or 3 roots (by Stickelberger)
    k = 0;

    zp2_div2 (t, b, p);
    zp2_sub (u, d, t, p);
    if ( ! zp2_cbrt (u, u, p) ) goto done;
    zp2_add (v, u, u, p); zp2_add (v, v, u, p); zp2_invert (v, v, p);  zp2_mul (v, v, a, p);  zp2_neg (v, v, p);            // v = -a/3u

    // test z3^i u+v for i=0,1,2 to find a root
    // relocate a and b in case they overlap with r
    zp2_set (d, a);  zp2_set (w, b);
    for ( i = 0 ; i < 3 ; i++ ) {
        zp2_add (g, u, v, p);                                                               // test root g = u+v
        zp2_sqr (t, g, p);  zp2_add (t, t, d, p);  zp2_mul (t, t, g, p); zp2_add (t, t, w, p);                  // t = f(g) = g^3+a*g+b
        if ( zp2_is_zero (t) ) { zp2_set (r[k++], g); break; }
        zp2_mul (u, u, p->z3, p);
    }
    if ( k ) {
        // the other 2 roots are z3^2u + z3v and z3u + z3^2v
        zp2_mul (u, u, p->z3, p);  zp2_mul (v, v, p->z3, p);
        zp2_mul (g, u, p->z3, p);  zp2_add (g, g, v, p);
        zp2_sqr (t, g, p);  zp2_add (t, t, d, p);  zp2_mul (t, t, g, p); zp2_add (t, t, w, p);                  // t = f(g) = g^3+a*g+b
        if ( zp2_is_zero (t) ) { zp2_set (r[k++], g); } else { gmp_printf ("Fatal error in zp_poly_cubic_roots, only found %d of 3 roots of x^3+%Zd*x+%Zd mod %Zd\n", k, a, b, p); abort(); }
        zp2_mul (g, v, p->z3, p);  zp2_add (g, g, u, p);
        zp2_sqr (t, g, p);  zp2_add (t, t, d, p);  zp2_mul (t, t, g, p); zp2_add (t, t, w, p);                  // t = f(g) = g^3+a*g+b
        if ( zp2_is_zero (t) ) { zp2_set (r[k++], g); } else { gmp_printf ("Fatal error in zp_poly_cubic_roots, only found %d of 3 roots of x^3+%Zd*x+%Zd mod %Zd\n", k, a, b, p); abort(); }
    }
done:
//if ( k ) gmp_printf ("(%Zd*z+%Zd), (%Zd*z+%Zd), (%Zd*z+%Zd)\n", r[0][1], r[0][0], r[1][1], r[1][0], r[2][1], r[2][0]);
    zp2_clear (d);  zp2_clear (g); zp2_clear (w); zp2_clear (u);  zp2_clear (v); zp2_clear (t);
    return k;
}

// given f(x)=x^3+f[2]x^2+f[1]x+f[0] computes a, b, and s such that f(x-s) = x^3+ax+b
void zp2_poly_depress_monic_cubic (zp2_t a, zp2_t b, zp2_t s, zp2_t f[3], zp_mod_t p)
{
//gmp_printf ("Depressing cubic x^3 + (%Zd*z+%Zd)*x^2 + (%Zd*z+%Zd)*x + (%Zd*z+%Zd)\n", f[2][1], f[2][0], f[1][1], f[1][0], f[0][1], f[0][0]);
    zp2_div3 (s, f[2], p);                  // store shift in s
    zp2_add (a, s, s, p);  zp2_add (a, a, s, p);    // a = 3s
    zp2_add (b, f[2], f[2], p);                 // b = 2f2
    zp2_sub (a, a, b, p);
    zp2_mul (a, a, s, p);
    zp2_add (a, a, f[1], p);                    // a = f[1] + s(3s-2f[2])
    zp2_sub (b, f[2], s, p);
    zp2_mul (b, b, s, p);
    zp2_sub (b, f[1], b, p);
    zp2_mul (b, b, s, p);
    zp2_sub (b, f[0], b, p);                    // b = f[0] - s(f[1]-s(f[2]-s))
//gmp_printf ("a=%Zd*z+%Zd, b=%Zd*z+%Zd\n", a[1], a[0], b[1], b[0]);
}

// finds roots of a monic cubic, returns number of roots found (0, 1, or 3), r and f may coincide 
int zp2_poly_cubic_roots (zp2_t r[3], zp2_t f[3], zp_mod_t p)
{
    zp_t g[4], rr[3];
    zp2_t s, a, b;
    int i, k;
    
    // if f has any roots over Fp^2, at least one them must be defined over Fp (because Fp^2 is a degree 2 extension and a cubic without a root is irreducible)
    if ( zp2_is_in_zp(f[2]) &&  zp2_is_in_zp(f[1]) &&  zp2_is_in_zp(f[0]) ) {
        for ( i = 0 ; i < 3 ; i++ ) { mpz_init_set (g[i], f[i][0]); zp_init (rr[i], p); }
        k = zp_poly_cubic_roots (rr, g, p);
        if ( k == 1 ) {
            zp2_set_zp (r[0], rr[0]);
            mpz_init_set_ui (g[3],1);
            zp_poly_remove_root (g, g, 3, rr[0], p->p);
            mpz_clear (g[3]);
            zp2_set_zp (r[1], g[0]);  zp2_set_zp (r[2], g[1]);
            k += zp2_poly_quadratic_roots (r+1, r+1, p);
        } else if ( k == 3 ) {
            for ( i = 0 ; i < k ; i++ ) zp2_set_zp (r[i], rr[i]);
        }
        for ( i = 0 ; i < 3 ; i++ ) { zp_clear (g[i]); zp_clear (rr[i]); }
        return k;
    }
    if ( zp2_is_zero (f[2]) ) return  zp2_poly_depressed_cubic_roots (r, f[1], f[0], p);
        
    zp2_init (s, p); zp2_init (a, p); zp2_init (b, p);
    zp2_poly_depress_monic_cubic (a, b, s, f, p);
    k = zp2_poly_depressed_cubic_roots (r, a, b, p);
    for ( i = 0 ; i < k ; i++ ) zp2_sub (r[i], r[i], s, p);
    zp2_clear (s); zp2_clear (a); zp2_clear (b);
    return k;
}


// computes the parity of the # of irreducible factors of a monic cubic f via Stickelberger
int zp_poly_cubic_parity (zp_t f[3], zp_mod_t p)
{
    zp_t t1, t2, D;
    int parity;
    
    // compute D = disc(f)
    zp_init (t1, p); zp_init (t2, p); zp_init (D, p);
    zp_sqr (t2, f[1], p);                                                                                           // save t2=f1^2
    mpz_mul (D, t2, f[1]);  mpz_mul_2exp (D, D, 2);                                                                 // D = 4f1^3
    mpz_mul (t1, f[0], f[0]); mpz_mul_ui (t1, t1, 27); mpz_add (D, D, t1); mpz_mod (D, D, p->p);                                // D = 4f1^3 + 27f0^2
    zp_neg (D, D, p);                                                                                           // D = - 4f1^3 - 27f0^2
    if ( ! zp_is_zero (f[2]) ) {
        mpz_mul (t1, f[0], f[2]);  mpz_mul_2exp (t1, t1, 2); mpz_sub (t2, t2, t1); mpz_mod (t2, t2, p->p);                      // t2 = -4f0f2 + f1^2
        mpz_mul (t1, f[0], f[1]);  mpz_mul_ui (t1, t1, 18);  mpz_mul (t2, t2, f[2]); mpz_add(t2, t2, t1); mpz_mod (t2, t2, p->p);       // t2 = f2(-4f0f2 + f1^2) + 18f0f1
        zp_mul (t2, t2, f[2], p);                                                                                   // t2 = 18f0f1f2 - 4f0f2^3 +f1^2f2^2
        zp_add (D, D, t2, p);                                                                                   // D = -27*f0^2 + 18*f0*f1*f2 - 4*f0*f2^3 - 4*f1^3 + f1^2*f2^2
    }
    // apply Stickelberger
    parity = ( mpz_legendre (D, p->p) >= 0 ? 1 : 0 );                                                                   // f has an odd number of factors iff D is a QR mod p
    zp_clear (t1); zp_clear (t2); zp_clear (D);
    return parity;
}


// computes the parity of the # of irreducible factors of a monic cubic f via Stickelberger
int zp2_poly_cubic_parity (zp2_t f[3], zp_mod_t p)
{
    zp2_t t1, t2, D;
    int parity;

    if ( zp2_is_in_zp(f[0]) && zp2_is_in_zp(f[1]) && zp2_is_in_zp(f[2]) ) return 1;                                 // if D is in Fp then it is certainly a QR in Fp^2
    
    // compute D = disc(f)
    zp2_init (t1, p); zp2_init (t2, p); zp2_init (D, p);
    zp2_sqr (t2, f[1], p);                                                                          // save t2=f1^2
    zp2_mul (D, t2, f[1], p);  zp2_add(D, D, D, p); zp2_add (D, D, D, p);                                       // D = 4f1^3
    zp2_sqr (t1, f[0], p);  zp2_mul_ui (t1, t1, 27, p);  zp2_add (D, D, t1, p);                                     // D = 4f1^3 + 27f0^2
    zp2_neg (D, D, p);                                                                              // D = - 4f1^3 - 27f0^2
    if ( ! zp2_is_zero (f[2]) ) {
        zp2_mul (t1, f[0], f[2], p);  zp2_add (t1, t1, t1, p);  zp2_add (t1, t1, t1, p); zp2_sub (t2, t2, t1, p);           // t2 = -4f0f2 + f1^2
        zp2_mul (t1, f[0], f[1], p);  zp2_mul_ui (t1, t1, 18, p);  zp2_mul (t2, t2, f[2], p); zp2_add(t2, t2, t1, p);       // t2 = f2(-4f0f2 + f1^2) + 18f0f1
        zp2_mul (t2, t2, f[2], p);                                                                      // t2 = 18f0f1f2 - 4f0f2^3 +f1^2f2^2
        zp2_add (D, D, t2, p);                                                                      // D = -27*f0^2 + 18*f0*f1*f2 - 4*f0*f2^3 - 4*f1^3 + f1^2*f2^2
    }
    // apply Stickelberger
    zp2_norm (t1[0], D, p);                                                                         // D is a QR in Fp^2 iff N(D) is a QR in Fp
    parity = ( mpz_legendre (t1[0], p->p) >= 0 ? 1 : 0 );                                                   // f has an odd number of factors iff D is a QR mod p
    zp2_clear (t1); zp2_clear (t2); zp2_clear (D);
    return parity;
}

#define ROOTS       100

#ifdef ZP_POLY_CUBIC_ROOTS_PERFTEST

int main (int argc, char *argv[])
{
    clock_t start, end;
    zp_mod_t p;
    mpz_t r[3], f[4], s, x, y, e;
    zp2_t u, v, g;
    int b, i, j, k, n, minb, maxb, delta;
    
    if ( argc < 4 ) { printf ("zp_poly_cubic_roots_perftest minpbits maxpbits delta [iterations]\n"); return 0; }
    minb = atoi(argv[1]);  maxb = atoi(argv[2]);  delta = atoi(argv[3]);
    if ( argc > 4 ) n = atoi(argv[4]); else n = 1;
    mp_set_memory_functions (zp_malloc, zp_realloc, zp_free);

    for ( i = 0 ; i <= 3 ; i++ ) mpz_init2 (f[i], maxb); 
    for ( i = 0 ; i < 3 ; i++ ) mpz_init2 (r[i], maxb);
    mpz_init2(e,2*maxb); mpz_init2(x,maxb); mpz_init2(y,maxb); mpz_init(s);
    mpz_init2(u[0],maxb); mpz_init2(u[1],maxb); mpz_init2(v[0],maxb); mpz_init2(v[1],maxb); mpz_init2(g[0],maxb); mpz_init2(g[1],maxb);
    for ( b = minb ; b <= maxb  ; b+= delta ) {
        for ( k = 1; k <= 2; k++ ) {
            zp_random_prime (x, b);
            while ( ! mpz_congruent_ui_p (x, k, 3) ) mpz_nextprime (x, x);
            zp_mod_init (p, x);
            printf ("%d-bit p = %d mod 3, %ld mod 4\n", b, k, mpz_fdiv_ui (p->p, 4));
            zp_randomb (e, b);
            zp_random (x, p);
            start = clock();
            for ( i = 0 ; i < n ; i++ ) zp_invert (y, x, p);
            end = clock();
            printf ("%d iterations of %d-bit zp_invert took %.6f secs\n", n, b, (double)(end-start)/CLOCKS_PER_SEC);
            zp_randomb (e, b);
            start = clock();
            for ( i = 0 ; i < n ; i++ ) zp_exp (y, x, e, p);
            end = clock();
            printf ("%d iterations of %d-bit zp_exp took %.6f secs\n", n, b, (double)(end-start)/CLOCKS_PER_SEC);
            zp2_random (u, p);
            start = clock();
            for ( i = 0 ; i < n ; i++ ) zp2_invert (v, u, p);
            end = clock();
            printf ("%d iterations of %d-bit zp2_invert took %.6f secs\n", n, b, (double)(end-start)/CLOCKS_PER_SEC);
            start = clock();
            for ( i = 0 ; i < n ; i++ ) for ( j = 0 ; j < b ; j++ ) zp2_sqr (v, u, p);
            end = clock();
            printf ("%d iterations of %d calls to %d-bit zp2_sqr took %.6f secs\n", n, b, b, (double)(end-start)/CLOCKS_PER_SEC);
            zp_randomb (e, 2*b);
            start = clock();
            for ( i = 0 ; i < n ; i++ ) zp2_exp (v, u, e, p);
            end = clock();
            printf ("%d iterations of %d-bit zp2_exp took %.6f secs\n", n, b, (double)(end-start)/CLOCKS_PER_SEC);
            start = clock();
            for ( i = 0 ; i < n ; i++ ) zp2_sylow_gen (p->g3, p->z3, 3, p);
            end = clock();
            printf ("%d iterations of %d-bit zp2_sylow_gen(3) took %.6f secs\n", n, b, (double)(end-start)/CLOCKS_PER_SEC);
            start = clock();
            for ( i = 0 ; i < n ; i++ ) zp2_cbrt (v, u, p);
            end = clock();
            printf ("%d iterations of %d-bit zp2_cbrt took %.6f secs\n", n, b, (double)(end-start)/CLOCKS_PER_SEC);
            zp_set_one (f[1]); zp_random(f[0], p);
            zp_random (x, p);  zp_poly_add_root (f, f, 1, x, p->p);
            start = clock();
            for ( i = 0 ; i < n ; i++ ) zp_poly_quadratic_roots (r, f, p);
            end = clock();
            printf ("%d iterations of %d-bit zp_poly_quadratic_roots took %.6f secs\n", n, b, (double)(end-start)/CLOCKS_PER_SEC);
            zp_set_one (f[1]); zp_random(f[0], p);
            zp_random (x, p);  zp_poly_add_root (f, f, 1, x, p->p);
            zp_random (x, p);  zp_poly_add_root (f, f, 2, x, p->p);
            start = clock();
            for ( i = 0 ; i < n ; i++ ) if ( ! zp_poly_cubic_parity (f, p) ) { printf ("parity check failed!\n"); abort(); }
            end = clock();
            printf ("%d iterations of %d-bit zp_poly_cubic_parity took %.6f secs\n", n, b, (double)(end-start)/CLOCKS_PER_SEC);
            start = clock();
            for ( i = 0 ; i < n ; i++ ) zp_poly_cubic_roots (r, f, p);
            end = clock();
            printf ("%d iterations of %d-bit zp_poly_cubic_roots took %.6f secs\n", n, b, (double)(end-start)/CLOCKS_PER_SEC);
            zp_mod_clear (p);
        }
    }
    mpz_clear (e); mpz_clear (x); mpz_clear (y); mpz_clear (s);
    zp2_clear (u); zp2_clear (v); zp2_clear (g);
    for ( i = 0 ; i <= 3; i++ ) mpz_clear (f[i]);
    for ( i = 0 ; i < 3; i++ ) mpz_clear (r[i]);
    puts ("cleaning up...");
    zp_cleanup ();
    puts ("reporting...");
    zp_mem_report (1);
}

#endif

#ifdef ZP_POLY_CUBIC_ROOTS_TEST

int main (int argc, char *argv[])
{
    zp_mod_t p;
    mpz_t t, f[4], r[3], g[4], pp;
    int maxp;
    int i, d, k;
    
    if ( argc < 2 ) { printf ("zp_poly_cubic_roots maxp\n"); return 0; }
    maxp = atoi (argv[1]);
    for ( i = 0 ; i < 3 ; i++ ) { mpz_init (f[i]); mpz_init (r[i]); mpz_init (g[i]); }
    mpz_init_set_ui (f[3], 1); mpz_init (g[3]);
    mpz_init (t);
    mpz_init_set_ui (pp, 5);
    while ( mpz_cmp_ui (pp, maxp) <= 0 ) {
        zp_mod_init (p, pp);
        do {
            k = zp_poly_cubic_roots (r, f, p);
            if ( k != 0 && k != 1 && k != 3 ) { printf ("Invalid return value %d\n", k); abort(); }
            for ( i = 0 ; i < k ; i++ ) {
                zp_poly_eval (t, f, 3, r[i], p);
                if ( ! zp_is_zero(t) ) { gmp_printf ("zp_poly_cubic_roots returned invalid root %Zd for x^3+%Zdx^2+%Zdx+%Zd mod %Zd, evaluated to %Zd\n", r[i], f[2], f[1], f[0], p->p, t); abort(); }
            }
            if ( k < 3 ) {
                zp_poly_copy (g, f, 3);  d = 3;
                for ( i = 0 ; i < k ; i++ ) { zp_poly_remove_root (g, g, d, r[i], p->p);  d--; }
                if ( d > 0 ) {
                    // check for missed roots using brute force to avoid depending on code in zp_poly_exp (which might call zp_roots functions).  Note that p is small so this isn't so bad.
                    zp_set_zero (t);
                    do {
                        zp_poly_eval (r[0], g, d, t, p);
                        if ( zp_is_zero (r[0]) ) { gmp_printf ("zp_poly_cubic_roots missed root %Zd for x^3+%Zdx^2+%Zdx+%Zd mod %Zd\n", t, f[2], f[1], f[0], p->p); abort(); }
                        zp_inc (t, p);
                    } while ( ! zp_is_zero (t) );
                }
            }
        } while ( zp_poly_next (f, 2, p->p) );
        gmp_printf ("Successfully found roots of every cubic mod %Zd\n", p->p);
        zp_mod_clear (p);
        mpz_nextprime (pp, pp);
    }
}

#endif

#ifdef ZP_SQRT_TEST

int main (int argc, char *argv[])
{
    clock_t start, end;
    zp_mod_t p;
    zp_t a[ROOTS], r[ROOTS], t;
    int i, j ,n, min, max, delta;
    
    if ( argc < 2 ) { puts ("zp_sqrt_test minbits maxbits delta\n"); return 0; }
    min = atoi(argv[1]);
    if ( argc > 2 ) max = atoi(argv[2]); else max = min;
    if ( argc > 3 ) delta = atoi(argv[3]); else delta = 1;

    mpz_init (t);
    for ( i = 0 ; i < ROOTS ; i++ ) { mpz_init (a[i]); mpz_init (r[i]); }
    for ( n = min ; n <= max ; n+=delta ) {
        zp_random_prime (t, n);  zp_mod_init (p, t);
        printf ("k=%d\n", p->k);
        for ( i = 0 ; i < ROOTS ; i++ ) zp_random (a[i], p);
        zp_randomm (t, p->p);
        start = clock();
        for ( i = 0 ; i < ROOTS ; i++ ) zp_exp (r[i], a[i], t, p);
        end = clock();
        printf ("%d calls to zp_exp with %d-bit p took %.6f secs\n", ROOTS, n, (double)(end-start)/CLOCKS_PER_SEC);
        start = clock();
        for ( i = j = 0 ; i < ROOTS ; i++ ) if ( zp_sqrt (r[i], a[i], p) < 1 ) { mpz_set_ui (r[i],0); j++; }
        end = clock();
        for ( i = 0 ; i < ROOTS ; i++ ) {
            if ( ! zp_is_zero(r[i]) ) { zp_sqr (t, r[i], p);  if ( ! zp_equal (t, a[i]) ) {  gmp_printf ("zp_sqrt failed for a=%Zd mod %Zd!\n", a[i], p->p); abort(); } }
            else if ( mpz_legendre (a[i],p->p) > 0 ) {  gmp_printf ("zp_sqrt failed for a=%Zd mod %Zd!\n", a[i], p->p); abort(); }
        }
        printf ("%d calls (%d failed) to zp_sqrt with %d bit p took %.3f secs\n", ROOTS, j, n, (double)(end-start)/CLOCKS_PER_SEC);
        start = clock();
        for ( i = j = 0 ; i < ROOTS ; i++ ) if ( zp_sqrt_muller (r[i], a[i], p) < 1 ) { mpz_set_ui (r[i],0); j++; }
        end = clock();
        for ( i = 0 ; i < ROOTS ; i++ ) {
            if ( ! zp_is_zero(r[i]) ) { zp_sqr (t, r[i], p);  if ( ! zp_equal (t, a[i]) ) {  gmp_printf ("zp_sqrt_muller failed for a=%Zd mod %Zd!\n", a[i], p->p); abort(); } }
            else if ( mpz_legendre (a[i],p->p) > 0 ) {  gmp_printf ("zp_sqrt_muller failed for a=%Zd mod %Zd!\n", a[i], p->p); abort(); }
        }
        printf ("%d calls (%d failed) to zp_sqrt_muller with %d bit p took %.3f secs\n", ROOTS, j, n, (double)(end-start)/CLOCKS_PER_SEC);
        zp_mod_clear (p);
    }
}

#endif

#ifdef ZP_CBRT_TEST

int main (int argc, char *argv[])
{
    clock_t start, end;
    zp_mod_t p;
    mpz_t a[ROOTS], r[ROOTS], t, e1, e2;
    int i, j ,n, min, max, delta;
    
    if ( argc < 2 ) { puts ("zp_cbrt_test minbits maxbits delta\n"); return 0; }
    min = atoi(argv[1]);
    if ( argc > 2 ) max = atoi(argv[2]); else max = min;
    if ( argc > 3 ) delta = atoi(argv[3]); else delta = 1;

    mpz_init (t); mpz_init (e1); mpz_init (e2);
    for ( i = 0 ; i < ROOTS ; i++ ) { mpz_init (a[i]); mpz_init (r[i]); }
    for ( n = min ; n <= max ; n += delta ) {
        zp_random_prime (t, n);  zp_mod_init (p, t);
        mpz_sub_ui (e1, p->p, 1);  mpz_set_ui(r[0], 3);
        j = mpz_remove (e1, e1, r[0]);
        printf ("v_3(p-1) = %d\n", j);
        for ( i = 0 ; i < ROOTS ; i++ ) { zp_random (t, p);  zp_sqr (a[i], t, p);  zp_mul (a[i], a[i], t, p); }     // just test perfect cubes
        zp_randomb (e1, n);  zp_randomb (e2, n);
        start = clock();
        for ( i = 0 ; i < ROOTS ; i++ ) zp_exp (t, a[i], e1, p);
        end = clock();
        printf ("%d zp_exp calls with %d-bit p took %.6f secs\n", ROOTS, n, (double)(end-start)/CLOCKS_PER_SEC);
        start = clock();
        for ( i = 0 ; i < ROOTS ; i++ ) for ( j = 0 ; j < n ; j++ ) zp_sqr (r[i], r[i], p);
        end = clock();
        printf ("%d*%d sqr mods with %d-bit p took %.6f secs\n", ROOTS, n, n, (double)(end-start)/CLOCKS_PER_SEC);
        start = clock();
        for ( i = j = 0 ; i < ROOTS ; i++ ) if ( ! zp_cbrt (r[i], a[i], p) ) { mpz_set_ui (r[i],0); j++; }
        end = clock();
        for ( i = 0 ; i < ROOTS ; i++ ) { zp_sqr (t, r[i], p);  zp_mul (t, t, r[i], p);  if ( ! zp_equal (t, a[i]) ) {  gmp_printf ("mpz_rth_root_mod failed with r=3 for a=%Zd mod %Zd!\n", a[i], p); abort(); } }
        printf ("%d calls (%d failed) to zp_cbrt with %d-bit p took %.6f secs\n", ROOTS, j, n, (double)(end-start)/CLOCKS_PER_SEC);
    }
}

#endif


#ifdef ZP2_SQRT_TEST

int main (int argc, char *argv[])
{
    mpz_t pp;
    zp_mod_t p;
    zp2_t a, b;
    mpz_t w;
    long cnt, maxp;
    
    if ( argc < 2 ) { puts ("zp2_sqrt maxp"); return 0; }
    maxp = atoi(argv[1]);
    mpz_init_set_ui (pp, 5);  mpz_init (w);
    while ( mpz_cmp_ui (pp, maxp) <= 0 ) {
        zp_mod_init (p, pp);
        zp2_init (a, p); zp2_init (b, p);
//      gmp_printf ("z^2=%ld mod %Zd\n", p->s, p->p);
        cnt = 0;
        do {
            do {
                if ( zp2_sqrt (b, a, p) ) {
                    zp2_sqr (b, b, p);
                    if ( ! zp2_equal (a,b) ) { gmp_printf ("Error, invalid sqrt of %Zd*z+%Zd with z^2=%ld mod %Zd\n", a[1], a[0], p->s, p->p); abort(); }
                    cnt++;
                } else {
                }
                zp_inc (a[0], p);
            } while ( ! zp_is_zero(a[0]) );
            zp_inc (a[1], p);
        } while ( ! zp_is_zero (a[1]) );
        mpz_mul(w,p->p,p->p);  mpz_add_ui(w,w,1); mpz_div_2exp (w,w,1);
        if ( mpz_cmp_ui(w,cnt) ) { gmp_printf ("Found %ld of %Zd expected residues\n", cnt, w); abort(); }
        gmp_printf ("Successfully computed square-roots over all of Fp^2 for p=%Zd\n", p->p);
        zp2_clear (a); zp2_clear (b);
        zp_mod_clear (p);
        mpz_nextprime (pp, pp);
    }
    mpz_clear (pp);  mpz_clear (w);
}

#endif


#ifdef ZP2_RTH_ROOT_TEST

int main (int argc, char *argv[])
{
    zp_mod_t p;
    zp2_t a, b, c, g;
    mpz_t w, pp;
    long cnt, maxp;
    int r;
    
    if ( argc < 3 ) { puts ("zp2_rth_root_test r p"); return 0; }
    r = atoi(argv[1]);
    maxp = atoi(argv[2]);
    mpz_init_set_ui (pp, 5);  mpz_init (w);
    while ( mpz_cmp_ui (pp, maxp) <= 0 ) {
        zp_mod_init (p, pp);
//      gmp_printf ("z^2=%ld mod %Zd\n", p->s, p->p);
        zp2_init(a,p); zp2_init (b,p); zp2_init(c,p); zp2_init(g,p);
        zp2_sylow_gen (g, c, r, p);
        cnt = 0;
        do {
            do {
    //          gmp_printf ("Testing %Zd*z+%Zd\n", a[1], a[0]);
                if ( zp2_rth_root (b, a, r, g, p) ) {
    //              gmp_printf ("rth root is %Zd*z+%Zd\n", b[1], b[0]);
                    zp2_exp_ui (c, b, r, p);
                    if ( ! zp2_equal (a, c) ) { gmp_printf ("Error, %Zd*z+%Zd is not a %d-root of %Zd*z+%Zd with z^2=%Zd mod %Zd\n",  b[1], b[0], r, a[1], a[0], p->s, p->p); abort(); }
                    cnt++;
                } else {
                }
                zp_inc (a[0], p);
            } while ( ! zp_is_zero(a[0]) );
            zp_inc (a[1], p);
        } while ( ! zp_is_zero (a[1]) );
        mpz_mul(w,p->p,p->p);  mpz_sub_ui(w,w,1);
        if ( mpz_divisible_ui_p (w, r) ) mpz_divexact_ui (w,w,r);
        mpz_add_ui (w,w,1);
        if ( mpz_cmp_ui(w,cnt) ) { gmp_printf ("Found %ld of %Zd expected residues\n", cnt, w); abort(); }
        gmp_printf ("Successfully computed %d-roots over all of Fp^2 for p=%Zd\n", r, p->p);
        zp2_clear (a); zp2_clear (b); zp2_clear (c); zp2_clear (g);
        zp_mod_clear (p);
        mpz_nextprime (pp, pp);
    }
    mpz_clear (pp);  mpz_clear (w);
}

#endif
