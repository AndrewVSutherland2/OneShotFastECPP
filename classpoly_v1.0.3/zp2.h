/*

zp2.h: support for some basic operations in (Z/pZ)[X]/(X^2-s) ~ Fp^2
	   primarliy used to solve cubics in Z/pZ[X] via Cardona's method

Copyright (C) 2011 Andrew V. Sutherland

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

#ifndef _ZP2_INCLUDE_
#define _ZP2_INCLUDE_

#include "zp_poly_util.h"


typedef zp_t zp2_t[2];

static inline void zp2_init (zp2_t t, zp_mod_t p) { zp_init (t[0], p); zp_init (t[1], p); }
static inline void zp2_clear (zp2_t t) { zp_clear (t[0]); zp_clear (t[1]); }
static inline void zp2_random (zp2_t r, zp_mod_t p) { zp_random (r[0], p);  zp_random (r[1], p); }

// all of the functions below work in the field (Z/pZ)[X]/(X^2-s), where s is a given non-residue in (Z/pZ)*

static inline int zp2_is_in_zp (zp2_t t) { return zp_is_zero(t[1]); }
static inline void zp2_bar (zp2_t b, zp2_t a, zp_mod_t p) { zp_neg (b[1], a[1], p); zp_set (b[0], a[0]); }

// N = norm(a)
static inline void zp2_norm (zp_t N, zp2_t a, zp_mod_t p)
{
	zp_sqr (p->t[0], a[1], p);
	mpz_mul_ui (p->t[0], p->t[0], p->s);
	mpz_mul (p->t[1], a[0], a[0]);
	mpz_sub (p->t[1], p->t[1], p->t[0]);
	mpz_mod (N, p->t[1], p->p);		// N(a) = a[0]^2 - s*a[1]^2
}

// b= a^2 (overlap ok)
static inline void zp2_sqr (zp2_t b, zp2_t a, zp_mod_t p)
{
// use Karatsuba to get 3S rather than 2S + 1M -- in testing this is actually slower so don't use
//	mpz_mul (p->t[1], a[1], a[1]); mpz_mul (p->t[0], a[0], a[0]);
//	mpz_add (b[1], a[0], a[1]); mpz_mul (b[1], b[1], b[1]); mpz_sub (b[1], b[1], p->t[1]);  mpz_sub (b[1], b[1], p->t[0]); mpz_mod (b[1], b[1], p->p);
//	mpz_mul_ui (b[0], p->t[1], p->s); mpz_add (b[0], b[0], p->t[0]); mpz_mod (b[0], b[0], p->p);
	// be careful to just use 2 mod reductions
	mpz_mul (p->t[1], a[1], a[1]); 
	mpz_mul (p->t[0], a[0], a[1]);  mpz_mul_2exp (p->t[0], p->t[0], 1);  mpz_mod (b[1], p->t[0], p->p);									// b1 = 2*a0*a1
	mpz_mul (p->t[0], a[0], a[0]);  mpz_mul_ui (p->t[1], p->t[1], p->s); mpz_add (p->t[0], p->t[0], p->t[1]);  mpz_mod(b[0], p->t[0], p->p);	// b0 = a0^2 + s*a1^2
}

// c= a * b (overlap ok)
static inline void zp2_mul (zp2_t c, zp2_t b, zp2_t a, zp_mod_t p)
{
	mpz_add (p->w[0] ,a[0],a[1]); mpz_add (p->w[1], b[0], b[1]); mpz_mul (p->w[0], p->w[0], p->w[1]);				// w[0] = (a[0]+a[1])*(b[0]+b[1])
	mpz_mul (p->t[0], a[0],b[0]);  mpz_mul (p->t[1], a[1], b[1]);											// t[0] = a[0]*b[0], t[1]=a[1]*b[1]
	mpz_sub (p->w[1], p->w[0], p->t[0]); mpz_sub (p->w[1], p->w[1], p->t[1]);								// w[1] = w[0]-t[0]-t[1] = a[0]b[1]+a[1]b[0]
	mpz_mul_ui (p->w[0], p->t[1], p->s); mpz_add (p->w[0], p->w[0], p->t[0]);					 			// w[0] = a[0]b[0]+a[1]b[1]s
	mpz_mod (c[0], p->w[0], p->p);  mpz_mod (c[1], p->w[1], p->p);										// reduce into c
}

// c = a * b
static inline void zp2_mul_zp (zp2_t c, zp2_t a, mpz_t b, zp_mod_t p)
	{ mpz_mul (p->t[0],a[0],b);  mpz_mod(c[0],p->t[0],p->p);  mpz_mul (p->t[1],a[1],b);  mpz_mod(c[1],p->t[1],p->p); }

// c = a * b
static inline void zp2_mul_ui (zp2_t c, zp2_t a, unsigned long b, zp_mod_t p)
	{ mpz_mul_ui (p->t[0],a[0],b);  mpz_mod(c[0],p->t[0],p->p);  mpz_mul_ui (p->t[1],a[1],b);  mpz_mod(c[1],p->t[1],p->p); }

// c = a + b
static inline void zp2_add (zp2_t c, zp2_t a, zp2_t b, zp_mod_t p)
	{ zp_add (c[0],a[0],b[0],p);  zp_add (c[1],a[1],b[1],p); }

// c = a + b
static inline void zp2_add_zp (zp2_t c, zp2_t a, zp_t b, zp_mod_t p)
	{ zp_add (c[0],a[0],b,p); }

// c = a + b
static inline void zp2_add_ui (zp2_t c, zp2_t a, unsigned long b, zp_mod_t p)
	{ zp_set (c[1],a[1]); zp_add_ui (c[0],a[0],b,p); }

// c = a - b
static inline void zp2_sub (zp2_t c, zp2_t a, zp2_t b, zp_mod_t p)
	{ zp_sub (c[0],a[0],b[0],p);  zp_sub (c[1],a[1],b[1],p); }

// c = a - b
static inline void zp2_sub_ui (zp2_t c, zp2_t a, unsigned long b, zp_mod_t p)
	{ zp_set (c[1],a[1]); zp_sub_ui (c[0],a[0],b,p); }

// c = a - b
static inline void zp2_sub_zp (zp2_t c, zp2_t a, mpz_t b,zp_mod_t p)
	{ zp_sub (c[0],a[0],b,p); }

// b = a/2
static inline void zp2_div2 (zp2_t b, zp2_t a, zp_mod_t p)
	{ zp_div2 (b[0],a[0],p);  zp_div2 (b[1],a[1],p); }

// b = a/2
static inline void zp2_div3 (zp2_t b, zp2_t a, zp_mod_t p)
	{ zp_div3 (b[0],a[0],p);  zp_div3 (b[1],a[1],p); }

// b = -a
static inline void zp2_neg (zp2_t b, zp2_t a, zp_mod_t p)
	{ zp_neg (b[0],a[0],p);  zp_neg (b[1],a[1],p); }
	
// a++
static inline void zp2_inc (zp2_t a, zp_mod_t p) { zp_inc (a[0], p); }

// a--
static inline void zp2_dec (zp2_t a, zp_mod_t p) { zp_dec (a[0], p); }

// a==b
static inline int zp2_equal (zp2_t a, zp2_t b) { return ( zp_equal(a[0],b[0]) && zp_equal(a[1],b[1]) ); }
	
// a=b
static inline void zp2_set (zp2_t b, zp2_t a) { zp_set(b[0],a[0]);  zp_set(b[1],a[1]); }
static inline void zp2_set_zp (zp2_t b, zp_t a) { zp_set(b[0],a);  zp_set_zero(b[1]); }
static inline void zp2_set_ui (zp2_t b, unsigned long a, zp_mod_t p) { zp_set_ui(b[0],a, p);  zp_set_zero(b[1]); }

// a=0
static inline void zp2_set_zero (zp2_t a) { zp_set_zero(a[0]);  zp_set_zero(a[1]); }

// a==0
static inline int zp2_is_zero (zp2_t a) { return ( zp_is_zero(a[0]) && zp_is_zero(a[1]) ); }

// a=1
static inline void zp2_set_one (zp2_t a) { zp_set_one(a[0]);  zp_set_zero(a[1]); }

// a==1
static inline int zp2_is_one (zp2_t a) { return ( zp_is_one(a[0]) && zp_is_zero(a[1]) ); }

// b = 1/a = bar(a)/N(a)
static inline void zp2_invert (zp2_t b, zp2_t a, zp_mod_t p)
	{ zp2_norm (p->w[0], a, p);  zp_invert (p->w[0], p->w[0], p); zp2_bar(b, a, p); zp2_mul_zp (b, b, p->w[0], p); }		// note that norm and mul_zp use t but not w

// b = a^e for nonnegative e
void zp2_exp (zp2_t b, zp2_t a, mpz_t e, zp_mod_t p);
static inline void zp2_exp_ui (zp2_t b, zp2_t a, unsigned long e, zp_mod_t p)
    { mpz_t t;  mpz_init_set_ui (t, e); zp2_exp (b, a, t, p); mpz_clear (t); }											// we can't use temps in p across the call to zp2_exp!
void zp2_exp2 (zp2_t b1, zp2_t b2, zp2_t a, mpz_t e1, mpz_t e2, zp_mod_t p);
    
// b = sqrt(a) , returns 0 if no such b exists
int zp2_sqrt (zp2_t b, zp2_t a, zp_mod_t p);

// randomized algorithm to find a generator g for the r-Sylow subgroup of (Z/pZ)*
// also computes a primitive rth root of unity z, if one exists in (Z/pZ)* (ow z is set to zero)
// returns the r-adic valuation of order of the r-Sylow subgroup
int zp2_sylow_gen (zp2_t g, zp2_t z, int r, zp_mod_t p);

// b = a^(1/4), given r-sylow generator g.  returns 0 if no such b exits
int zp2_rth_root (zp2_t b, zp2_t a, int r, zp2_t g, zp_mod_t p);
	
// given monic quadratic x^2 + b*x + c computes roots r1 and r2 or returns 0 if none exist.
// r1=b and r2=c is ok
static inline int zp2_quadratic_roots (zp2_t r1, zp2_t r2, zp2_t b, zp2_t c, zp_mod_t p)
{
	zp2_sqr (p->w, b, p);									// note that z2_sqr only uses t, not w
	zp2_add (r2, c, c, p);  zp2_add (r2, r2, r2, p);
	zp2_sub (r2, p->w, r2, p);								// r2 = b^2-4c
	if ( ! zp2_sqrt (r2, r2, p) ) return 0;						// r2 = sqrt(b^2-4c) (can't use w across the call to zp2_sqrt)
	zp2_sub (r1, r2, b, p);
	zp2_div2 (r1, r1, p);										// r2 = (-b+sqrt(b^2-4c)) / 2
	zp2_sub (r2, r1, r2, p);									// r1 = (-b-sqrt(b^2-4c)) / 2
	return 2;
}

#endif
