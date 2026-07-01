/*

zp2.h: support for some basic operations in (Z/pZ)[X]/(X^2-s) ~ Fp^2
	   primarliy used to solve cubics in Z/pZ[X] via Cardano's method

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

#include <stdlib.h>
#include "gmp.h"
#include "zp_poly_util.h"
#include "zp_roots.h"
#include "zp2.h"

// note that we are NOT ALLOWED to use p->w or p->t, these are reserved for zp2.h

void _zp2_exp_sw (zp2_t b, zp2_t a, mpz_t e, int k, zp_mod_t p);

// compute b = a^e in Fp[x]/(x^2-s) for e >= 0 using standard LR binary exponentiation, a and b may coincide
// optimized to take advantage of a^(p+1)=N(a)
void zp2_exp (zp2_t b, zp2_t a, mpz_t e, zp_mod_t p)
{
	zp2_t c, w;
	int i;

	assert ( mpz_sgn(e) >= 0 );
	// handle trivial cases first: the exponents 0, 1, 2, 3, 4 arise often enough to check for immediately
	if ( mpz_cmp_ui(e,4) <= 0 ) {
		if ( mpz_is_zero(e) ) { zp2_set_one (b); return; }
		if ( mpz_is_one(e) ) { zp2_set (b, a); return; }
		if ( mpz_cmp_ui(e,2) == 0 ) { zp2_sqr (b, a, p); return; }
		if ( mpz_cmp_ui(e,4) == 0 ) { zp2_sqr (b, a, p); zp2_sqr (b, b, p); return; }
		zp2_init (c, p);  zp2_sqr (c, a, p); zp2_mul (b, a, c, p); zp2_clear (c); return;
	}
	
	if ( zp2_is_in_zp (a) ) {
		zp_exp (b[0], a[0], e, p);
		zp_set_zero (b[1]); 
		return;
	} else if ( zp_is_zero (a[0]) ) {
		mpz_t e1;
		
		// (a[1]*z)^e = (a[1]^2*s)^floor(e/2) * (a[1]*z)^(e mod 2)
		mpz_init_set (e1, e);  mpz_div_2exp (e1, e1, 1); 
		zp_sqr (b[0], a[1], p);  zp_mul_ui (b[0], b[0], p->s, p);
		zp_exp (b[0], b[0], e1, p);
		if ( mpz_tstbit (e,0) ) { zp_mul(b[1], b[0], a[1], p); zp_set_zero (b[0]); } else zp_set_zero (b[1]);		// note a[1] is preserved above even if a and b coincide
		mpz_clear (e1);
		return;
	}
	if ( mpz_cmp(e, p->p) > 0 ) {
		mpz_t e1, e2, n;

		mpz_init_set (e1, e);  mpz_init2 (e2, 2*zp_bits(p));  mpz_init_modp (n, p->p);
		// reduce e mod p^2-1 first then compute two digits mod p+1
		mpz_mul (e2, p->p, p->p);  mpz_sub_ui (e2, e2, 1);  mpz_mod (e2, e1, e2);
		mpz_add_ui (n, p->p, 1); mpz_fdiv_qr (e2, e1, e2, n);
		// use the fact that a^(p+1) = N(a)
		zp2_norm (n, a, p);  zp_exp (n, n, e2, p);
		zp2_exp (b, a, e1, p);		// recursive call will return because e1 <= p
		zp2_mul_zp (b, b, n, p);
		mpz_clear (e1); mpz_clear (e2); mpz_clear(n);
		return;
	}
	_zp2_exp_sw (b, a, e, 5, p);
	return;
	
	zp2_init (w, p);  zp2_init (c, p);	
	zp2_set (c, a);		// copy in case a == b
	zp2_set (b, a);
	for ( i = mpz_sizeinbase (e,2) - 2 ; i >= 0 ; i-- ) {
		zp2_sqr (b, b, p);
		if ( mpz_tstbit (e, i) ) zp2_mul (b, b, c, p);
	}
	zp2_clear (w);  zp2_clear (c);
}

// standard sliding window exponentiation (Alg 9.10 in the HEHECC, p. 150)
void _zp2_exp_sw (zp2_t b, zp2_t a, mpz_t e, int k, zp_mod_t p)
{
	zp2_t *A;
	int i, j, m, n, s;
		
	// precompute A[i] = a^(2i-1) for 0 <= i <= 2^(k-1)
	n = 1<<(k-1);
	A = zp_malloc (n*sizeof(*A));
	for ( i = 0 ; i < n ; i++ ) zp2_init (A[i], p);
	zp2_set (A[0], a);  zp2_sqr (b, a, p);
	for ( i = 1 ; i < n ; i++ ) zp2_mul (A[i],A[i-1], b, p);
	
	zp2_set_one (b);
	i = mpz_sizeinbase (e, 2) -1;
	while ( i >= 0 ) {
		if ( ! mpz_tstbit(e,i) ) { zp2_sqr (b, b, p); i--; continue; }		// slide
		s = i-k+1;  if ( s < 0 ) s = 0;
		while ( ! mpz_tstbit (e,s) ) s++;
		for ( j = 0 ; j < i-s+1 ; j++ ) zp2_sqr (b, b, p);
		for ( m = 1, j = i-1 ; j >= s ; j-- ) { m <<= 1;  m += mpz_tstbit(e,j); }
		zp2_mul (b, b, A[(m-1)>>1], p);
		i = s-1;
	}
	
	for ( i = 0 ; i < n ; i++ ) zp2_clear (A[i]);
	zp_free (A, n*sizeof(*A));
	return;
}


// dual binary RL exponentiation: b1 = a^e2, b2=a^e2 for nonneg e1 and e2
void zp2_exp2 (zp2_t b1, zp2_t b2, zp2_t a, mpz_t e1, mpz_t e2, zp_mod_t p)
{
	zp2_t x, b3, w;
	int i, j, n;

	if ( zp2_is_in_zp (a) ) { zp_set_zero (b1[1]);  zp_set_zero (b2[1]);  zp_exp2 (b1[0], b2[0], a[0], e1, e2, p);  return; }
	if ( mpz_cmp (e1,p->p) > 0 ) {
		mpz_t e11, e12, n1, n2;
		mpz_init_set (e11, e1); mpz_init2 (e12, 2*zp_bits(p)); mpz_init_modp (n1, p->p); mpz_init_modp (n2, p->p);
		mpz_mul (e12, p->p, p->p);  mpz_sub_ui (e12, e12, 1); mpz_mod (e12, e11, e12);
		mpz_add_ui (n1, p->p, 1);  mpz_fdiv_qr (e12, e11, e12, n1);
		if ( mpz_cmp (e2, p->p) > 0 ) {
			mpz_t e21, e22;
			mpz_init_set (e21, e2); mpz_init2 (e22, 2*zp_bits(p));
			mpz_mul (e22, p->p, p->p);  mpz_sub_ui (e22, e22, 1);  mpz_mod (e22, e21, e22);
			mpz_fdiv_qr (e22, e21, e22, n1);
			zp2_norm (n1, a, p);
			zp_exp2 (n1, n2, n1, e12, e22, p);
			zp2_exp2 (b1, b2, a, e11, e21, p);
			zp2_mul_zp (b1, b1, n1, p);  zp2_mul_zp (b2, b2, n2, p);
			mpz_clear (e21); mpz_clear(e22);
		} else {
			zp2_norm (n1, a, p);
			zp_exp (n1, n1, e12, p);
			zp2_exp2 (b1, b2, a, e11, e2, p);
			zp2_mul_zp (b1, b1, n1, p);
		}
		mpz_clear (e11); mpz_clear(e12); mpz_clear (n1); mpz_clear (n2);
		return;
	}
	if ( mpz_cmp (e2, p->p) > 0 ) { zp2_exp2 (b2, b1, a, e2, e1, p); return; }

	zp2_init (x, p);  zp2_init (b3, p); zp2_init (w, p);
	i = mpz_sizeinbase (e1, 2);
	n = mpz_sizeinbase (e2,2);
	if ( i > n ) n = i;
	zp2_set (x, a);
	zp2_set_one (b1);  zp2_set_one (b2); zp2_set_one (b3);
	for ( i = 0 ; i < n ; i++ ) {
		j = mpz_tstbit(e1,i) + (mpz_tstbit(e2,i)<<1);
		switch (j) {
		case 1: zp2_mul (b1, b1, x, p); break;
		case 2: zp2_mul (b2, b2, x, p); break;
		case 3: zp2_mul (b3, b3, x, p); break;
		}
		if ( i < n-1 ) zp2_mul (x, x, x, p);
	}
	zp2_mul (b1, b1, b3, p);
	zp2_mul (b2, b2, b3, p);
	zp2_clear (b3); zp2_clear (x); zp2_clear (w);
}
