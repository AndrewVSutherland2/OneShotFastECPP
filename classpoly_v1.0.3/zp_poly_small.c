/*

zp_poly_small.c: fast basic arithmetic functions for low degree cases

Copyright (C) 2010-2012, David Harvey and Andrew V. Sutherland

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

#include "zp_poly_util.h"
#include "zp_poly_mul.h"
#include "zp_poly_main.h"

// h = f1*f2.  Assumes either d1 or d2 (or both) is very small.
void zp_poly_mul_mod_xn_small (mpz_t h[], mpz_t f1[], int d1, mpz_t f2[], int d2, int n, zp_p _p)
{
	zp_p_to_mpz p=_p;
	mpz_t *f,*F, w;
	register int i, j, k, d, D;
	
	// make F of deg D the large poly, and f of deg d the smaller one
	if ( d1 < d2 ) { f = f1;  d = d1;  F = f2; D = d2; } else { f = f2; d = d2; F = f1; D = d1; }
	if ( d < 0 ) { zp_poly_zero_fill(h,n-1); return; }			// fast mult by 0, but we need to zero fill
	if ( zp_poly_is_one(f,d) )								// fast mult by 1
		{ zp_poly_copy (h, F, zp_imin(D,n-1)); zp_poly_zero_pad(h,D+1,n-1); return; }

	// work top down (necessary in case h =f1 or f2)
	mpz_init_modp_w (w, p);
	for ( k = n-1 ; k > D+d ; k-- ) mpz_set_ui(h[k], 0);			// pad to length n with leading zeros if needed
	for ( ; k > D ; k-- ) {						 		// high triangle
		i = k-d;  j = d;
		mpz_mul(w,F[i++],f[j--]);
		while ( i <= D ) mpz_addmul (w,F[i++],f[j--]);
		mpz_mod(h[k], w, p);
	}
	for ( ; k >= d ; k-- ) {								// middle rectangle
		i = k-d; j = d;
		mpz_mul(w,F[i++],f[j--]);
		while ( j >= 0 ) mpz_addmul (w,F[i++],f[j--]);
		mpz_mod(h[k], w, p);
	}
	for ( ; k >= 0 ; k-- ) {								// low triangle
		i = 0; j = k;
		mpz_mul(w,F[i++],f[j--]);
		while ( j >= 0 ) mpz_addmul (w,F[i++],f[j--]);
		mpz_mod(h[k], w, p);
	}
	mpz_clear (w);
}

// h += f1*f2
void zp_poly_addmul_small (mpz_t h[], int *dh, mpz_t f1[], int d1, mpz_t f2[], int d2, zp_p _p)
{
	zp_p_to_mpz p=_p;
	mpz_t *f,*F, w;
	register int i, j, k, d, D;
	
	// make F of deg D the large poly, and f of deg d the smaller one
	if ( d1 < d2 ) { f = f1;  d = d1;  F = f2; D = d2; } else { f = f2; d = d2; F = f1; D = d1; }
	if ( d < 0 ) return;	// no need to zero-fill, degree is known
	if ( zp_poly_is_one(f,d) ) { *dh = zp_poly_add (h, h, *dh, F, D, p); return; }
	
	// work top down (necessary in case h =f1 or f2)
	mpz_init_modp (w, p);
	for ( k = D+d ; k > D ; k-- ) {		// high triangle
		i = k-d;  j = d;
		mpz_mul(w,F[i++],f[j--]);
		while ( i <= D ) mpz_addmul (w,F[i++],f[j--]);
		if ( k <= *dh ) mpz_add(w,w,h[k]);
		mpz_mod(h[k], w, p); 
	}
	for ( ; k >= d ; k-- ) {			// middle rectangle
		i = k-d; j = d;
		mpz_mul(w,F[i++],f[j--]);
		while ( j >= 0 ) mpz_addmul (w,F[i++],f[j--]);
		if ( k <= *dh ) mpz_add(w,w,h[k]);
		mpz_mod(h[k], w, p);
	}
	for ( ; k >= 0 ; k-- ) {			// low triangle
		i = 0; j = k;
		mpz_mul(w,F[i++],f[j--]);
		while ( j >= 0 ) mpz_addmul (w,F[i++],f[j--]);
		if ( k <= *dh ) mpz_add(w,w,h[k]);
		mpz_mod(h[k], w, p);
	}
	mpz_clear (w);
	if ( d+D >= *dh ) *dh = zp_poly_degree (h,d+D);
}

// h -= f1*f2
void zp_poly_submul_small (mpz_t h[], int *dh, mpz_t f1[], int d1, mpz_t f2[], int d2, zp_p _p)
{
	zp_p_to_mpz p=_p;
	mpz_t *f,*F, w;
	register int i, j, k, d, D;

	// make F of deg D the large poly, and f of deg d the smaller one
	if ( d1 < d2 ) { f = f1;  d = d1;  F = f2; D = d2; } else { f = f2; d = d2; F = f1; D = d1; }
	if ( d < 0 ) return;	// no need to zero-fill, degree is known
	if ( zp_poly_is_one(f,d) ) { *dh = zp_poly_sub (h, h, *dh, F, D, p); return; }

	// work top down (necessary in case h =f1 or f2)
	mpz_init_modp_w (w, p);
	for ( k = D+d ; k > D ; k-- ) {		// high triangle
		i = k-d;  j = d;
		mpz_mul(w,F[i++],f[j--]);
		while ( i <= D ) mpz_addmul (w,F[i++],f[j--]);
		if ( k <= *dh ) mpz_sub(w,w,h[k]);
		mpz_neg (w,w); mpz_mod(h[k], w, p); 
	}
	for ( ; k >= d ; k-- ) {			// middle rectangle
		i = k-d; j = d;
		mpz_mul(w,F[i++],f[j--]);
		while ( j >= 0 ) mpz_addmul (w,F[i++],f[j--]);
		if ( k <= *dh ) mpz_sub(w,w,h[k]);
		mpz_neg (w,w); mpz_mod(h[k], w, p);
	}
	for ( ; k >= 0 ; k-- ) {			// low triangle
		i = 0; j = k;
		mpz_mul(w,F[i++],f[j--]);
		while ( j >= 0 ) mpz_addmul (w,F[i++],f[j--]);
		if ( k <= *dh ) mpz_sub(w,w,h[k]);
		mpz_neg (w,w); mpz_mod(h[k], w, p);
	}
	mpz_clear (w);
	if ( d+D >= *dh ) *dh = zp_poly_degree (h,d+D);
}

// q = a div b for monic b.  Assumes da-db is very small.  q may equal a but not b
void zp_poly_div_small_monic (mpz_t q[], mpz_t a[], int da, mpz_t b[], int db, zp_p _p)
{
	zp_p_to_mpz p=_p;
	mpz_t w;
	register int i, j, k, dq;
	
	assert ( zp_poly_is_monic(b,db) && q != b );
	dq = da - db;
	if ( dq < 0 ) return;			// quotient is zero, nothing to do
	// compute quotient from high to low
	mpz_init_modp (w, p);
	for ( k = 0 ; k <= dq ; k++ ) {
		mpz_set(w,a[da-k]);
		i = dq; j = db-k;  if ( j < 0 ) { i += j; j = 0; }
		while ( j < db ) mpz_submul(w,q[i--],b[j++]);
		mpz_mod (q[dq-k], w, p);
	}
	mpz_clear (w);
}

// q = a div b with b not nesc. monic (faster if it is).  Assumes da-db is very small.  q may alias a but not b
void zp_poly_div_small  (mpz_t q[], mpz_t a[], int da, mpz_t b[], int db, zp_p _p)
{
	zp_p_to_mpz p=_p;
	mpz_t c, w;
	register int i, j, k, dq;
	
	if ( zp_poly_is_monic(b,db) ) { zp_poly_div_small_monic (q, a, da, b, db, p);  return; }
	dq = da - db;
	if ( dq < 0 ) return;			// quotient is zero, nothing to do
	assert ( q != b );
	mpz_init_modp (c, p); mpz_init_modp_w (w, p);
	if ( ! mpz_invert (c,b[db],p) ) { fprintf (stderr, "division by zero in zp_poly_div_small\n"); abort(); }
	// compute quotient from high to low
	for ( k = 0 ; k <= dq ; k++ ) {
		mpz_set(w,a[da-k]);
		i = dq; j = db-k;  if ( j < 0 ) { i += j; j = 0; }
		while ( j < db ) mpz_submul(w,q[i--],b[j++]);
		mpz_mod (w, w, p);
		mpz_mul (w, w, c);
		mpz_mod (q[dq-k], w, p);
	}
	mpz_clear (c); mpz_clear (w);
}

// r = a - qb where q is the Euclidean quotient of a/b and da >= db > dr
// r is filled to deg db-1, the degree of r is returned
int zp_poly_rem_small (mpz_t r[], mpz_t a[], int da, mpz_t b[], int db, mpz_t q[], zp_p _p)
{
	zp_p_to_mpz p=_p;
	mpz_t w;
	register int i, k, d, dq;
	
	assert ( da >= db && db >= 0 );
	if ( da < 0 ) { for ( i = 0 ; i < db ; i++ ) mpz_set_zero (r[i]); return -1; }
	mpz_init_modp (w, p);
	dq = da-db;
	for ( k = db-1 ; k >= 0 ; k-- ) {
		d = ( k < dq ? k : dq );
		mpz_set(w,a[k]);
		for ( i = 0 ; i <= d ; i++ ) mpz_submul (w, q[i], b[k-i]);
		mpz_mod (r[k],w,p);
	}
	mpz_clear (w);
	return zp_poly_degree (r, db-1);
}

// g = f * (x - r) where f is a poly in Z/pZ[x].  f and g may point to the same place.
void zp_poly_add_root (mpz_t g[], mpz_t f[], int d, mpz_t r, zp_p _p)
{
	zp_p_to_mpz p=_p;
	mpz_t w;
	register int i;
	
	if ( d < 0 ) { mpz_set_one(g[1]); mpz_neg_mod(g[0],r,p); return; }
	mpz_init_modp_w (w, p);
	mpz_set(g[d+1],f[d]);
	for ( i = d ; i > 0 ; i-- ) {
		mpz_mul_mod (w, r, f[i], p);
		mpz_sub_mod (g[i],f[i-1],w,p);
	}
	mpz_mul_mod (w, r, f[0], p);
	mpz_neg_mod (g[0], w, p);
	mpz_clear (w);
}

// computes g=f/(x-r), assuming f(r) = 0.  f and g may point to the same place
void zp_poly_remove_root (mpz_t g[], mpz_t f[], int d, mpz_t r, zp_p _p)
{
	zp_p_to_mpz p=_p;
	mpz_t t1, t2;
	register int i;
	
	assert (d > 0);
	mpz_init_modp (t1, p); mpz_init_modp (t2, p);
	mpz_set (t1, f[d]);
	for ( i = d-1 ; i > 0 ; i-- ) {
		mpz_set (t2,f[i]);
		mpz_set (g[i],t1);
		mpz_mul_mod (t1, t1, r, p);
		mpz_add_mod (t1, t1, t2, p);
	}
	mpz_set(g[0],t1);
	mpz_clear (t1); mpz_clear (t2);
}

// replace monic f(x) by f(x-f[d-1]/d), killing x^(d-1) coefficient
// sets s to f[d-1]/d, if r is a root of the new f then r-s is a root of the old f
void zp_poly_depress_monic_inplace (mpz_t f[], int d, mpz_t s, zp_p _p)
{
	zp_p_to_mpz p=_p;
	mpz_t t;
	register int i, j;

	if ( mpz_is_zero(f[d-1]) ) { mpz_set_zero(s); return; }
	mpz_init_modp(t,p);
	mpz_set_ui(s,d);
	if (! mpz_invert(t,s,p) ) { gmp_printf ("Cannot depress degree %d poly mod %Zd, division by zero\n", d, p); abort(); }
	mpz_mul_mod(s,t,f[d-1],p);				// s = f[d-1]/d
	mpz_neg_mod(s,s,p);					// We will replace f by f(x-f[d-1]/d)
	for ( i = d ; i ; i-- ) {
		for ( j = i ; j < d ; j++ ) { mpz_mul_mod(t,s,f[j],p);  mpz_add_mod(f[j-1],f[j-1],t,p); }
		mpz_add_mod(f[j-1],f[j-1],s,p);
	}
	mpz_neg_mod(s,s,p);
	mpz_clear (t);
}
