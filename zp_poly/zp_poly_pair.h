/*

zp_poly_pair.h: basic function for 2-vectors (pairs) and 2x2 matrices (quads)
			used to support gcd code in zp_poly

This module contains only the minimal functionality needed to support the gcd
code, it is not intended for general use.

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


#ifndef ZP_POLY_PAIR_H
#define ZP_POLY_PAIR_H

#include <gmp.h>
#include "zp_poly_util.h"
#include "zp_poly_main.h"

// a zp_poly_pair_t is a vector of length 2 over Z/pZ[X]
typedef struct zp_poly_pair_struct {
	mpz_t *f[2];
	int d[2];
} zp_poly_pair_t[1];

// a zp_poly_quad_t is a pair of column vectors of 2 polys representing a 2x2 matrix over Z/pZ[X]
typedef zp_poly_pair_t zp_poly_quad_t[2];

static inline int zp_poly_pair_mind(zp_poly_pair_t v) { return ( v->d[0] <= v->d[1] ? v->d[0] : v->d[1] ); }
static inline int zp_poly_pair_maxd(zp_poly_pair_t v) { return ( v->d[0] >= v->d[1] ? v->d[0] : v->d[1] ); }
static inline int zp_poly_pair_gap (zp_poly_pair_t v) { int gap = v->d[0]-v->d[1]; return ( gap < 0 ? -gap : gap ); }
static inline int zp_poly_pair_mini (zp_poly_pair_t v) { return ( v->d[0] <= v->d[1] ? 0 : 1 ); }
static inline int zp_poly_pair_maxi (zp_poly_pair_t v) { return ( v->d[0] >= v->d[1] ? 0 : 1 ); }
static inline int zp_poly_quad_maxd (zp_poly_quad_t A)
	{ int c1, c2;  c1 = zp_poly_pair_maxd(A[0]);  c2 = zp_poly_pair_maxd(A[1]);  return ( c1 > c2 ? c1 : c2 ); }

static inline void zp_poly_pair_copy (zp_poly_pair_t out, zp_poly_pair_t in)
{
	zp_poly_copy (out->f[0], in->f[0], in->d[0]);  out->d[0] = in->d[0];
	zp_poly_copy (out->f[1], in->f[1], in->d[1]);  out->d[1] = in->d[1];
}

static inline void zp_poly_pair_print (zp_poly_pair_t v)
	{ zp_poly_print(v->f[0], v->d[0]); printf ("    "); zp_poly_print(v->f[1], v->d[1]); puts (""); }

static inline void zp_poly_quad_print (zp_poly_quad_t A)
{
	printf ("    ");  zp_poly_print(A[0]->f[0], A[0]->d[0]); printf ("    "); zp_poly_print(A[1]->f[0], A[1]->d[0]); puts ("");
	printf ("    ");  zp_poly_print(A[0]->f[1], A[0]->d[1]); printf ("    "); zp_poly_print(A[1]->f[1], A[1]->d[1]); puts ("");
}

static inline void zp_poly_pair_init_on_stack (zp_poly_pair_t v, int d, zp_poly_stack_t stack)
	{ v->f[0] = zp_poly_stack_alloc(d, stack); v->f[1] = zp_poly_stack_alloc(d, stack); }

static inline void zp_poly_pair_copy_on_stack (zp_poly_pair_t out, zp_poly_pair_t in, zp_poly_stack_t stack)
	{ zp_poly_pair_init_on_stack (out, zp_poly_pair_maxd(in), stack); zp_poly_pair_copy(out, in); }
	
static inline void zp_poly_quad_init_on_stack (zp_poly_quad_t A, int d, zp_poly_stack_t stack)
	{ zp_poly_pair_init_on_stack (A[0], d, stack);  zp_poly_pair_init_on_stack (A[1], d, stack); }
	
static inline int zp_poly_pair_equal (zp_poly_pair_t a, zp_poly_pair_t b)
	{ return ( zp_poly_equal (a->f[0],a->d[0],b->f[0],b->d[0]) && zp_poly_equal (a->f[1],a->d[1],b->f[1],b->d[1]) ); }

// out may alias in
static inline void zp_poly_pair_shift_right (zp_poly_pair_t out, zp_poly_pair_t in, int d)
{
	register int i;
	
	out->d[0] = in->d[0] - d;
	if ( out->d[0] < 0 ) out->d[0] = -1;
	for ( i = 0 ; i <= out->d[0] ; i++ ) mpz_set(out->f[0][i],in->f[0][i+d]);
	out->d[1] = in->d[1] - d;
	if ( out->d[1] < 0 ) out->d[1] = -1;
	for ( i = 0 ; i <= out->d[1] ; i++ ) mpz_set(out->f[1][i],in->f[1][i+d]);
}

// A = identity matrix
static inline void zp_poly_quad_id (zp_poly_quad_t A)
	{ mpz_set_ui(A[0]->f[0][0],1); A[0]->d[0] = 0;  A[1]->d[0] = -1; A[0]->d[1] = -1; mpz_set_ui(A[1]->f[1][0],1); A[1]->d[1] = 0; }
		
// returns 1 of A == id, 0 if not
static inline int zp_poly_quad_is_id (zp_poly_quad_t A)
	{ if ( A[0]->d[0] || A[1]->d[1] || A[0]->d[1] >= 0 || A[1]->d[0] >= 0 ) return 0;  return ( mpz_is_one(A[0]->f[0][0]) && mpz_is_one(A[1]->f[1][0]) ? 1 : 0 ); }

// y = Ax  (y and x may coincide, and y may be a column of A)
static inline void zp_poly_quad_pair_mul_mod_xn (zp_poly_pair_t y, zp_poly_quad_t A, zp_poly_pair_t x, int n, mpz_t p, zp_poly_stack_t stack)
{
	zp_poly_mul42n (y->f[0], y->f[1], A[0]->f[0], A[0]->d[0], A[1]->f[0], A[1]->d[0], A[0]->f[1], A[0]->d[1], A[1]->f[1], A[1]->d[1], x->f[0], x->d[0], x->f[1], x->d[1], n, p);
	y->d[0] = zp_poly_degree (y->f[0], n-1);  y->d[1] = zp_poly_degree (y->f[1], n-1);
	return;
}

// C = AB (C may alias A or B, A and B may coincide)
static inline void zp_poly_quad_mul (zp_poly_quad_t C, zp_poly_quad_t A, zp_poly_quad_t B, mpz_t p, zp_poly_stack_t stack)
{
	int d;
	
       d = zp_imax(zp_imax(zp_imax(A[0]->d[0]+B[0]->d[0],A[1]->d[0]+B[0]->d[1]),zp_imax(A[0]->d[0]+B[1]->d[0],A[1]->d[0]+B[1]->d[1])),
			  zp_imax(zp_imax(A[0]->d[1]+B[0]->d[0],A[1]->d[1]+B[0]->d[1]),zp_imax(A[0]->d[1]+B[1]->d[0],A[1]->d[1]+B[1]->d[1])));
	zp_poly_mul4_mod_xn (C[0]->f[0],C[1]->f[0],C[0]->f[1],C[1]->f[1],
	                                    A[0]->f[0],A[0]->d[0],A[1]->f[0],A[1]->d[0],A[0]->f[1],A[0]->d[1],A[1]->f[1],A[1]->d[1],
	                                    B[0]->f[0],B[0]->d[0],B[1]->f[0],B[1]->d[0],B[0]->f[1],B[0]->d[1],B[1]->f[1],B[1]->d[1],d+1,p);
	C[0]->d[0] = zp_poly_degree (C[0]->f[0],d);  C[1]->d[0] = zp_poly_degree (C[1]->f[0],d); 
	C[0]->d[1] = zp_poly_degree (C[0]->f[1],d);  C[1]->d[1] = zp_poly_degree (C[1]->f[1],d);  
	return;
}

// A = E*A where E is the elementary matrix with -q in position (i,1-i) and 1's on the diagonal
static inline void zp_poly_quad_mul_lefte (zp_poly_quad_t A, mpz_t q[], int dq, int i, mpz_t p, zp_poly_stack_t stack)
{
	// q is expected to have low degree (typically 1), so don't try to optimize this
	zp_poly_submul (A[0]->f[i], A[0]->d+i, q, dq, A[0]->f[1-i], A[0]->d[1-i], p);
	zp_poly_submul (A[1]->f[i], A[1]->d+i, q, dq, A[1]->f[1-i], A[1]->d[1-i], p);
}

// A = A*E where E is the elementary matrix with -q in position (i,1-i) and 1's on the diagonal
static inline void zp_poly_quad_mul_righte (zp_poly_quad_t A, mpz_t q[], int dq, int i, mpz_t p, zp_poly_stack_t stack)
{
	// q is expected to have low degree (typically 1), so don't try to optimize this
	zp_poly_addmul (A[1-i]->f[0], A[1-i]->d+0, q, dq, A[i]->f[0], A[i]->d[0], p);
	zp_poly_addmul (A[1-i]->f[1], A[1-i]->d+1, q, dq, A[i]->f[1], A[i]->d[1], p);
}

#endif
