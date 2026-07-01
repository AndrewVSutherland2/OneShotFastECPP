/*

zp_roots.c: functions for computing rth roots and solving degree 2 and 3 equations in Fp=Z/pZ and Fp^2

Copyright (C) 2011, Andrew V. Sutherland

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

#ifndef _ZP_ROOTS_INCLUDE_
#define _ZP_ROOTS_INCLUDE_

#include "zp_poly_util.h"
#include "zp2.h"

// b = sqrt(a), returns 0 if no sqrt exist and -1 if p is found to be composite
int zp_sqrt (zp_t b, zp_t a, zp_mod_t p);

// functions for computing 4th and 8th roots  are optimized for p=3 mod 4 and will revert to iterated square roots otherwise
int zp_fourth_root (zp_t b, zp_t a, zp_mod_t p);
int zp_sixth_root (zp_t b, zp_t a, zp_mod_t p); 
int zp_eighth_root (zp_t b, zp_t a, zp_mod_t p);

// computes a generator for the r-sylow of (Z/pZ)* and an rth-root of unity (if one exists).  returns r-adic valuation of p-1
int zp_sylow_gen (zp_t g, zp_t z, int r, zp_mod_t p);

// b = a^(1/r), returns 0 if no rth root exists.  r must be an odd prime.  g is a generator of the r-sylow
int zp_rth_root (zp_t b, zp_t a, int r, zp_t g, zp_mod_t p);

// b = a^(1/3), returns 0 if no cube-root exists
static inline int zp_cbrt (zp_t b, zp_t a, zp_mod_t p)
{
	if ( p->p1mod3 && zp2_is_zero (p->g3) )
		zp_sylow_gen (p->g3[0], p->z3[0], 3, p);		// if p=1 mod 3 then 3-sylow of Fp^2 lies in Fp and ow we don't need it
	return zp_rth_root (b, a, 3, p->g3[0], p);
}

// all of the functions below assume the given poly is monic and do not look at the leading coefficient (so it need not be present and r and f may coincide)

// finds roots of a monic quadratic, returns number of roots found (0 or 2).  r and f may coincide
int zp_poly_quadratic_roots (zp_t r[2], zp_t f[2], zp_mod_t p);
int zp2_poly_quadratic_roots (zp2_t r[2], zp2_t f[2], zp_mod_t p);

// finds roots of a monic cubic, returns number of roots found (0, 1, or 3).  r and f may coincide
int zp_poly_cubic_roots (zp_t r[3], zp_t f[3], zp_mod_t p);
int zp2_poly_cubic_roots (zp2_t r[3], zp2_t f[3], zp_mod_t p);

// computes the parity of the # of irreducible factors of a monic cubic f via Stickelberger
int zp_poly_cubic_parity (zp_t f[3], zp_mod_t p);
int zp2_poly_cubic_parity (zp2_t f[3], zp_mod_t p);

// see zp_poly.h for zp_poly_quartic_roots definition which is currently not implemented via radicals

#endif
