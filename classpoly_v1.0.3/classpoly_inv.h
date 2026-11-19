#ifndef _CLASSPOLY_INV_INCLUDE_
#define _CLASSPOLY_INV_INCLUDE_

#include "ff_poly.h"
#include "class_inv.h"
#include "classpoly.h"

/*
    Copyright 2012 Andrew V. Sutherland

    This file is part of classpoly.

    classpoly is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    classpoly is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with classpoly.  If not, see <http://www.gnu.org/licenses/>.
*/

typedef struct classpoly_inv_struct {
	int inv;
	int p1, p2;
	int inverted, negated, units;
	int ambiguous_H, ambiguous_sign;
	mpz_t T;
	int use_trace;
	int signfree;					// set when the caller does not need the sign of an ambiguous_sign invariant fixed (the trace T is then not computed)
	long Nrep[IQ_MAX_GENS];
} classpoly_inv_t[1];

// true when the roots mod p of this invariant's class polynomial are determined only up to a common sign that takes the trace to fix
static inline int classpoly_inv_needs_trace (int inv, long enum_cnt)
	{ return inv_sqrt_invariant (inv) && ! ((enum_cnt&1) && inv_units (inv)); }

int classpoly_inv_setup_signfree (classpoly_inv_t Hinv, classpoly_t H, int inv, classpoly_crt_t crt, int signfree);
static inline int classpoly_inv_setup (classpoly_inv_t Hinv, classpoly_t H, int inv, classpoly_crt_t crt)
	{ return classpoly_inv_setup_signfree (Hinv, H, inv, crt, 0); }
void classpoly_inv_clear (classpoly_inv_t Hinv);

void _ff_classpoly_inv_from_j (ff_t *F, ff_t *J, classpoly_inv_t Hinv, classpoly_t H);
static inline void ff_classpoly_inv_from_j (ff_t *F, ff_t *J, classpoly_inv_t Hinv, classpoly_t H)
	{ if ( Hinv->p1 > 1 ) _ff_classpoly_inv_from_j (F,J,Hinv,H); else ff_inv_from_j (F,J,Hinv->inv); }
	
int _ff_classpoly_inv_disambiguate (classpoly_t H, classpoly_inv_t Hinv);
static inline int ff_classpoly_inv_disambiguate (classpoly_t H, classpoly_inv_t Hinv)
	{ if ( Hinv->ambiguous_H || Hinv->ambiguous_sign ) return _ff_classpoly_inv_disambiguate (H, Hinv); else return 1; }

#endif
