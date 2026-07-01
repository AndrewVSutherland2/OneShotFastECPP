#include <gmp.h>
#include "class_inv.h"
#include "zp_poly.h"

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

extern bipoly_mpz_t *Phi_fj;
extern int phi_fj_inv, phi_fj_terms, phi_fj_fdeg, phi_fj_jdeg;

static inline void mpz_inv_power (mpz_t X, mpz_t W, mpz_t P, int inv)
{
	switch (inv) {
	case INV_W2W3E1: mpz_mul(X,W,W); mpz_mod(X,X,P);  mpz_mul(X,X,W); mpz_mod(X,X,P); mpz_mul(X,X,X); mpz_mod(X,X,P); mpz_mul(X,X,X); mpz_mod(X,X,P); return;					// twelfth
	case INV_W3E2: case INV_W2W11E2: case INV_W2W3E2: case INV_W2W5E1: case INV_W3W3E1: case INV_W3W11E2:
			mpz_mul(X,W,W); mpz_mod(X,X,P);  mpz_mul(X,X,W); mpz_mod(X,X,P); mpz_mul(X,X,X); mpz_mod(X,X,P); return;													// sixth
	case INV_W2W7E1: mpz_mul(X,W,W); mpz_mod(X,X,P); mpz_mul(X,X,X); mpz_mod(X,X,P); return;																			// fourth
	case INV_W5E2: case INV_W3W5E1: case INV_W2W5E2: case INV_W3W3E2: case INV_W5W5E1: case INV_W2W17E1:
			mpz_mul(X,W,W); mpz_mod(X,X,P);  mpz_mul(X,X,W); mpz_mod(X,X,P); return;																				// cube
	case INV_W7E2: case INV_W2W7E2: case INV_W2W13E1:  case INV_W3W7E1: mpz_mul(X,W,W); mpz_mod(X,X,P); return;														// square
	default: mpz_set(X,W); return;
	}
}



void mpz_j_from_gamma2 (mpz_t J, mpz_t G2, mpz_t P)
{
	static mpz_t X;
	static int init;

	if ( ! init ) { mpz_init(X); init = 1; }
	mpz_mul(X,G2,G2); mpz_mod(X,X,P); mpz_mul(J,X,G2); mpz_mod(J,J,P);
}


void mpz_j_from_f8 (mpz_t J, mpz_t F8, mpz_t P)
{
	static mpz_t X,Y;
	static int init;

	if ( ! init ) { mpz_init(X); mpz_init(Y); init = 1; }
	mpz_invert(Y,F8,P);
	mpz_mul(X,F8,F8); mpz_mod(X,X,P);  mpz_mul(X,F8,X); mpz_mod(X,X,P);				// cube f^8 to get X=f^24
	mpz_sub_ui(X,X,16);
	if ( mpz_sgn(X) < 0 ) mpz_add(X,X,P);
	mpz_mul(Y,Y,X); mpz_mod(X,Y,P);												// gamma2 = (f^24-16)/f^8
	mpz_j_from_gamma2(J,X,P);
}

void mpz_j_from_f4 (mpz_t J, mpz_t F4, mpz_t P)
{
	static mpz_t X;
	static int init;

	if (! init ) { mpz_init(X); init  = 1; }
	mpz_mul(X,F4,F4); mpz_mod(X,X,P);
	mpz_j_from_f8 (J,X,P);
}

void mpz_j_from_f3 (mpz_t J, mpz_t F3, mpz_t P)
{
	static mpz_t F, E;
	static int init;

	if (! init ) { mpz_init(E); mpz_init(F); init  = 1; }
	if ( ! mpz_congruent_ui_p (P,2,3) ) {printf ("Unable to compute j from invariant f^3 for P not 2 mod 3\n"); abort (); }
	mpz_mul_2exp (E,P,1); mpz_sub_ui(E,E,1);  mpz_divexact_ui (E,E,3);	// E = (2p-1)/3
	mpz_powm (F,F3,E,P);										// compute F = cube root of F^3
	mpz_j_from_f8 (J,F,P);
}

void mpz_j_from_f2 (mpz_t J, mpz_t F2, mpz_t P)
{
	static mpz_t X;
	static int init;

	if (! init ) { mpz_init(X); init  = 1; }
	mpz_mul(X,F2,F2); mpz_mod(X,X,P); mpz_mul(X,X,X); mpz_mod(X,X,P);
	mpz_j_from_f8 (J,X,P);
}

void mpz_j_from_f (mpz_t J, mpz_t F, mpz_t P)
{
	static mpz_t X;
	static int init;

	if (! init ) { mpz_init(X); init  = 1; }
	mpz_mul(X,F,F); mpz_mod(X,X,P); mpz_mul(X,X,X); mpz_mod(X,X,P); mpz_mul(X,X,X); mpz_mod(X,X,P);
	mpz_j_from_f8 (J,X,P);
}

void mpz_j_from_t6 (mpz_t J, mpz_t T6, mpz_t P)
{
	static mpz_t X,Y,Z;
	static int init;

	if ( ! init ) { mpz_init(X); mpz_init(Y); mpz_init(Z); init = 1; }
	mpz_invert(Y,T6,P);
	mpz_mul(X,T6,T6);  mpz_mod(X,X,P);
	mpz_mul_ui(Z,T6,6); mpz_mod(Z,Z,P);
	mpz_sub(X,X,Z); mpz_sub_ui(X,X,27); while ( mpz_sgn(X) < 0 ) mpz_add(X,X,P);
	mpz_mul(X,X,Y); mpz_mod(X,X,P);
	mpz_j_from_gamma2(J,X,P);
}

void mpz_j_from_t2 (mpz_t J, mpz_t T2, mpz_t P)
{
	static mpz_t X;
	static int init;

	if (! init ) { mpz_init(X); init  = 1; }
	mpz_mul(X,T2,T2); mpz_mod(X,X,P); mpz_mul(X,X,T2); mpz_mod(X,X,P);
	mpz_j_from_t6 (J,X,P);
}

void mpz_j_from_t (mpz_t J, mpz_t T, mpz_t P)
{
	static mpz_t X;
	static int init;

	if (! init ) { mpz_init(X); init  = 1; }
	mpz_mul(X,T,T); mpz_mod(X,X,P); mpz_mul(X,X,T); mpz_mod(X,X,P);  mpz_mul(X,X,X);  mpz_mod(X,X,P);
	mpz_j_from_t6 (J,X,P);
}

void mpz_j_from_u8 (mpz_t J, mpz_t U8, mpz_t P)
{
	static mpz_t T0,T1,T2,T3;
	static int init;

	if ( ! init ) { mpz_init(T0); mpz_init(T1); mpz_init(T2); init = 1; }
	mpz_mul(T2,U8,U8);
	mpz_set_ui(T0,1); mpz_sub(T0,T0,U8);  mpz_add(T3,T0,T2); mpz_sub(T0,T3,U8); mpz_mul(T1,T2,T0);
	mpz_invert(T0,T1,P);
	mpz_mul(T1,T3,T3); mpz_mul(T2,T1,T3); mpz_mul(T1,T0,T2); mpz_mul_ui(J,T1,256);
}

void mpz_j_from_u2 (mpz_t J, mpz_t U, mpz_t P)
{
	static mpz_t X;
	static int init;

	if (! init ) { mpz_init(X); init  = 1; }
	mpz_mul(X,U,U); mpz_mod(X,X,P); mpz_mul(X,X,X); mpz_mod(X,X,P);
	mpz_j_from_u8 (J,X,P);
}

void mpz_j_from_u (mpz_t J, mpz_t U, mpz_t P)
{
	static mpz_t X;
	static int init;

	if (! init ) { mpz_init(X); init  = 1; }
	mpz_mul(X,U,U); mpz_mod(X,X,P); mpz_mul(X,X,X); mpz_mod(X,X,P);  mpz_mul(X,X,X);  mpz_mod(X,X,P);
	mpz_j_from_u8 (J,X,P);
}

void _mpz_j_from_inv (mpz_t J, mpz_t A, mpz_t P, int inv)
{
	zp_mod_t mod;
	mpz_t x, *f;
	int d;

	inv_load_phi_fj (inv);
	mpz_init (x);
	f = zp_poly_alloc (phi_fj_jdeg, P);
	mpz_inv_power (x, A, P, inv);
	bipoly_eval_mod_mpz (f, phi_fj_jdeg, Phi_fj, phi_fj_terms, 0, x, P);
	d = zp_poly_degree (f, phi_fj_jdeg);
	zp_poly_make_monic (f, d, x, P);
	mpz_clear (x);
	zp_mod_init (mod, P);
	if ( ! zp_poly_find_root (J, f, d, mod) ) { gmp_printf ("No roots of Phi_fj(A,Y) mod P for A=%Zd\nP=%Zd\n", A, P); abort (); }
	zp_mod_clear (mod);
	zp_poly_free (f, phi_fj_jdeg);
}


void mpz_j_from_inv (mpz_t J, mpz_t X, mpz_t P, int inv)
{
	switch (inv) {
	case INV_J: mpz_set(J,X); break;
	case INV_G2: mpz_j_from_gamma2 (J,X,P); break;
	case INV_F: mpz_j_from_f (J,X,P); break;
	case INV_F8: mpz_j_from_f8 (J,X,P); break;
	case INV_F4: mpz_j_from_f4 (J,X,P); break;
	case INV_F2: mpz_j_from_f2 (J,X,P); break;
	case INV_T6: mpz_j_from_t6 (J,X,P); break;
	case INV_T2: mpz_j_from_t2 (J,X,P); break;
	case INV_F3: mpz_j_from_f3 (J,X,P); break;
	case INV_T: mpz_j_from_t (J,X,P); break;
	case INV_U8: mpz_j_from_u8 (J,X,P); break;
	case INV_U2: mpz_j_from_u2 (J,X,P); break;
	case INV_U: mpz_j_from_u (J,X,P); break;
	case INV_W3E2: case INV_W5E2: case INV_W7E2: _mpz_j_from_inv (J,X,P,inv); break;
	case INV_W3W5E1: case INV_W2W3E2: case INV_W2W5E2: case INV_W2W7E2:  case INV_W3W3E2: case INV_W3W3E1: case INV_W5W5E1:
	case INV_W2W3E1: case INV_W2W5E1: case INV_W2W7E1: case INV_W2W11E2: case INV_W2W13E1: case INV_W2W17E1: case INV_W3W7E1: case INV_W3W11E2: _mpz_j_from_inv (J,X,P,inv); break;
	default:
		if ( inv >= INV_ATKIN && inv <= INV_ATKIN_END ) {
			_mpz_j_from_inv (J, X, P, inv);
		} else if ( inv >= INV_DOUBLE_ETA && inv <= INV_DOUBLE_ETA_END ) {
			_mpz_j_from_inv (J, X, P, inv);
		} else {
			printf("Don't know how to handle invariant I=%d\n", inv); exit(0);
		}
	}
}
