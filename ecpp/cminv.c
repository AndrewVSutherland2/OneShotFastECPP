#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "class_inv.h"          // INV_* codes, inv_string()
#include "fproot.h"
#include "invj.h"
#include "cminv.h"

/* ---- pure-mpz hardwired invariant->j formulas ----------------------------- *
 * Ported verbatim from A. V. Sutherland's classpoly class_inv_mpz.c (GPL).    */

static void cm_j_from_gamma2 (mpz_t J, mpz_t G2, mpz_t P)
{ mpz_t X; mpz_init (X); mpz_mul (X,G2,G2); mpz_mod (X,X,P); mpz_mul (J,X,G2); mpz_mod (J,J,P); mpz_clear (X); }

static void cm_j_from_f8 (mpz_t J, mpz_t F8, mpz_t P)
{
	mpz_t X,Y; mpz_init (X); mpz_init (Y);
	mpz_invert (Y,F8,P);
	mpz_mul (X,F8,F8); mpz_mod (X,X,P);  mpz_mul (X,F8,X); mpz_mod (X,X,P);      // X = f^24
	mpz_sub_ui (X,X,16);  if ( mpz_sgn (X) < 0 ) mpz_add (X,X,P);
	mpz_mul (Y,Y,X); mpz_mod (X,Y,P);                                           // gamma2 = (f^24-16)/f^8
	cm_j_from_gamma2 (J,X,P);
	mpz_clear (X); mpz_clear (Y);
}

static void cm_j_from_f4 (mpz_t J, mpz_t F4, mpz_t P)
{ mpz_t X; mpz_init (X); mpz_mul (X,F4,F4); mpz_mod (X,X,P); cm_j_from_f8 (J,X,P); mpz_clear (X); }

static void cm_j_from_f3 (mpz_t J, mpz_t F3, mpz_t P)
{
	mpz_t F,E; mpz_init (F); mpz_init (E);
	if ( ! mpz_congruent_ui_p (P,2,3) ) { printf ("cannot compute j from f^3 for P not 2 mod 3\n"); abort (); }
	mpz_mul_2exp (E,P,1); mpz_sub_ui (E,E,1); mpz_divexact_ui (E,E,3);          // (2p-1)/3
	mpz_powm (F,F3,E,P);                                                        // cube root of f^3
	cm_j_from_f8 (J,F,P);
	mpz_clear (F); mpz_clear (E);
}

static void cm_j_from_f2 (mpz_t J, mpz_t F2, mpz_t P)
{ mpz_t X; mpz_init (X); mpz_mul (X,F2,F2); mpz_mod (X,X,P); mpz_mul (X,X,X); mpz_mod (X,X,P); cm_j_from_f8 (J,X,P); mpz_clear (X); }

static void cm_j_from_f (mpz_t J, mpz_t F, mpz_t P)
{ mpz_t X; mpz_init (X); mpz_mul (X,F,F); mpz_mod (X,X,P); mpz_mul (X,X,X); mpz_mod (X,X,P); mpz_mul (X,X,X); mpz_mod (X,X,P); cm_j_from_f8 (J,X,P); mpz_clear (X); }

static void cm_j_from_t6 (mpz_t J, mpz_t T6, mpz_t P)
{
	mpz_t X,Y,Z; mpz_init (X); mpz_init (Y); mpz_init (Z);
	mpz_invert (Y,T6,P);
	mpz_mul (X,T6,T6); mpz_mod (X,X,P);
	mpz_mul_ui (Z,T6,6); mpz_mod (Z,Z,P);
	mpz_sub (X,X,Z); mpz_sub_ui (X,X,27); while ( mpz_sgn (X) < 0 ) mpz_add (X,X,P);
	mpz_mul (X,X,Y); mpz_mod (X,X,P);
	cm_j_from_gamma2 (J,X,P);
	mpz_clear (X); mpz_clear (Y); mpz_clear (Z);
}

static void cm_j_from_t2 (mpz_t J, mpz_t T2, mpz_t P)
{ mpz_t X; mpz_init (X); mpz_mul (X,T2,T2); mpz_mod (X,X,P); mpz_mul (X,X,T2); mpz_mod (X,X,P); cm_j_from_t6 (J,X,P); mpz_clear (X); }

static void cm_j_from_t (mpz_t J, mpz_t T, mpz_t P)
{ mpz_t X; mpz_init (X); mpz_mul (X,T,T); mpz_mod (X,X,P); mpz_mul (X,X,T); mpz_mod (X,X,P); mpz_mul (X,X,X); mpz_mod (X,X,P); cm_j_from_t6 (J,X,P); mpz_clear (X); }

static void cm_j_from_u8 (mpz_t J, mpz_t U8, mpz_t P)
{
	mpz_t T0,T1,T2,T3; mpz_init (T0); mpz_init (T1); mpz_init (T2); mpz_init (T3);
	mpz_mul (T2,U8,U8);
	mpz_set_ui (T0,1); mpz_sub (T0,T0,U8);  mpz_add (T3,T0,T2); mpz_sub (T0,T3,U8); mpz_mul (T1,T2,T0);
	mpz_invert (T0,T1,P);
	mpz_mul (T1,T3,T3); mpz_mul (T2,T1,T3); mpz_mul (T1,T0,T2); mpz_mul_ui (J,T1,256); mpz_mod (J,J,P);
	mpz_clear (T0); mpz_clear (T1); mpz_clear (T2); mpz_clear (T3);
}

static void cm_j_from_u2 (mpz_t J, mpz_t U, mpz_t P)
{ mpz_t X; mpz_init (X); mpz_mul (X,U,U); mpz_mod (X,X,P); mpz_mul (X,X,X); mpz_mod (X,X,P); cm_j_from_u8 (J,X,P); mpz_clear (X); }

static void cm_j_from_u (mpz_t J, mpz_t U, mpz_t P)
{ mpz_t X; mpz_init (X); mpz_mul (X,U,U); mpz_mod (X,X,P); mpz_mul (X,X,X); mpz_mod (X,X,P); mpz_mul (X,X,X); mpz_mod (X,X,P); cm_j_from_u8 (J,X,P); mpz_clear (X); }

// power transform on the raw invariant value before applying Phi_inv (etas)
static void cm_inv_power (mpz_t X, mpz_t W, mpz_t P, int inv)
{
	switch (inv) {
	case INV_W2W3E1:
		mpz_mul (X,W,W); mpz_mod (X,X,P); mpz_mul (X,X,W); mpz_mod (X,X,P); mpz_mul (X,X,X); mpz_mod (X,X,P); mpz_mul (X,X,X); mpz_mod (X,X,P); return;   // twelfth
	case INV_W3E2: case INV_W2W11E2: case INV_W2W3E2: case INV_W2W5E1: case INV_W3W3E1: case INV_W3W11E2:
		mpz_mul (X,W,W); mpz_mod (X,X,P); mpz_mul (X,X,W); mpz_mod (X,X,P); mpz_mul (X,X,X); mpz_mod (X,X,P); return;                                    // sixth
	case INV_W2W7E1:
		mpz_mul (X,W,W); mpz_mod (X,X,P); mpz_mul (X,X,X); mpz_mod (X,X,P); return;                                                                    // fourth
	case INV_W5E2: case INV_W3W5E1: case INV_W2W5E2: case INV_W3W3E2: case INV_W5W5E1: case INV_W2W17E1:
		mpz_mul (X,W,W); mpz_mod (X,X,P); mpz_mul (X,X,W); mpz_mod (X,X,P); return;                                                                    // cube
	case INV_W7E2: case INV_W2W7E2: case INV_W2W13E1: case INV_W3W7E1:
		mpz_mul (X,W,W); mpz_mod (X,X,P); return;                                                                                                      // square
	default: mpz_set (X,W); return;
	}
}

/* ---- modular-polynomial families (Atkin, single/double eta) --------------- *
 * X' = inv_power(X); j = root of Phi_inv(X', J) over F_p, via invj+fproot.     */
static int cm_bipoly (mpz_t *J, int maxj, mpz_t X, mpz_t P, int inv, const char *phidir)
{
	mpz_t Xp;  mpz_init (Xp);  cm_inv_power (Xp, X, P, inv);
	fp_ctx C;  fp_init (&C, P);
	mp_limb_t Xm[64];  fp_set_mpz (&C, Xm, Xp);
	bipoly Phi;
	if ( ! invj_load (&Phi, phidir, inv_string (inv)) ) { mpz_clear (Xp); fp_clear (&C); return -1; }
	mp_limb_t jr[8*64];
	int cap = maxj < 8 ? maxj : 8;
	int nj = invj_jroots (&C, &Phi, Xm, jr, cap, 0x5eed17 + (unsigned) inv);
	for ( int i = 0 ; i < nj ; i++ ) fp_get_mpz (&C, J[i], jr + (size_t) i * C.s);
	invj_clear (&Phi);  mpz_clear (Xp);  fp_clear (&C);
	return nj;
}

int cm_j_from_inv (mpz_t *J, int maxj, mpz_t X, mpz_t P, int inv, const char *phidir)
{
	switch (inv) {
	case INV_J:  mpz_set (J[0], X);            return 1;
	case INV_G2: cm_j_from_gamma2 (J[0],X,P);  return 1;
	case INV_F:  cm_j_from_f  (J[0],X,P);       return 1;
	case INV_F8: cm_j_from_f8 (J[0],X,P);       return 1;
	case INV_F4: cm_j_from_f4 (J[0],X,P);       return 1;
	case INV_F2: cm_j_from_f2 (J[0],X,P);       return 1;
	case INV_F3: cm_j_from_f3 (J[0],X,P);       return 1;
	case INV_T6: cm_j_from_t6 (J[0],X,P);       return 1;
	case INV_T2: cm_j_from_t2 (J[0],X,P);       return 1;
	case INV_T:  cm_j_from_t  (J[0],X,P);       return 1;
	case INV_U8: cm_j_from_u8 (J[0],X,P);       return 1;
	case INV_U2: cm_j_from_u2 (J[0],X,P);       return 1;
	case INV_U:  cm_j_from_u  (J[0],X,P);       return 1;
	default:     return cm_bipoly (J, maxj, X, P, inv, phidir);   // Atkin / single+double eta
	}
}
