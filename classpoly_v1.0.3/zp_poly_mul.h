/*

zp_poly_mul.h: exported private parameters and datatypes used by zp_poly_mul.c
users of zp_poly should not rely on any of the definitions contained herein.

Copyright (C) 2010, David Harvey and Andrew V. Sutherland

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

#ifndef ZP_POLY_MUL_H
#define ZP_POLY_MUL_H


#include <gmp.h>
#include "zp_poly_util.h"

// crude crossover for Strassen matrix mult
// use Strassen whenever log2(p) >= MIN_PBITS
// the degree doesn't seem to matter much
// set to -1 to completely disable automatic use of Strassen
#define ZP_POLY_STRASSEN_MIN_PBITS		1000

#define ZP_POLY_TFT_MIN_PBITS			500			// use naive KS for small moduli
#define ZP_POLY_TFT_MAX_DEG_PBIT_RATIO	12			// switch to naive KS when degree is large relative to pbits
//#define ZP_POLY_TFT_MIN_PBITS			0			// never switch to KS for small moduli
//#define ZP_POLY_TFT_MAX_DEG_PBIT_RATIO	100000000000	// never switch to naive KS when degree is large relative to pbits

// define for better performance when using the KGZ patch
#define ZP_POLY_USE_KGZ	0

#if ZP_POLY_USE_KGZ

extern int
mpn_mul_fft_int (mp_ptr op, mp_size_t pl,
                 mp_srcptr n, mp_size_t nl,
                 mp_srcptr m, mp_size_t ml,
                 int k);

// argggh name mangling??

#define mpn_fft_best_k __gmpn_fft_best_k

extern int
mpn_fft_best_k (mp_size_t n, int sqr_b);


// taken from core 2 tuning file
#define ZP_POLY_KGZ_MUL_THRESHOLD 238
#define ZP_POLY_KGZ_SQR_THRESHOLD 109


#endif


/*
  fermat_t = pointer to an integer modulo B^m + 1.
  Stored as array of length m+2, represents the value
     B^x[0] * (x[1] + B*x[2] + ... + B^(m-1)*x[m] + B^m*x[m+1]).
  Always semi-normalised so that 0 <= x[m+1] <= 1.
  The first value x[0] is called the bias, must have 0 <= x[0] < 2*m.

  We take sqrt(2) = B^(m/4) - B^(3m/4), and assume m is EVEN.
  (Only place this assumption is used is in the multiply-by-(1 - B^(m/2))
  code in mul_sqrt2exp().)
*/

typedef mp_limb_t* fermat_t;


typedef struct _fermat_vec_struct
{
  int n;            // length of x
  mp_limb_t* buf;   // enough space for n+1 coefficients
  fermat_t* x;      // pointers into buf (can get permuted)
  fermat_t* x2;     // temp buffer used in TFT code
  fermat_t t;       // one temp buffer (can get permuted with those in x)
  size_t m;         // as in fermat_t
  int best_k[2];    // for mpn_mul_fft_int (entries for mul and sqr)
} fermat_vec_struct;

typedef fermat_vec_struct fermat_vec_t[1];


/*
  A vector of fermat_t's
*/
typedef struct _zp_poly_mod_struct
{
  int algo;
  int d;
  mpz_t* g;
  mpz_t p;

  // algo == STUPID
  mpz_t* ginv;    // length d, ginv[d-1] = 1, (g * ginv) >> d = x^(d-1)

  // algo == BARRETT, MONTGOMERY
  int k, K, m, order;            // FFT parameters
  fermat_vec_struct* g_tft;      // transform of g

  // algo == BARRETT
  fermat_vec_struct* ginv_tft;   // transform of ginv

  // algo == MONTGOMERY
  fermat_vec_struct* R_tft;      // transform of R
  fermat_vec_struct* S_tft;      // transform of S
  mpz_t* Q_mod_g;                // Q mod g

  // algo == KRONECKER
  long cbits;				// make this a long because d*cbits may not fit in an int
  mpz_t g_ks, ginv_ks;		// g_ks, ginv_ks are (big) integers representing g and ginv
  mpz_t N1, N2, N3;			// large work variables (polys)
  mpz_t n1, n2;			// small work variables (coefficients)
}
zp_poly_mod_struct;

/*
  h =  ((h mod x^|m|) + sign(m) * f1 * f2) mod x^n, where sign(0) = 1.
  use m=0 for h = f1 * f2 mod x^n, m<0 for submul, and m>0 for addmul
  n1 and n2 are the lengths of f1 and f2 respectively.
  faster squaring happens automatically if f1 == f2 and d1 == d2

  As in all ***_mod_xn functions below, the output length is n and will be
  zero padded if necessary.

  WARNING: performance will become HORRIBLE if degree is much bigger than
  coefficient size in bits (when algo is BARRETT or MONTGOMERY).
  You have been warned.  Using AUTO will automatically switch to KRONECKER when needed.
*/
void zp_array_mul_mod_xn (mpz_t h[], mpz_t f1[], int n1, mpz_t f2[], int n2,
                          int m, int n, mpz_t p);

// h = f1 * f2, output length n1+n2-1
static inline
void zp_array_mul (mpz_t h[], mpz_t f1[], int n1, mpz_t f2[], int n2, mpz_t p)
  { zp_array_mul_mod_xn (h, f1, n1, f2, n2, 0, n1+n2-1, p); }

// h1 = f0 * f1 mod x^n,  h2 = f0 * f2 mod x^n (zero padded to length n)
  void zp_array_mul12_mod_xn (mpz_t h1[], mpz_t h2[], mpz_t f0[], int d0,
                            mpz_t f1[], int d1, mpz_t f2[], int d2, int n, mpz_t p);
  
// h = f1*f2 + f3*f4 mod x^n, clever about squaring, and zeros
// output length is zero padded to n
void zp_array_mul2_mod_xn (mpz_t h[], mpz_t f1[], int n1, mpz_t f2[], int n2,
			   mpz_t f3[], int n3, mpz_t f4[], int n4, int n, mpz_t p);

// h = f1*f2 + f3*f4
static inline
void zp_array_mul2 (mpz_t h[], mpz_t f1[], int n1, mpz_t f2[], int n2,
		    mpz_t f3[], int n3, mpz_t f4[], int n4, mpz_t p)
  { zp_array_mul2_mod_xn (h, f1, n1, f2, n2, f3, n3, f4, n4, zp_imax(n1+n2-1,n3+n4-1), p); }

// c0 = a00*b0 + a01*b1 mod x^n, c1 = a10*b0 + a11*b1 mod x^n
void zp_array_mul42_mod_xn (mpz_t c0[], mpz_t c1[],
  mpz_t a00[], int d00, mpz_t a01[], int d01, mpz_t a10[], int d10, mpz_t a11[], int d11,
  mpz_t b0[], int e0, mpz_t b1[], int e1, int n, mpz_t p);
  
// c0 = a00*b0 + a01*b1 and c1 = a10*b0 + a11*b1 of degree known to be at most n
// output is zero padded to length n
void zp_array_mul42n (mpz_t c0[], mpz_t c1[],
  mpz_t a00[], int d00, mpz_t a01[], int d01, mpz_t a10[], int d10, mpz_t a11[], int d11,
  mpz_t b0[], int e0, mpz_t b1[], int e1, int n, mpz_t p);
    
// C = A*B mod x^n where A and B are 2x2 poly matrices
// A has entries aij of deg dij and B has entried bij of deg eij
// C may alias A or B
void zp_array_mul4_mod_xn (mpz_t c00[], mpz_t c01[], mpz_t c10[], mpz_t c11[],
  mpz_t a00[], int d00, mpz_t a01[], int d01, mpz_t a10[], int d10, mpz_t a11[], int d11,
  mpz_t b00[], int e00, mpz_t b01[], int e01, mpz_t b10[], int e10, mpz_t b11[], int e11,
  int n, mpz_t p);
  
//same as above but using Strassen to trade 1 mul for 11 adds
void zp_array_mul4_mod_xn_strassen (mpz_t c00[], mpz_t c01[], mpz_t c10[], mpz_t c11[],
  mpz_t a00[], int d00, mpz_t a01[], int d01, mpz_t a10[], int d10, mpz_t a11[], int d11,
  mpz_t b00[], int e00, mpz_t b01[], int e01, mpz_t b10[], int e10, mpz_t b11[], int e11,
  int n, mpz_t p);
  

// h = (h mod x^|m|) + sign(m)*f) mod x^n, where sign(0)=1 (zero fill to n)
static inline void zp_array_addsub_modp (mpz_t h[], mpz_t f[], int d, int m, int n, mpz_t p)
{
  int i, e;

  i = 0;
  if (d > n) d = n;
  if (m > 0)
    {
      e = zp_imin (m, d);
      for (; i < e; i++) mpz_add_mod (h[i],h[i],f[i], p);
    }
  else if (m < 0)
    {
      e = zp_imin (-m, d);
      for (; i < e; i++) mpz_sub_mod (h[i], h[i], f[i], p);
      for (; i < d; i++) mpz_neg(h[i],f[i]);
    }
  if (h != f) for (; i < d; i++) mpz_set(h[i],f[i]); else i = d;
  if (i < m) i = m;
  if (i < -m ) i = m;
  for (; i < n; i++) mpz_set_ui(h[i], 0);
}

#endif
