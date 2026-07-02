/*

zp_poly_mul.c: multiplication module for zp_poly

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


#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <gmp.h>
#include "zp_poly_util.h"
#include "zp_poly_mul.h"
#include "zp_poly_main.h"


int bit_reverse (int x, int k)
{
  return k ? (bit_reverse (x >> 1, k - 1) | ((x & 1) << (k - 1))) : 0;
}


// ==========================================================================
// fermat_t stuff, all PRIVATE


void fermat_swap (fermat_t* x, fermat_t* y)
{
  fermat_t z;
  z = *x;
  *x = *y;
  *y = z;
}


// y = x
void fermat_set (fermat_t y, fermat_t x, size_t m)
{
  size_t i;
  for (i = 0; i <= m + 1; i++)
    y[i] = x[i];
}


// x = 0
void fermat_zero (fermat_t x, size_t m)
{
  size_t i;
  for (i = 0; i <= m + 1; i++)
    x[i] = 0;
}


// reduce x to be in [0, B^m]
// i.e. make x[m+1] = 1 only if x == -1 mod B^m + 1, else 0
void fermat_normalise (fermat_t x, size_t m)
{
  if (x[m+1])
    {
      x[m+1] = 0;
      if (mpn_sub_1 (x + 1, x + 1, m, 1))
    {
      // oops, reduced too far
      size_t i;
      for (i = 1; i <= m; i++)
        x[i] = 0;
      x[m+1] = 1;
    }
    }
}


// z = x + B^b * y
// PRETENDS that all bias fields are zero
// must have 0 <= b < 2*m
// z may alias x, but NOT y
void fermat_add_raw (fermat_t z, fermat_t x, fermat_t y, size_t b, size_t m)
{
  mp_limb_t cy;

  assert (b < 2*m);
  assert (z != y);

  if (b == 0)
    {
      mpn_add_n (z + 1, x + 1, y + 1, m + 1);
    }
  else if (b < m)
    {
      z[m+1] = x[m+1] + mpn_add_n (z + b + 1, x + b + 1, y + 1, m - b);
      cy = y[m+1] + mpn_sub_n (z + 1, x + 1, y + m - b + 1, b);
      mpn_sub_1 (z + b + 1, z + b + 1, m - b + 1, cy);
    }
  else if (b == m)
    {
      mpn_sub_n (z + 1, x + 1, y + 1, m + 1);
    }
  else
    {
      b -= m;
      z[m+1] = x[m+1] - mpn_sub_n (z + b + 1, x + b + 1, y + 1, m - b);
      cy = y[m+1] + mpn_add_n (z + 1, x + 1, y + m - b + 1, b);
      mpn_add_1 (z + b + 1, z + b + 1, m - b + 1, cy);
    }

  // semi-normalise
  if ((mp_limb_signed_t) z[m+1] < 0)
    z[m+1] = mpn_add_1 (z + 1, z + 1, m, -z[m+1]);
  else if (z[m+1] >= 2)
    z[m+1] = 1 - mpn_sub_1 (z + 1, z + 1, m, z[m+1] - 1);
}


// z = x + y
// z may alias x, but NOT y
void fermat_add (fermat_t z, fermat_t x, fermat_t y, size_t m)
{
  ptrdiff_t b;

  assert (z != y);

  b = (ptrdiff_t) y[0] - (ptrdiff_t) x[0];
  if (b < 0)
    b += 2*m;
  fermat_add_raw (z, x, y, b, m);
  z[0] = x[0];
}


// z = x - y
// z may alias x, but NOT y
void fermat_sub (fermat_t z, fermat_t x, fermat_t y, size_t m)
{
  ptrdiff_t b;

  assert (z != y);

  b = (ptrdiff_t) y[0] - (ptrdiff_t) x[0] - m;
  if (b < 0)
    b += 2*m;
  if (b < 0)
    b += 2*m;
  fermat_add_raw (z, x, y, b, m);
  z[0] = x[0];
}


// y = x, but with bias changed to zero
// x may not alias y
void fermat_remove_bias (fermat_t y, fermat_t x, size_t m)
{
  size_t i, b;

  assert (x != y);

  b = x[0];
  y[0] = 0;

  if (b == 0)
    {
      fermat_set (y, x, m);
    }
  else if (b < m)
    {
      for (i = 0; i < b; i++)
    y[i + 1] = ~x[i + m - b + 1];
      for (i = 0; i < m - b; i++)
    y[i + b + 1] = x[i + 1];
      y[m+1] = -mpn_sub_1 (y + b + 1, y + b + 1, m - b, x[m+1] + 1) - 1;
    }
  else if (b == m)
    {
      for (i = 0; i <= m; i++)
    y[i + 1] = ~x[i + 1];
      y[m+1] -= 1;
    }
  else
    {
      b -= m;
      for (i = 0; i < b; i++)
        y[i + 1] = x[i + m - b + 1];
      for (i = 0; i < m - b; i++)
    y[i + b + 1] = ~x[i + 1];
      y[m+1] = mpn_add_1 (y + b + 1, y + b + 1, m - b, x[m+1] + 1) - 1;
    }

  // semi-normalise
  if ((mp_limb_signed_t) y[m+1] < 0)
    y[m+1] = mpn_add_1 (y + 1, y + 1, m, -y[m+1]);
  else if (y[m+1] >= 2)
    y[m+1] = 1 - mpn_sub_1 (y + 1, y + 1, m, y[m+1] - 1);
}



// x *= sqrt2^e
// (possibly swapping buffers with t if required)
// 0 <= e < 4 * m * GMP_NUMB_BITS
void fermat_mul_sqrt2exp (fermat_t* x, fermat_t* t, size_t e, size_t m)
{
  int f;
  size_t r, k;
  fermat_t y;

  assert (x[0] != t[0]);
  assert (e < 4 * m * GMP_NUMB_BITS);

  // decompose sqrt2^e as (1 - B^(m/2))^f * B^k * 2^r
  // where  0 <= f <= 1,  0 <= r < GMP_NUMB_BITS,  0 <= k < 2*m
  f = e & 1;
  e >>= 1;

  if (f)
    {
      // multiply by (1 - B^(m/2)) into t
      fermat_add_raw (t[0], x[0], x[0], 3*m/2, m);
      t[0][0] = x[0][0];
      fermat_swap (x, t);

      // incorporate additional B^(m/4) into e
      e += GMP_NUMB_BITS * m / 4;
      if (e >= 2 * GMP_NUMB_BITS * m)
    e -= 2 * GMP_NUMB_BITS * m;
    }

  r = e % GMP_NUMB_BITS;
  k = e / GMP_NUMB_BITS;

  y = x[0];

  // multiply by 2^r
  if (r > 0)
    {
      mpn_lshift (y + 1, y + 1, m + 1, r);
      y[m+1] = -mpn_sub_1 (y + 1, y + 1, m, y[m+1]);
      if (y[m+1])
    y[m+1] = mpn_add_1 (y + 1, y + 1, m, 1);
    }

  // multiply by B^k
  y[0] += k;
  if (y[0] >= 2 * m)
    y[0] -= 2 * m;
}


// z = x * y
// t should be 2*m limbs scratch space
// TODO: currently this makes zp_poly_mod_mulx non-threadsafe, since it
// can modify one of the precomputed transforms. We should fix this.
// Perhaps insist that the input already be normalised, make the caller take
// care of this?
void fermat_mul (fermat_t z, fermat_t x, fermat_t y, mp_limb_t* t,
         size_t m, int* best_k)
{
  size_t i;

  fermat_normalise (x, m);
  fermat_normalise (y, m);

  if (x[m+1])
    {
      // special case: x == -1
      for (i = 0; i <= m; i++)
    z[i+1] = y[i+1];
      z[0] = x[0] + y[0] + m;
      if (z[0] >= 2 * m)
    z[0] -= 2 * m;
    }
  else if (y[m+1])
    {
      // special case: y == -1
      for (i = 0; i <= m; i++)
    z[i+1] = x[i+1];
      z[0] = x[0] + y[0] + m;
      if (z[0] >= 2 * m)
    z[0] -= 2 * m;
    }
  else
    {
#if ZP_POLY_USE_KGZ
      int sqr;
      size_t threshold;
        
      sqr = (x == y);
      threshold = sqr ? ZP_POLY_KGZ_SQR_THRESHOLD : ZP_POLY_KGZ_MUL_THRESHOLD;
      if (m >= threshold)
    {
      z[m+1] = mpn_mul_fft_int (z + 1, m, x + 1, m, y + 1, m, best_k[sqr]);
    }
      else
#endif
    {
      // multiply
      mpn_mul (t, x + 1, m, y + 1, m);

      // reduce mod B^m + 1
      z[m+1] = -mpn_sub_n (z + 1, t, t + m, m);
      if (z[m+1])
        z[m+1] = mpn_add_1 (z + 1, z + 1, m, 1);
    }

      // add biases
      z[0] = x[0] + y[0];
    }

  if (z[0] >= 2 * m)
    z[0] -= 2 * m;
}


// convert x to fermat_t.
// Assumes 0 <= x <= B^m.
void fermat_import (fermat_t y, mpz_t x, size_t m)
{
  size_t count, i;
  y[0] = 0;
  mpz_export (y + 1, &count, -1, sizeof (mp_limb_t), 0, 0, x);
  for (i = count; i <= m; i++)
    y[i + 1] = 0;
}


// convert x to mpz_t, in range 0 <= y <= B^m.
// t is scratch space
void fermat_export (mpz_t y, fermat_t x, fermat_t t, size_t m)
{
  fermat_remove_bias (t, x, m);
  fermat_normalise (t, m);
  mpz_import (y, m + 1, -1, sizeof (mp_limb_t), 0, 0, t + 1);
}


// convert x to mpz_t, in range -B^m/2 - 1 <= y < B^m/2.
// t is scratch space
void fermat_export_signed (mpz_t y, fermat_t x, fermat_t t, size_t m)
{
  size_t i;

  fermat_remove_bias (t, x, m);
  fermat_normalise (t, m);

  if (t[m+1] == 1)
    {
      // special case x == B^m
      mpz_set_si (y, -1);
      return;
    }

  if (((mp_limb_signed_t) t[m]) >= 0)
    {
      // positive
      mpz_import (y, m, -1, sizeof (mp_limb_t), 0, 0, t + 1);
    }
  else
    {
      // negative
      for (i = 0; i < m; i++)
    t[i + 1] = ~t[i + 1];
      mpz_import (y, m, -1, sizeof (mp_limb_t), 0, 0, t + 1);
      mpz_add_ui (y, y, 2);
      mpz_neg (y, y);
    }
}

// for debugging
void fermat_print (fermat_t x, size_t m)
{
  ptrdiff_t i;

  printf ("[%3d]", (int) x[0]);
  for (i = m; i >= 0; i--)
    printf (" %016lx", x[i+1]);
}


// ==========================================================================
// fermat_vec_t stuff, all PRIVATE

// limb allocation (use zp_malloc to track memory usage)
mp_limb_t *mp_limb_vec_alloc (int n)
  { return (mp_limb_t *) zp_malloc (sizeof (mp_limb_t) * n); }
void mp_limb_vec_free (mp_limb_t *t, int n)
  { zp_free (t, sizeof (mp_limb_t) * n); }


/*
  NOTE: all the FFT routines accept a `order' parameter. The i-th bit
  (starting at the LSB) encodes the type of butterfly used at recursion
  level i (top level is i == 0). If the transform length is 2^k, bits beyond
  index k are ignored.

  If the bit is 0, the usual butterfly is used, i.e.

  (x, y) -> (x + y, zeta * (x - y)).

  If the bit is 1, the "reversed" butterfly is used, i.e.

  (x, y) -> (zeta * (x - y), x + y).
*/

void fermat_vec_init (fermat_vec_t vec, int n, size_t m)
{
  int i;

  assert (m % 2 == 0);

  vec->m = m;
  vec->n = n;

#if ZP_POLY_USE_KGZ
  vec->best_k[0] = mpn_fft_best_k (m, 0);
  vec->best_k[1] = mpn_fft_best_k (m, 1);
#endif

  vec->x = (fermat_t*) zp_malloc (sizeof (fermat_t) * n);
  vec->x2 = (fermat_t*) zp_malloc (sizeof (fermat_t) * n);
  vec->buf = mp_limb_vec_alloc ((m + 2) * (n + 1));
  vec->t = vec->buf;
  for (i = 0; i < n; i++)
    vec->x[i] = vec->buf + (i + 1) * (m + 2);
}


void fermat_vec_clear (fermat_vec_t vec)
{
  zp_free (vec->x2, sizeof(fermat_t) * vec->n);
  zp_free (vec->x, sizeof(fermat_t) * vec->n);
  mp_limb_vec_free (vec->buf, (vec->m + 2) * (vec->n + 1));
}


// imports f[0], ... f[d-1] into this vector.
// expects all coefficients in [0, B^m].
void fermat_vec_import (fermat_vec_t vec, mpz_t f[], int d)
{
  int i;

  for (i = 0; i < d; i++)
    fermat_import (vec->x[i], f[i], vec->m);
}


// exports f[s], ... f[s+d-1] to f, reducing every coefficient mod p.
// uses fermat_export_signed if sign is nonzero, otherwise fermat_export
void fermat_vec_export_modp (mpz_t f[], fermat_vec_t vec,
                 int s, int d, int sign, mpz_t p)
{
  int i;
  mpz_t t;

  mpz_init2 (t,2*zp_bits(p)+ZP_PBITS_FUDGE);    // avoid reallocating
  for (i = 0; i < d; i++)
    {
      if (sign)
    fermat_export_signed (t, vec->x[s + i], vec->t, vec->m);
      else
    fermat_export (t, vec->x[s + i], vec->t, vec->m);
      mpz_mod (f[i], t, p);
    }
  mpz_clear (t);
}

// exports f[0], ... f[min(d,n)-1] to f, reducing mod p
// Assumes f[i] unsigned, zero pads to length n
void fermat_vec_export_xn_modp (mpz_t f[], fermat_vec_t vec,
                                int d, int n, mpz_t p)
{
  int i;
  mpz_t t;

  if (d > n) d = n; // truncate mod x^n
  mpz_init2 (t,2*zp_bits(p)+ZP_PBITS_FUDGE);    // avoid reallocating
  for (i = 0; i < d; i++)
    {
      fermat_export (t, vec->x[i], vec->t, vec->m);
      mpz_mod (f[i], t, p);
    }
  mpz_clear (t);
  for (; i < n; i++) mpz_set_ui (f[i],0);
}

// exports f[0], ... f[min(d,n)-1] to f, reducing mod p
// handles signed values, zero pads to length n
void fermat_vec_exports_xn_modp (mpz_t f[], fermat_vec_t vec,
                                int d, int n, mpz_t p)
{
  int i;
  mpz_t t;

  if (d > n) d = n; // truncate mod x^n
  mpz_init2 (t,2*zp_bits(p)+ZP_PBITS_FUDGE);    // avoid reallocating
  for (i = 0; i < d; i++)
    {
      fermat_export_signed (t, vec->x[i], vec->t, vec->m);
      mpz_mod (f[i], t, p);
    }
  mpz_clear (t);
  for (; i < n; i++) mpz_set_ui (f[i],0);
}

// exports f[0], ... f[min(d,n)-1] to f, reducing mod p
// if m is nonzero the output is added (m>0) or subtracted (m<0)
// to the current contents of f[i] (automatically zero-extended
// when m is less than d).  Assumes f[i] unsigned.
// If n is greater than d and m, output is zero padded to length n
void fermat_vec_export_addsub_xn_modp (mpz_t f[], fermat_vec_t vec,
                           int d, int m, int n, mpz_t p)
{
  int i, e;
  mpz_t t;

  i = 0;
  if (d > n) d = n;   // truncate mod x^n
  mpz_init2 (t,2*zp_bits(p)+ZP_PBITS_FUDGE);    // avoid reallocating
  if (m > 0)
    {
      e = zp_imin (m, d);
      for (; i < e; i++ )
        {
      fermat_export (t, vec->x[i], vec->t, vec->m);
      mpz_add (t, t, f[i]);  mpz_mod (f[i], t, p);
        }
    }
  else if (m < 0)
    {
      e = zp_imin (-m, d);
      for (; i < e; i++ )
        {
      fermat_export (t, vec->x[i], vec->t, vec->m);
      mpz_sub (t, t, f[i]); mpz_neg (t,t); mpz_mod (f[i], t, p);
        }
    }
  for (; i < d; i++ )
    {
      fermat_export (t, vec->x[i], vec->t, vec->m);
      if (m < 0) mpz_neg(t,t);
      mpz_mod (f[i], t, p);
    }
  mpz_clear (t);
  if (i < m) i = m;
  if (i < -m) i = -m;
  for (; i < n; i++) mpz_set_ui (f[i],0);
}

// returns integer r such that coefficient size (m) must be divisible
// by r for transforms of length K = 2^k
int get_mdiv (int K)
{
  if (K <= 4 * GMP_NUMB_BITS)
    return 2;
  return K / 4 / GMP_NUMB_BITS;
}



// FFT on first 2^k coefficients
void fermat_vec_fft (fermat_vec_t vec, int k, int order)
{
  int i, K;
  size_t w, r;

  if (k == 0)
    return;

  K = 1 << (k - 1);

  // sqrt2^w = k-th root of unity
  w = (4 * GMP_NUMB_BITS * vec->m) >> k;
  r = 0;

  // butterflies
  for (i = 0; i < K; i++, r += w)
    {
      fermat_sub (vec->t, vec->x[i], vec->x[i+K], vec->m);
      fermat_add (vec->x[i], vec->x[i], vec->x[i+K], vec->m);
      fermat_swap (&vec->t, &vec->x[i+K]);
      fermat_mul_sqrt2exp (&vec->x[i+K], &vec->t, r, vec->m);

      if (order & 1)
    fermat_swap (&vec->x[i], &vec->x[i+K]);
    }

  // recurse
  fermat_vec_fft (vec, k - 1, order >> 1);
  vec->x += K;
  fermat_vec_fft (vec, k - 1, order >> 1);
  vec->x -= K;
}



// transposes first K1*K2 entries of vec->x
// from K2 columns and K1 rows
// to K1 columns and K2 rows
// (swaps buffers with vec->x2)
void fermat_vec_transpose (fermat_vec_t vec, int K1, int K2)
{
  int i, j;
  fermat_t* temp;

  for (i = 0; i < K2; i++)
    for (j = 0; j < K1; j++)
      vec->x2[i*K1 + j] = vec->x[j*K2 + i];

  temp = vec->x;
  vec->x = vec->x2;
  vec->x2 = temp;
}



// truncated FFT
// 2^k = transform length
// z, n as in "A cache-friendly truncated FFT"
// sqrt2^y is twisting parameter (\zeta in the paper)
// must have 0 <= y < 4 * GMP_NUMB_BITS * m / 2^k
// (algorithm basically is copied out of pseudo-code from the paper)
void fermat_vec_tft (fermat_vec_t vec, int k, size_t y, int z, int n,
             int order)
{
  int i, k1, k2, K1, K2, n1, n1p, n2, z1, z2p, z2;
  size_t w;
  fermat_t* save_x;

  assert (z <= (1 << k) && 1 <= n && 1 <= (1 << k));
  assert (y < ((4 * GMP_NUMB_BITS * vec->m) >> k));

  if (z == 0)
    {
      for (i = 0; i < n; i++)
    fermat_zero (vec->x[i], vec->m);
      return;
    }

  if (k == 0)
    return;

  if (k == 1)
    {
      // base case

      if (n == 2)
    {
      if (z == 2)
        {
              fermat_sub (vec->t, vec->x[0], vec->x[1], vec->m);
          fermat_add (vec->x[0], vec->x[0], vec->x[1], vec->m);
          fermat_swap (&vec->x[1], &vec->t);
          fermat_mul_sqrt2exp (&vec->x[1], &vec->t, y, vec->m);
          if (order & 1)
        fermat_swap (&vec->x[0], &vec->x[1]);
        }
      else  // z == 1
        {
          fermat_set (vec->x[1], vec->x[0], vec->m);
          if (order & 1)
        fermat_mul_sqrt2exp (&vec->x[0], &vec->t, y, vec->m);
          else
        fermat_mul_sqrt2exp (&vec->x[1], &vec->t, y, vec->m);
        }
    }
      else  // n == 1
    {
      if (order & 1)
        {
          if (z == 2)
        fermat_sub (vec->x[0], vec->x[0], vec->x[1], vec->m);
          fermat_mul_sqrt2exp (&vec->x[0], &vec->t, y, vec->m);
        }
      else
        {
          if (z == 2)
        fermat_add (vec->x[0], vec->x[0], vec->x[1], vec->m);
        }
    }

      return;
    }

  // recursive case

  k1 = k / 2;
  k2 = k - k1;
  K1 = 1 << k1;
  K2 = 1 << k2;

  n2 = n & (K2 - 1);
  n1 = n >> k2;
  n1p = ((n - 1) >> k2) + 1;
  z2 = z & (K2 - 1);
  z1 = z >> k2;
  z2p = z1 ? K2 : z2;

  w = (4 * GMP_NUMB_BITS * vec->m) >> k;

  fermat_vec_transpose (vec, K1, K2);

  // column transforms
  save_x = vec->x;
  for (i = 0; i < z2; i++, vec->x += K1)
    fermat_vec_tft (vec, k1, y + w * i, z1 + 1, n1p, order);
  for (; i < z2p; i++, vec->x += K1)
    fermat_vec_tft (vec, k1, y + w * i, z1, n1p, order);
  vec->x = save_x;

  fermat_vec_transpose (vec, K2, K1);

  // row transforms
  save_x = vec->x;
  for (i = 0; i < n1; i++, vec->x += K2)
    fermat_vec_tft (vec, k2, y << k1, z2p, K2, order >> k1);
  if (n2)
    fermat_vec_tft (vec, k2, y << k1, z2p, n2, order >> k1);
  vec->x = save_x;
}



// inverse FFT on first 2^k coefficients
void fermat_vec_ifft (fermat_vec_t vec, int k, int order)
{
  int i, K;
  size_t w, r;

  if (k == 0)
    return;

  K = 1 << (k - 1);

  // recurse
  fermat_vec_ifft (vec, k - 1, order >> 1);
  vec->x += K;
  fermat_vec_ifft (vec, k - 1, order >> 1);
  vec->x -= K;

  // sqrt2^w = k-th root of unity
  w = (4 * GMP_NUMB_BITS * vec->m) >> k;
  r = 2 * GMP_NUMB_BITS * vec->m;

  // butterflies
  for (i = 0; i < K; i++, r -= w)
    {
      if (order & 1)
    fermat_swap (&vec->x[i], &vec->x[i+K]);

      fermat_mul_sqrt2exp (&vec->x[i+K], &vec->t, r, vec->m);
      fermat_add (vec->t, vec->x[i], vec->x[i+K], vec->m);
      fermat_sub (vec->x[i], vec->x[i], vec->x[i+K], vec->m);
      fermat_swap (&vec->t, &vec->x[i+K]);
    }
}



// truncated IFFT
// notation same as fermat_vec_tft
void fermat_vec_itft (fermat_vec_t vec, int k, size_t y, int z, int n, int f,
              int order)
{
  int i, k1, k2, K1, K2, n1, n2, z1, z2, z2p, fp, m1, m2;
  size_t w;
  fermat_t* save_x;

  assert (1 <= z && z <= (1 << k));
  assert (z >= n);
  assert (f == 0 || f == 1);
  assert (1 <= n + f && n + f <= (1 << k));
  assert (y < ((4 * GMP_NUMB_BITS * vec->m) >> k));

  if (k == 0)
    return;

  if (k == 1)
    {
      // base case

      if (n == 2)
    {
      if (order & 1)
        fermat_swap (&vec->x[0], &vec->x[1]);
      y = 2 * GMP_NUMB_BITS * vec->m - y;
      fermat_mul_sqrt2exp (&vec->x[1], &vec->t, y, vec->m);
      fermat_swap (&vec->t, &vec->x[1]);
      fermat_add (vec->x[1], vec->x[0], vec->t, vec->m);
      fermat_sub (vec->x[0], vec->x[0], vec->t, vec->m);
    }
      else if (n == 1)
    {
      if (f)
        {
          if (z == 2)
        {
          if (order & 1)
            {
              if (y > 0)
            y = 4 * GMP_NUMB_BITS * vec->m - y;
              fermat_mul_sqrt2exp (&vec->x[0], &vec->t, y, vec->m);
              fermat_add (vec->x[1], vec->x[1], vec->x[0], vec->m);
              fermat_add (vec->x[0], vec->x[0], vec->x[1], vec->m);
            }
          else
            {
              y = y + 2 * GMP_NUMB_BITS * vec->m;
              fermat_sub (vec->x[1], vec->x[1], vec->x[0], vec->m);
              fermat_sub (vec->x[0], vec->x[0], vec->x[1], vec->m);
              fermat_mul_sqrt2exp (&vec->x[1], &vec->t, y, vec->m);
            }

        }
          else  // z == 1
        {
          if (order & 1)
            {
              if (y > 0)
            y = 4 * GMP_NUMB_BITS * vec->m - y;
              fermat_mul_sqrt2exp (&vec->x[0], &vec->t, y, vec->m);
              fermat_set (vec->x[1], vec->x[0], vec->m);
              fermat_mul_sqrt2exp (&vec->x[0], &vec->t, 2, vec->m);
            }
          else
            {
              fermat_set (vec->x[1], vec->x[0], vec->m);
              fermat_mul_sqrt2exp (&vec->x[0], &vec->t, 2, vec->m);
              fermat_mul_sqrt2exp (&vec->x[1], &vec->t, y, vec->m);
            }
        }
        }
      else  // f == 0
        {
          if (z == 2)
        {
          if (order & 1)
            {
              y = 2 * GMP_NUMB_BITS * vec->m - y + 2;
              fermat_mul_sqrt2exp (&vec->x[0], &vec->t, y, vec->m);
              fermat_sub (vec->x[1], vec->x[1], vec->x[0], vec->m);
              fermat_swap (&vec->x[0], &vec->x[1]);
            }
          else
            {
              fermat_mul_sqrt2exp (&vec->x[0], &vec->t, 2, vec->m);
              fermat_sub (vec->x[0], vec->x[0], vec->x[1], vec->m);
            }
        }
          else  // z == 1
        {
          if (order & 1)
            {
              y = (y <= 2) ? (2 - y) :
            (4 * GMP_NUMB_BITS * vec->m + 2 - y);
              fermat_mul_sqrt2exp (&vec->x[0], &vec->t, y, vec->m);
            }
          else
            {
              fermat_mul_sqrt2exp (&vec->x[0], &vec->t, 2, vec->m);
            }
        }
        }
    }
      else  // n == 0
    {
      if (order & 1)
        {
          if (z == 2)
        fermat_sub (vec->x[0], vec->x[0], vec->x[1], vec->m);
          y = (y >= 2) ? (y - 2) : (4 * GMP_NUMB_BITS * vec->m + y - 2);
          fermat_mul_sqrt2exp (&vec->x[0], &vec->t, y, vec->m);
        }
      else
        {
          if (z == 2)
        fermat_add (vec->x[0], vec->x[0], vec->x[1], vec->m);
          y = 4 * GMP_NUMB_BITS * vec->m - 2;
          fermat_mul_sqrt2exp (&vec->x[0], &vec->t, y, vec->m);
        }
    }

      return;
    }

  // recursive case

  k1 = k / 2;
  k2 = k - k1;
  K1 = 1 << k1;
  K2 = 1 << k2;

  n2 = n & (K2 - 1);
  n1 = n >> k2;
  z2 = z & (K2 - 1);
  z1 = z >> k2;
  fp = (n2 + f > 0);
  z2p = z1 ? K2 : z2;
  m1 = (n2 < z2) ? n2 : z2;
  m2 = (n2 < z2) ? z2 : n2;

  w = (4 * GMP_NUMB_BITS * vec->m) >> k;

  // row transforms
  save_x = vec->x;
  for (i = 0; i < n1; i++, vec->x += K2)
    fermat_vec_itft (vec, k2, y << k1, K2, K2, 0, order >> k1);
  vec->x = save_x;

  fermat_vec_transpose (vec, K1, K2);

  // rightmost column transforms
  save_x = vec->x;
  for (i = n2, vec->x += n2 * K1; i < m2; i++, vec->x += K1)
    fermat_vec_itft (vec, k1, y + w * i, z1 + 1, n1, fp, order);
  for (; i < z2p; i++, vec->x += K1)
    fermat_vec_itft (vec, k1, y + w * i, z1, n1, fp, order);
  vec->x = save_x;

  // last row transform
  if (fp)
    {
      fermat_vec_transpose (vec, K2, K1);

      save_x = vec->x;
      vec->x += n1 * K2;
      fermat_vec_itft (vec, k2, y << k1, z2p, n2, f, order >> k1);
      vec->x = save_x;

      fermat_vec_transpose (vec, K1, K2);
    }

  // leftmost column transforms
  save_x = vec->x;
  for (i = 0; i < m1; i++, vec->x += K1)
    fermat_vec_itft (vec, k1, y + w * i, z1 + 1, n1 + 1, 0, order);
  for (; i < n2; i++, vec->x += K1)
    fermat_vec_itft (vec, k1, y + w * i, z1, n1 + 1, 0, order);
  vec->x = save_x;

  fermat_vec_transpose (vec, K2, K1);
}



// for debugging
void fermat_vec_print (fermat_vec_t vec)
{
  int i;

  for (i = 0; i < vec->n; i++)
    {
      fermat_print (vec->x[i], vec->m);
      printf ("\n");
    }
}


// ==========================================================================

// PRIVATE
// allocate array of n mpz_t's, all mpz_init'd
mpz_t* mpz_vec_alloc (int n)
{
  int i;
  mpz_t* v;

  v = (mpz_t*) zp_malloc (sizeof (mpz_t) * n);
  for (i = 0; i < n; i++)
    mpz_init (v[i]);

  return v;
}


// PRIVATE
void mpz_vec_free (mpz_t* v, int n)
{
  int i;
  for (i = 0; i < n; i++)
    mpz_clear (v[i]);
  zp_free (v, sizeof(mpz_t) * n);
}

size_t choose_m (size_t bits, int K)
{
  size_t m, mdiv;

  m = bits / GMP_NUMB_BITS + 1;
  mdiv = get_mdiv (K);

  // round up to multiple of mdiv
  m = ((m - 1) / mdiv + 1) * mdiv;

  // round up to make GMP happy
#if ZP_POLY_USE_KGZ
  int best_k;
  best_k = mpn_fft_best_k (m, 0);
  m = (((m - 1) >> best_k) + 1) << best_k;

  best_k = mpn_fft_best_k (m, 1);
  m = (((m - 1) >> best_k) + 1) << best_k;
#endif

  return m;
}


// ==========================================================================
// plain polynomial multiplication exported for drew's pleasure :-)

// a few handy helpers to make the code below slightly easier to read
static inline
int is1 (mpz_t f[], int d) { return (d==1 && mpz_cmp_ui(f[0],1)==0); }

static inline
void zp_array_copy (mpz_t h[], mpz_t f[], int d)
  { int i; if ( h != f ) for (i = 0; i < d; i++) mpz_set(h[i],f[i]); }

static inline
void zp_array_zero (mpz_t h[], int d)
  { int i; for (i = 0; i < d; i++) mpz_set_ui(h[i], 0); }

static inline
void fermat_vec_import_tft (fermat_vec_t v, mpz_t f[], int d, int k, int n, int o)
  { fermat_vec_import (v, f, d);  fermat_vec_tft (v, k, 0, d, n, o); }

// vec3 = vec1*vec2
static inline
void fermat_vec_mul (fermat_vec_t vec3, fermat_vec_t vec1, fermat_vec_t vec2,
             int d, int k, size_t m, mp_limb_t *temp)
{
  int i;

  for (i = 0; i < d; i++)
    fermat_mul (vec3->x[i], vec1->x[i], vec2->x[i], temp, m, vec3->best_k);
}

// vec3 = vec1*vec2
static inline
void fermat_vec_sqrt2exp (fermat_vec_t vec1, int d, int k, size_t m)
{
  int i;
  size_t s;

  s = (4 * GMP_NUMB_BITS * m - 2 * k) % (4 * GMP_NUMB_BITS * m);
  for (i = 0; i < d; i++)
    fermat_mul_sqrt2exp (&vec1->x[i], &vec1->t, s, m);
}

static inline
void fermat_vec_sqrt2exp_itft (fermat_vec_t v, int d, int k, size_t m, int o)
  { fermat_vec_sqrt2exp (v, d, k, m);  fermat_vec_itft (v, k, 0, d, d, 0, o); }
  
// vec3 = vec1 + vec2
// vec3 may alias vec1 but not vec2
static inline
void fermat_vec_add (fermat_vec_t vec3, fermat_vec_t vec1, fermat_vec_t vec2, int d, size_t m)
{
  int i;

  assert (vec3 != vec2);
  for (i = 0; i < d; i++) fermat_add (vec3->x[i], vec1->x[i], vec2->x[i], m);
}

// vec1 = vec1 + vec2
static inline
void fermat_vec_addto (fermat_vec_t vec1, fermat_vec_t vec2, int d, size_t m)
{
  int i;

  assert (vec1 != vec2);
  for (i = 0; i < d; i++) fermat_add (vec1->x[i], vec1->x[i], vec2->x[i], m);
}

// vec3 = vec1 - vec2
// vec3 may alias vec1 but not vec2
static inline
void fermat_vec_sub (fermat_vec_t vec3, fermat_vec_t vec1, fermat_vec_t vec2,
                     int d, size_t m)
{
  int i;
    
  assert (vec3 != vec2);
  for (i = 0; i < d; i++) fermat_sub (vec3->x[i], vec1->x[i], vec2->x[i], m);
}


// ==========================================================================
// kronecker substitution code, all PRIVATE

static inline void mpz_trim (mpz_t N) { while ( N->_mp_size && ! N->_mp_d[N->_mp_size-1] ) N->_mp_size--; }

static inline void ks_unpack ( mpz_t n, mpz_t N, int off, int len)
{
    register int i, j, k, m, e;

    // deal with cases where off+len extends beyond the length of N
    k = N->_mp_size*GMP_NUMB_BITS;
    if ( off >= k ) { n->_mp_size = 0; return; }
    if ( off+len > k ) len = k-off;
    
    j = off / GMP_NUMB_BITS;  k = off % GMP_NUMB_BITS;  e = (off + len ) / GMP_NUMB_BITS;

    // make sure n is big enough
    if ( len > n->_mp_alloc*GMP_NUMB_BITS ) mpz_realloc2 (n, len);

    if ( ! k ) {
        i = 0;
        while ( j < e ) n->_mp_d[i++] = N->_mp_d[j++];
        k = (off+len) % GMP_NUMB_BITS;
        if ( k ) { m = GMP_NUMB_BITS - k;  n->_mp_d[i++] = (N->_mp_d[j] << m) >> m; }
    } else if ( j == e ) {
        m = GMP_NUMB_BITS - ((off + len) % GMP_NUMB_BITS);
        n->_mp_d[0] = (N->_mp_d[j]<<m) >> (k+m);
        i = 1;
    } else {
        m = GMP_NUMB_BITS - k;
        n->_mp_d[0] = N->_mp_d[j++] >> k;
        for ( i = 1 ; j < e ; i++, j++ ) { n->_mp_d[i-1] |= N->_mp_d[j] << m;  n->_mp_d[i] = N->_mp_d[j] >> k; }
        e = (off+len) % GMP_NUMB_BITS;
        if ( e <= k ) {
              m = GMP_NUMB_BITS - e;
            if ( e ) n->_mp_d[i-1] |= (N->_mp_d[j] << m) >> (k-e);
        } else {
            n->_mp_d[i-1] |= N->_mp_d[j] << m;
            m = GMP_NUMB_BITS - (e-k);  e = GMP_NUMB_BITS - e;
            n->_mp_d[i++] = (N->_mp_d[j] << e) >> m;
        }
    }
    n->_mp_size = i;
    mpz_trim (n);
}

static inline void ks_pack (mpz_t N, mpz_t n, int off, int len)
{
    register int i, j, k, m, e;

//  assert ( len >= GMP_NUMB_BITS* n->_mp_size && len >= GMP_NUMB_BITS );

    j = off / GMP_NUMB_BITS;  k = off % GMP_NUMB_BITS;  m = GMP_NUMB_BITS - k;  e = (off + len ) / GMP_NUMB_BITS;

    // assume N is big enough, caller should check this
    //if ( off+len > N->_mp_alloc*GMP_NUMB_BITS ) mpz_realloc2 (N, off+len);
    
    // zero extend if needed
    // IMPORTANT: we don't want N to be trimmed during packing, otherwise packing zeros takes quadratic time!
    for ( i = N->_mp_size ; i < e ; i++ ) N->_mp_d[i] = 0;
    if ( i==e && ((off+len)%GMP_NUMB_BITS) ) N->_mp_d[i] = 0;
    
    // if off is aligned on a limb boundary, life is pretty simple
    if ( !k ) {
        for ( i = 0 ; i < n->_mp_size ; i++ ) N->_mp_d[j++] = n->_mp_d[i];
    } else if ( ! n->_mp_size ) {
        N->_mp_d[j] = ((N->_mp_d[j]<<m)>>m);  j++;
    } else {
        N->_mp_d[j] = ((N->_mp_d[j]<<m)>>m) | (n->_mp_d[0] << k);  j++;
        for ( i = 1 ; i < n->_mp_size ; i++ ) N->_mp_d[j++] = (n->_mp_d[i-1] >> m) | (n->_mp_d[i] << k);
        if ( j >= N->_mp_size || j < e ) N->_mp_d[j++] = n->_mp_d[i-1] >> m;
        else {
            k = (off+len) % GMP_NUMB_BITS;
            N->_mp_d[j] = ((N->_mp_d[j]>>k)<<k) | (n->_mp_d[i-1] >> m);  j++;
        }
    }
    // zero fill to len bits if needed
    if ( j <= e ) {
        while ( j < e ) N->_mp_d[j++] = 0;
        k = (off+len) % GMP_NUMB_BITS;
        if ( k ) { N->_mp_d[j] = (N->_mp_d[j]>>k)<<k; j++; }
    }
    // adjust size, but don't trim leading zeros(!), let caller deal with this
    if ( j > N->_mp_size ) N->_mp_size = j;
    //mpz_trim(N);
}

// pack coefficient array a of length n (mod p) into integer N
static inline void ks_import (mpz_t N, mpz_t a[], int n, int bits)
{
    register int i, off;

    // this assumption could be removed, but it is convenient and will hold for all p > 2^31
    assert ( bits >= GMP_NUMB_BITS );

    N->_mp_size = 0;
    for ( i = off = 0 ; i < n ; i++, off += bits ) if ( mpz_sgn (a[i]) ) ks_pack (N, a[i], off, bits);
    mpz_trim (N);
}

// unpack integer N into coefficient array a of length n, reducing mod p
// a nonzero value of m indicates an addition or subtraction to the existing contents of a
// effectively sets a to (a mod x^|m| + sign(m) * N) mod x^n (where sign(m) == 1 if m is 0)
static inline void ks_export_addsub_xn_modp (mpz_t a[], mpz_t N, int m, int n, int bits, mpz_t p)
{
    mpz_t w;
    register int i,am, off, len;

    off = 0;
    len = bits;
    mpz_init2 (w, len);
    am = ( m < 0 ? -m : m);
    for ( i = 0 ; i < n && off < N->_mp_size*GMP_NUMB_BITS ; i++, off += len ) {
        ks_unpack (w, N, off, len);
        if ( m < 0 ) {
            if ( i < am ) mpz_sub (w, a[i], w); else mpz_neg(w, w);
        } else if ( m > 0 ) {
            if ( i < am ) mpz_add (w, a[i], w);
        }
        mpz_mod(a[i],w,p);
    }
    mpz_clear (w);
    if ( i < am ) i = am;
    while ( i < n ) mpz_set_ui(a[i++],0);
}

// reduce each packed coefficient of N mod p (in place)
static inline void ks_reduce (mpz_t N, int n, int bits, mpz_t p)
{
    mpz_t w;
    register int i, off;

    mpz_init2 (w, bits);
    for ( i = off = 0 ; i < n && off < N->_mp_size*GMP_NUMB_BITS ; i++, off += bits ) {
        ks_unpack (w, N, off, bits);
        if ( mpz_cmp(w,p) >= 0 ) {
            mpz_mod (w, w, p);
            ks_pack (N, w, off, bits);
        }
    }
    mpz_clear (w);
    mpz_trim (N);
}

// h = ((h mod x^|m|) + sign(m) * f1 * f2) mod x^n.  sign(0)=1.
// be clever about squaring and mult by 0 or 1
// interface is identical to zp_array_mul_mod_xn below
void zp_array_mul_mod_xn_ks (mpz_t h[], mpz_t f1[], int d1, mpz_t f2[], int d2, int m, int n, mpz_t p)
{
    mpz_t v1, v2, v3;
    int i, d3, sqr;
    long cbits;

    // truncate input lengths and compute output length
    if (d1 > n) d1 = n;
    if (d2 > n) d2 = n;
    d3 = d1 + d2 - 1;

    // if f1*f2 = 0, zero-fill the output and return
    if (!d1 || !d2) { for (i = (m<0?-m:m); i < n; i++) mpz_set_ui(h[i],0);  return; }
        
    // if f1 or f2 is 1, just add or subtract
    if (is1(f1,d1)) { zp_array_addsub_modp (h, f2, d2, m, n, p); return; }
    if (is1(f2,d2)) { zp_array_addsub_modp (h, f1, d1, m, n, p); return; }

    for ( i = 1 ; (1<<i ) <= d3 ; i++ );            // actually max(d1,d2) would suffice
    cbits = i + 2*zp_bits(p) + 1;
    if ( cbits < GMP_NUMB_BITS ) cbits = GMP_NUMB_BITS;

    // import inputs into mpz's and multiply them
    mpz_init2 (v1, d1*cbits);
    ks_import (v1, f1, d1, cbits);
    mpz_init2 (v3, d3*cbits);
    
    if ( f1 != f2 || d1 != d2 ) {
        mpz_init2 (v2, d2*cbits);
        ks_import (v2, f2, d2, cbits);
        mpz_mul (v3, v1, v2);
        sqr = 1;
    } else {
        mpz_mul (v3, v1, v1);
        sqr = 0;
    }
        
    // export output, adding/subtracting and reducing mod p in the process
    ks_export_addsub_xn_modp (h, v3, m, n, cbits, p);
    
    // cleanup
    mpz_clear (v1);
    if ( sqr ) mpz_clear (v2);
    mpz_clear (v3);
}

// h1 = f0*f1 mod x^n, h2 = f0*f2 mod x^n (zero-padded to length n)
// called by zp_array_mul12_mod_xn below for suitable p and d
// h1 and h2 may alias any of the inputs
void zp_array_mul12_mod_xn_ks (mpz_t h1[], mpz_t h2[], mpz_t f0[], int d0,
                               mpz_t f1[], int d1, mpz_t f2[], int d2, int n, mpz_t p)
{
    mpz_t v, v0, v1;
    int i, d, cbits;

    d = zp_imax (d0+d1, d0+d2) - 1;
    for ( i = 1 ; (1<<i ) <= d ; i++ );         // actually max(d0,d1,d2) would suffice

    cbits = i + 2*zp_bits(p) + 1;
    if ( cbits < GMP_NUMB_BITS ) cbits = GMP_NUMB_BITS;

    mpz_init2 (v0, d0*cbits);  ks_import (v0, f0, d0, cbits);
    mpz_init2 (v1, zp_imax(d1,d2)*cbits);  ks_import (v1, f1, d1, cbits);
    mpz_init2 (v, d*cbits);
    mpz_mul (v, v0, v1);
    ks_import (v1, f2, d2, cbits);
    ks_export_addsub_xn_modp (h1, v, 0, n, cbits, p);
    mpz_mul (v, v0, v1);
    ks_export_addsub_xn_modp (h2, v, 0, n, cbits, p);
    mpz_clear (v); mpz_clear (v0); mpz_clear (v1);
}

// h = (f1 * f2 + f3 * f4) mod x^n (output zero-padded to length n)
// called by zp_array_mul2_mod_xn or zp_array_mul2n below for suitable p and d
// h may alias any input
void zp_array_mul2_mod_xn_ks (mpz_t h[], mpz_t f1[], int d1, mpz_t f2[], int d2,
                  mpz_t f3[], int d3, mpz_t f4[], int d4, int n, mpz_t p)
{
    mpz_t v, v0, v1;
    int i, d, cbits;

    d = zp_imax (d1+d2, d3+d4) - 1;
    for ( i = 1 ; (1<<i ) <= d ; i++ );         // actually max(d1,d2,d3,d4) would suffice

    cbits = i + 2*zp_bits(p) + 1;
    if ( cbits < GMP_NUMB_BITS ) cbits = GMP_NUMB_BITS;

    mpz_init2 (v0, zp_imax(d1,d3)*cbits);  ks_import (v0, f1, d1, cbits);
    mpz_init2 (v1, zp_imax(d2,d4)*cbits);  ks_import (v1, f2, d2, cbits);
    mpz_init2 (v, d*cbits);
    mpz_mul (v, v0, v1);
    ks_import (v0, f3, d3, cbits);  ks_import (v1, f4, d4, cbits);
    ks_export_addsub_xn_modp (h, v, 0, n, cbits, p);
    mpz_mul (v, v0, v1);
    ks_export_addsub_xn_modp (h, v, n, n, cbits, p);
    mpz_clear (v); mpz_clear (v0); mpz_clear (v1);
}

// c0 = (a00 * b0 + a01 * b1) mod x^n, c1 = (a10 * b0 + a11 * b1) mod x^n
// i.e. c=a*b where a is a 2x2 matrix and b and c are 2-vectors.
// called by zp_array_mul2_mod_xn and zp_array_mul42n below for suitable p and d
// c may alias b or a column of a
void zp_array_mul42_mod_xn_ks (mpz_t c0[], mpz_t c1[],
  mpz_t a00[], int d00, mpz_t a01[], int d01, mpz_t a10[], int d10, mpz_t a11[], int d11,
  mpz_t b0[], int e0, mpz_t b1[], int e1, int n, mpz_t p)
{
    mpz_t v, v0, v1, v2, v3;
    int i, d, cbits;

    d = zp_imax(zp_imax(d00+e0,d01+e1),zp_imax(d10+e0,d11+e1)) - 1;
    for ( i = 1 ; (1<<i ) <= d ; i++ );

    cbits = i + 2*zp_bits(p) + 1;
    if ( cbits < GMP_NUMB_BITS ) cbits = GMP_NUMB_BITS;

    mpz_init2 (v0, e0*cbits);  ks_import (v0, b0, e0, cbits);
    mpz_init2 (v1, e1*cbits);  ks_import (v1, b1, e1, cbits);
    mpz_init2 (v2, zp_imax(d00,d10)*cbits);  ks_import (v2, a00, d00, cbits);
    mpz_init2 (v3, zp_imax(d01,d11)*cbits);  ks_import (v3, a01, d01, cbits);
    mpz_init2 (v, d*cbits);
    mpz_mul (v, v0, v2);
    ks_import (v2, a10, d10, cbits);
    ks_export_addsub_xn_modp (c0, v, 0, n, cbits, p);
    mpz_mul (v, v1, v3);
    ks_export_addsub_xn_modp (c0, v, n, n, cbits, p);
    ks_import (v3, a11, d11, cbits);
    mpz_mul (v, v0, v2);
    ks_export_addsub_xn_modp (c1, v, 0, n, cbits, p);
    mpz_mul (v, v1, v3);
    ks_export_addsub_xn_modp (c1, v, n, n, cbits, p);
    mpz_clear (v); mpz_clear (v0); mpz_clear (v1); mpz_clear (v2); mpz_clear (v3);
}

// C = A*B mod x^n where A and B are 2x2 poly matrices
// A has entries aij of deg dij and B has entried bij of deg eij
// C may alias A or B but A and B must be distinct
// TODO: add support for Strassen
void zp_array_mul4_mod_xn_ks (mpz_t c00[], mpz_t c01[], mpz_t c10[], mpz_t c11[],
  mpz_t a00[], int d00, mpz_t a01[], int d01, mpz_t a10[], int d10, mpz_t a11[], int d11,
  mpz_t b00[], int e00, mpz_t b01[], int e01, mpz_t b10[], int e10, mpz_t b11[], int e11,
  int n, mpz_t p)
{
    mpz_t v, v00, v01, v10, v11, w0, w1;
    int i, d, cbits;
    
    d = zp_imax(zp_imax(zp_imax(d00+e00,d01+e10),zp_imax(d00+e01,d01+e11)),
                zp_imax(zp_imax(d10+e00,d11+e10),zp_imax(d10+e01,d11+e11)))- 1;
    for ( i = 1 ; (1<<i ) <= d ; i++ );

    cbits = i + 2*zp_bits(p) + 1;
    if ( cbits < GMP_NUMB_BITS ) cbits = GMP_NUMB_BITS;

    mpz_init2 (v00, d00*cbits);  ks_import (v00, a00, d00, cbits);
    mpz_init2 (v01, d01*cbits);  ks_import (v01, a01, d01, cbits);
    mpz_init2 (v10, d10*cbits);  ks_import (v10, a10, d10, cbits);
    mpz_init2 (v11, d11*cbits);  ks_import (v11, a11, d11, cbits);
    mpz_init2 (w0, zp_imax(e00,e01)*cbits);  ks_import (w0, b00, e00, cbits);
    mpz_init2 (w1, zp_imax(e10,e11)*cbits);  ks_import (w1, b10, e10, cbits);
    mpz_init2 (v, d*cbits);
    mpz_mul (v, v00, w0);
    ks_export_addsub_xn_modp (c00, v, 0, n, cbits, p);
    mpz_mul (v, v01, w1);
    ks_export_addsub_xn_modp (c00, v, n, n, cbits, p);
    mpz_mul (v, v10, w0);
    ks_export_addsub_xn_modp (c10, v, 0, n, cbits, p);
    mpz_mul (v, v11, w1);
    ks_export_addsub_xn_modp (c10, v, n, n, cbits, p);
    ks_import (w0, b01, e01, cbits);  ks_import (w1, b11, e11, cbits);
    mpz_mul (v, v00, w0);
    ks_export_addsub_xn_modp (c01, v, 0, n, cbits, p);
    mpz_mul (v, v01, w1);
    ks_export_addsub_xn_modp (c01, v, n, n, cbits, p);
    mpz_mul (v, v10, w0);
    ks_export_addsub_xn_modp (c11, v, 0, n, cbits, p);
    mpz_mul (v, v11, w1);
    ks_export_addsub_xn_modp (c11, v, n, n, cbits, p);
    mpz_clear (v); mpz_clear (v00); mpz_clear (v01); mpz_clear (v10); mpz_clear (v11);  mpz_clear (w0);  mpz_clear (w1);
}

// ==========================================================================

// h = ((h mod x^|m|) + sign(m) * f1 * f2) mod x^n.  Here sign(0)=1.
// be clever about squaring and mult by 0 or 1
void zp_array_mul_mod_xn (mpz_t h[], mpz_t f1[], int d1, mpz_t f2[], int d2,
                          int m, int n, mpz_t p)
{
  fermat_vec_t vec1, vec2;
  fermat_vec_struct *v;
  mp_limb_t *temp;
  size_t fft_m;
  int i, k, K, d3;
//mpz_t *h2;

  // if f1*f2 = 0, zero-fill the output and return
  if (!d1 || !d2) { for (i = (m<0?-m:m); i < n; i++) mpz_set_ui(h[i],0);  return; }

  // truncate input lengths and compute output length
  if (d1 > n) d1 = n;
  if (d2 > n) d2 = n;
  d3 = d1 + d2 - 1;

//zp_array_mul_mod_xn_ks (h, f1, d1, f2, d2, m, n, p); return;

//h2 = mpz_vec_alloc (n);
//if ( m ) zp_poly_copy (h2, h, (m < 0 ?  -1-m : m-1));
//zp_array_mul_mod_xn_ks (h2, f1, d1, f2, d2, m, n, p);

  // use naive KS if p is small or if  degree is small or if degree is much larger than the size of p
  k = zp_bits(p);
  if ( k < ZP_POLY_TFT_MIN_PBITS || d3 / (2*k) > ZP_POLY_TFT_MAX_DEG_PBIT_RATIO ) { zp_array_mul_mod_xn_ks (h, f1, d1, f2, d2, m, n, p); return; }

  if (is1(f1,d1)) { zp_array_addsub_modp (h, f2, d2, m, n, p); return; }
  if (is1(f2,d2)) { zp_array_addsub_modp (h, f1, d1, m, n, p); return; }

  // setup, import, and transform f1 into vec1  
  for (k = 0, K = 1; K < d3; k++, K <<= 1);
  fft_m = choose_m (2 * zp_bits(p) + k + 3, K);
  fermat_vec_init (vec1, K, fft_m);
  fermat_vec_import_tft (vec1, f1, d1, k, d3, 0);

  // and similarly for f2 and vec2, but be clever about squaring
  if ((f1 == f2) && (d1 == d2))
      v = vec1;
  else
    {
      fermat_vec_init (vec2, K, fft_m);
      fermat_vec_import_tft (vec2, f2, d2, k, d3, 0);
      v = vec2;
    }

  // do the multiplication into vec1
  temp = mp_limb_vec_alloc (2 * fft_m);
  fermat_vec_mul (vec1, vec1, v, d3, k, fft_m, temp);
  mp_limb_vec_free (temp, 2 * fft_m);

  // inverse transform and export the result
  fermat_vec_sqrt2exp_itft (vec1, d3, k, fft_m, 0);
  fermat_vec_export_addsub_xn_modp (h, vec1, d3, m, n, p);
  fermat_vec_clear (vec1);    if (v != vec1) fermat_vec_clear (v);

//for ( i = 0 ; i < n ; i++ ) if ( mpz_cmp(h[i],h2[i]) ) break;
//if ( i < n ){ gmp_printf ("zp_array_mul_mod_xn_ks (%d,%d,%d,%d) failed at coefficient %d with p=%Zd\n", d1, d2, m, n, i, p); gmp_printf ("%Zd != %Zd\n", h[i], h2[i]); printf("h = "); zp_poly_print (h, n-1); puts (""); printf ("h2 = "); zp_poly_print (h2, n-1); puts (""); abort(); }
//else //gmp_printf ("zp_array_mul_mod_xn_ks (%d,%d,%d,%d) verified with p=%Zd\n", d1, d2, m, n, p); 
//mpz_vec_free (h2,n);
}


// h1 = f0*f1 mod x^n, h2 = f0*f2 mod x^n (zero-padded to length n)
// we are not clever about squaring, but do handle mult by 0
// if f1 and f2 differ greatly in size, this function should not be used.
// h1 and h2 may alias any of the inputs
// This function uses 3 tfts, 2 muls, 2 itfts, and 2 mods
// *** not currently used ***
void zp_array_mul12_mod_xn (mpz_t h1[], mpz_t h2[], mpz_t f0[], int d0,
                            mpz_t f1[], int d1, mpz_t f2[], int d2, int n, mpz_t p)
{
  fermat_vec_t vec0, vec1, vec2;
  mp_limb_t *temp;
  size_t m;
  int d, k, K;

  // handle zero cases explicitly
  if (!d0) { zp_array_zero (h1, n); zp_array_zero (h2, n); return; }
  if (!d1) { zp_array_mul_mod_xn (h2, f0, d0, f2, d2, 0, n, p); zp_array_zero (h1, n); return; }

if (!d2) { zp_array_mul_mod_xn (h1, f0, d0, f1, d1, 0, n, p); zp_array_zero (h2, n); return; }
  
  // truncate input lengths, compute output length d
  if (d0 > n) d0 = n;
  if (d1 > n) d1 = n;
  if (d2 > n) d2 = n;
  d = zp_imax (d0+d1, d0+d2) - 1;

  // revert to zp_array_mul_mod_xn_ks if p is small or degree/pbits ratio is large  
  k = zp_bits(p);
  if ( k < ZP_POLY_TFT_MIN_PBITS || d / (2*k) > ZP_POLY_TFT_MAX_DEG_PBIT_RATIO )
    { zp_array_mul12_mod_xn_ks (h1, h2, f0, d0, f1, d1, f2, d2, n, p); return; }
  
  for (k = 0, K = 1; K < d; k++, K <<= 1);
  m = choose_m (2 * zp_bits(p) + k + 3, K);

  // import/transform fs into vecs (note we need 3 vecs to handle aliasing)
  fermat_vec_init (vec0, K, m);  fermat_vec_import_tft (vec0, f0, d0, k, d, 0);
  fermat_vec_init (vec1, K, m);  fermat_vec_import_tft (vec1, f1, d1, k, d, 0);
  fermat_vec_init (vec2, K, m);  fermat_vec_import_tft (vec2, f2, d2, k, d, 0);

  // multiply/itf/export vecs into hs
  temp = mp_limb_vec_alloc (2 * m);
  fermat_vec_mul (vec1, vec1, vec0, d, k, m, temp);
  fermat_vec_sqrt2exp_itft (vec1, d, k, m, 0);
  fermat_vec_export_xn_modp (h1, vec1, d, n, p);
  fermat_vec_mul (vec2, vec2, vec0, d, k, m, temp);
  fermat_vec_sqrt2exp_itft (vec2, d, k, m, 0);
  fermat_vec_export_xn_modp (h2, vec2, d, n, p);
  mp_limb_vec_free (temp, 2 * m);

  // cleanup
  fermat_vec_clear (vec0); fermat_vec_clear (vec1);  fermat_vec_clear (vec2);
}


// h = (f1 * f2 + f3 * f4) mod x^n (output zero-padded to length n)
// we are clever about squaring, and multiplication by 0 (but not 1)
// if f1*f2 and f3*f4 differ greatly in size, this function should not be used.
// This function uses 4 tfts, 2 muls, 1 itft, and 1 mod
void zp_array_mul2_mod_xn (mpz_t h[], mpz_t f1[], int d1, mpz_t f2[], int d2,
               mpz_t f3[], int d3, mpz_t f4[], int d4, int n, mpz_t p)
{
  fermat_vec_t vec1, vec2, vec3;
  fermat_vec_struct *v;
  mp_limb_t *temp;
  size_t m;
  int d, k, K, init2;
    
  // if either product is zero, revert to zp_array_mul_mod_xn
  if (!d1 || !d2) { zp_array_mul_mod_xn (h, f3, d3, f4, d4, 0, n, p); return; }
  if (!d3 || !d4) { zp_array_mul_mod_xn (h, f1, d1, f2, d2, 0, n, p); return; }

  // truncate input lengths, compute the output length d
  if (d1 > n) d1 = n;
  if (d2 > n) d2 = n;
  if (d3 > n) d3 = n;
  if (d4 > n) d4 = n;
  d = zp_imax (d1+d2, d3+d4) - 1;

  // revert to zp_array_mul_mod_xn_ks if p is small or degree/pbits ratio is large  
  k = zp_bits(p);
  if ( k < ZP_POLY_TFT_MIN_PBITS || d / (2*k) > ZP_POLY_TFT_MAX_DEG_PBIT_RATIO )
    { zp_array_mul2_mod_xn_ks (h, f1, d1, f2, d2, f3, d3, f4, d4, n, p); return; }
  
  for (k = 0, K = 1; K < d; k++, K <<= 1);
  m = choose_m (2 * zp_bits(p) + k + 3, K);

  // import and transform f1 into vec1
  fermat_vec_init (vec1, K, m);
  fermat_vec_import_tft (vec1, f1, d1, k, d, 0);

  // and similarly for f2 and vec2, but be clever about squaring
  init2 = 0;
  if ((f1 == f2) && (d1 == d2))
    v = vec1;
  else
    {
      fermat_vec_init (vec2, K, m);  init2 = 1;
      fermat_vec_import_tft (vec2, f2, d2, k, d, 0);
      v = vec2;
    }
    
  // multiply transformed f1 and f2 into vec1
  temp = mp_limb_vec_alloc (2 * m);
  fermat_vec_mul (vec1, vec1, v, d, k, m, temp);

  // setup, import and transform f3 into vec3
  fermat_vec_init (vec3, K, m);
  fermat_vec_import_tft (vec3, f3, d3, k, d, 0);

  // import and transform f4 into vec2, but be clever about squaring
  if ((f3 == f4) && (d3 == d4))
    v = vec3;
  else
    {
      if (!init2) { fermat_vec_init (vec2, K, m); init2 = 1; }
      fermat_vec_import_tft (vec2, f4, d4, k, d, 0);
      v = vec2;
    }

  // multiply transformed f3 and f4 into vec3, add into vec1, and itft
  fermat_vec_mul (vec3, vec3, v, d, k, m, temp);
  fermat_vec_addto (vec1, vec3, d, m);
  fermat_vec_sqrt2exp_itft (vec1, d, k, m, 0);

  // export to zp_poly, truncating/padding to n as required
  fermat_vec_export_xn_modp (h, vec1, d, n, p);
    
  // clean up (some of this could be done earlier)
  mp_limb_vec_free (temp, 2 * m);
  fermat_vec_clear (vec1);  fermat_vec_clear (vec3);
  if ( init2 ) fermat_vec_clear (vec2);
}

// new version of zp_array_mul2_mod_xn
// REQUIRES that (f1*f2 + f3*f4) has degree < n
// *** this function is not currently used ***
void zp_array_mul2n (mpz_t h[], mpz_t f1[], int d1, mpz_t f2[], int d2,
             mpz_t f3[], int d3, mpz_t f4[], int d4, int n, mpz_t p)
{
  fermat_vec_t vec1, vec2, vec3;
  fermat_vec_struct *v;
  mp_limb_t *temp;
  size_t m;
  int d, k, K, init2, order;

  // if either product is zero, revert to zp_array_mul_mod_xn
  if (!d1 || !d2) { zp_array_mul_mod_xn (h, f3, d3, f4, d4, 0, n, p); return; }
  if (!d3 || !d4) { zp_array_mul_mod_xn (h, f1, d1, f2, d2, 0, n, p); return; }

  // compute notional output length d, and k, K, and m
  d = zp_imax (d1+d2, d3+d4) - 1;

  // revert to zp_array_mul_mod_xn_ks if p is small or degree/pbits ratio is large  
  k = zp_bits(p);
  if ( k < ZP_POLY_TFT_MIN_PBITS || d / (2*k) > ZP_POLY_TFT_MAX_DEG_PBIT_RATIO )
    { zp_array_mul2_mod_xn_ks (h, f1, d1, f2, d2, f3, d3, f4, d4, n, p); return; }
  
  for (k = 0, K = 1; K < d || K < n; k++, K <<= 1);
  m = choose_m (2 * zp_bits(p) + 2 * k + 5, K);

  // select butterfly orders so that P(x) is defined over Z
  order = bit_reverse (n, k);

  // import and transform f1 into vec1
  fermat_vec_init (vec1, K, m);
  // (split vec_import_tft into import and tft since the former doesn't
  // have "order" parameter... hint hint...)
  fermat_vec_import (vec1, f1, d1);
  fermat_vec_tft (vec1, k, 0, d1, n, order);

  // and similarly for f2 and vec2, but be clever about squaring
  init2 = 0;
  if ((f1 == f2) && (d1 == d2))
    v = vec1;
  else
    {
      fermat_vec_init (vec2, K, m);  init2 = 1;
      // (ditto, see above)
      fermat_vec_import (vec2, f2, d2);
      fermat_vec_tft (vec2, k, 0, d2, n, order);
      v = vec2;
    }

  // multiply transformed f1 and f2 into vec1
  temp = mp_limb_vec_alloc (2 * m);
  fermat_vec_mul (vec1, vec1, v, n, k, m, temp);

  // setup, import and transform f3 into vec3
  fermat_vec_init (vec3, K, m);
  // (ditto, see above)
  fermat_vec_import (vec3, f3, d3);
  fermat_vec_tft (vec3, k, 0, d3, n, order);

  // import and transform f4 into vec2, but be clever about squaring
  if ((f3 == f4) && (d3 == d4))
    v = vec3;
  else
    {
      if (!init2) { fermat_vec_init (vec2, K, m); init2 = 1; }
      // (ditto, see above)
      fermat_vec_import (vec2, f4, d4);
      fermat_vec_tft (vec2, k, 0, d4, n, order);
      v = vec2;
    }

  // multiply transformed f3 and f4 into vec3, add into vec1, and itft
  fermat_vec_mul (vec3, vec3, v, n, k, m, temp);
  fermat_vec_addto (vec1, vec3, n, m);
  // TODO: the mul_sqrt2exp calls should be factored out of fermat_vec_mul
  // and performed HERE instead. Although this will reduce cache-friendliness
  // unless merged with addto loop above... argggh whatever.
  fermat_vec_sqrt2exp_itft (vec1, n, k, m, order);

  // export to zp_poly
  fermat_vec_exports_xn_modp (h, vec1, n, n, p);
    
  // clean up (some of this could be done earlier)
  mp_limb_vec_free (temp, 2 * m);
  fermat_vec_clear (vec1);  fermat_vec_clear (vec3);
  if ( init2 ) fermat_vec_clear (vec2);
}

// c0 = (a00 * b0 + a01 * b1) mod x^n, c1 = (a10 * b0 + a11 * b1) mod x^n
// i.e. c=a*b where a is a 2x2 matrix and b and c are 2-vectors.
// dij is the length of aij, ei is the length of bi
// This function uses 6 tfts, 4 muls, 2 itfts, 2mods
// *** gcd code uses zp_array_mul42n below, not this function ***
void zp_array_mul42_mod_xn (mpz_t c0[], mpz_t c1[],
  mpz_t a00[], int d00, mpz_t a01[], int d01, mpz_t a10[], int d10, mpz_t a11[], int d11,
  mpz_t b0[], int e0, mpz_t b1[], int e1, int n, mpz_t p)
{
  fermat_vec_t vec0, vec1, vec2, vec3;
  mp_limb_t *temp;
  size_t m;
  int d, k, K;

  // truncate input lengths, compute the output length d, and k, K, and m
  if (d00 > n ) d00 = n;
  if (d01 > n ) d01 = n;
  if (d10 > n ) d10 = n;
  if (d11 > n ) d11= n;
  if (e0 > n ) e0 = n;
  if (e1 > n ) e1 = n;
  d = zp_imax(zp_imax(d00+e0,d01+e1),zp_imax(d10+e0,d11+e1)) - 1;
  if ( d <= 0 ) { zp_array_zero (c0, n); zp_array_zero (c1, n); return; }

  // revert to zp_array_mul_mod_xn_ks if p is small or degree/pbits ratio is large  
  k = zp_bits(p);
  if ( k < ZP_POLY_TFT_MIN_PBITS || d / (2*k) > ZP_POLY_TFT_MAX_DEG_PBIT_RATIO )
    { zp_array_mul42_mod_xn_ks (c0, c1, a00, d00, a01, d01, a10, d10, a11, d11, b0, e0, b1, e1, n, p); return; }
  
  for (k = 0, K = 1; K < d; k++, K <<= 1);
  m = choose_m (2 * zp_bits(p) + k + 5, K);

  // setup, import and transform b0 and b1 into vec0 and vec1
  fermat_vec_init (vec0, K, m);
  fermat_vec_import_tft (vec0, b0, e0, k, d, 0);
  fermat_vec_init (vec1, K, m);
  fermat_vec_import_tft (vec1, b1, e1, k, d, 0);

  // import and transform a00 and a01 into vec2 and vec3
  fermat_vec_init (vec2, K, m);
  fermat_vec_import_tft (vec2, a00, d00, k, d, 0);
  fermat_vec_init (vec3, K, m);
  fermat_vec_import_tft (vec3, a01, d01, k, d, 0);
    
  // compute vec2*vec0 + vec3*vec1 and export to c0
  // without modifying vec0 or vec1
  temp = mp_limb_vec_alloc (2 * m);
  fermat_vec_mul (vec2, vec2, vec0, d, k, m, temp);
  fermat_vec_mul (vec3, vec3, vec1, d, k, m, temp);
  fermat_vec_addto (vec2, vec3, d, m);
  fermat_vec_sqrt2exp_itft (vec2, d, k, m, 0);
  fermat_vec_export_xn_modp (c0, vec2, d, n, p);

  // import and transform a10 and a11 into vec2 and vec3
  fermat_vec_import_tft (vec2, a10, d10, k, d, 0);
  fermat_vec_import_tft (vec3, a11, d11, k, d, 0);

  // compute vec2*vec0 + vec3*vec1 and export to c1
  fermat_vec_mul (vec2, vec2, vec0, d, k, m, temp);
  fermat_vec_mul (vec3, vec3, vec1, d, k, m, temp);
  fermat_vec_addto (vec2, vec3, d, m);
  fermat_vec_sqrt2exp_itft (vec2, d, k, m, 0);
  fermat_vec_export_xn_modp (c1, vec2, d, n, p);

  // clean up
  mp_limb_vec_free (temp, 2 * m);
  fermat_vec_clear (vec0);  fermat_vec_clear (vec1);
  fermat_vec_clear (vec2);  fermat_vec_clear (vec3);
}

// c0 = a00 * b0 + a01 * b1 and c1 = a10 * b0 + a11 * b1
// where the output lengths of c0 and c1 are known to be at most n
// dij is the length of aij, ei is the length of bi
// This function uses 6 tfts, 4 muls, 2 itfts, 2mods
void zp_array_mul42n (mpz_t c0[], mpz_t c1[],
  mpz_t a00[], int d00, mpz_t a01[], int d01, mpz_t a10[], int d10, mpz_t a11[], int d11,
  mpz_t b0[], int e0, mpz_t b1[], int e1, int n, mpz_t p)
{
  fermat_vec_t v0, v1, v2, v3;
  mp_limb_t *temp;
  size_t m;
  int d, k, K, o;

  // compute the nominal output length d (may be larger than n), and k, K, and m
  d = zp_imax(zp_imax(d00+e0,d01+e1),zp_imax(d10+e0,d11+e1)) - 1;
  if ( d <= 0 ) { zp_array_zero (c0, n); zp_array_zero (c1, n); return; }

  // revert to zp_array_mul_mod_xn_ks if p is small or degree/pbits ratio is large  
  k = zp_bits(p);
  if ( k < ZP_POLY_TFT_MIN_PBITS || d / (2*k) > ZP_POLY_TFT_MAX_DEG_PBIT_RATIO )
    { zp_array_mul42_mod_xn_ks (c0, c1, a00, d00, a01, d01, a10, d10, a11, d11, b0, e0, b1, e1, n, p); return; }

  for (k = 0, K = 1; K < d || K < n; k++, K <<= 1);
  m = choose_m (2 * zp_bits(p) + 2 * k + 5, K);
  
  // select butterfly orders so that P(x) is defined over Z
  o = bit_reverse (n, k);

  // init, import and transform b0 and b1 into v0 and v1
  fermat_vec_init (v0, K, m);  
  fermat_vec_import_tft (v0, b0, e0, k, n, o);
  fermat_vec_init (v1, K, m);
  fermat_vec_import_tft (v1, b1, e1, k, n, o);

  // init, import and transform a00 and a01 into v2 and v3
  fermat_vec_init (v2, K, m);
  fermat_vec_import_tft (v2, a00, d00, k, n, o);
  fermat_vec_init (v3, K, m);
  fermat_vec_import_tft (v3, a01, d01, k, n, o);
    
  // compute v2*v0 + v3*v1 and export to c0
  // without modifying v0 or v1
  temp = mp_limb_vec_alloc (2 * m);
  fermat_vec_mul (v2, v2, v0, n, k, m, temp);
  fermat_vec_mul (v3, v3, v1, n, k, m, temp);
  fermat_vec_addto (v2, v3, n, m);
  fermat_vec_sqrt2exp_itft (v2, n, k, m, o);
  fermat_vec_exports_xn_modp (c0, v2, n, n, p);

  // import and transform a10 and a11 into v2 and v3
  fermat_vec_import_tft (v2, a10, d10, k, n, o);
  fermat_vec_import_tft (v3, a11, d11, k, n, o);

  // compute v2*v0 + v3*v1 and export to c1
  fermat_vec_mul (v2, v2, v0, n, k, m, temp);
  fermat_vec_mul (v3, v3, v1, n, k, m, temp);
  fermat_vec_addto (v2, v3, n, m);
  fermat_vec_sqrt2exp_itft (v2, n, k, m, o);
  fermat_vec_exports_xn_modp (c1, v2, n, n, p);

  // clean up
  mp_limb_vec_free (temp, 2 * m);
  fermat_vec_clear (v0);  fermat_vec_clear (v1);
  fermat_vec_clear (v2);  fermat_vec_clear (v3);
}

// C = A*B mod x^n where A and B are 2x2 poly matrices
// A has entries aij of deg dij and B has entried bij of deg eij
// C may alias A or B but A and B must be distinct
// This function uses 8 tfts, 8 muls, 4 itfts and 4 mods
void zp_array_mul4_mod_xn (mpz_t c00[], mpz_t c01[], mpz_t c10[], mpz_t c11[],
  mpz_t a00[], int d00, mpz_t a01[], int d01, mpz_t a10[], int d10, mpz_t a11[], int d11,
  mpz_t b00[], int e00, mpz_t b01[], int e01, mpz_t b10[], int e10, mpz_t b11[], int e11,
  int n, mpz_t p)
{
  fermat_vec_t v00, v01, v10, v11, t0, t1, t2, t3;
  mp_limb_t *temp;
  size_t m;
  int d, k, K;

  // truncate input lengths and compute the output size d
  if (d00 > n ) d00 = n;
  if (d01 > n ) d01 = n;
  if (d10 > n ) d10 = n;
  if (d11 > n ) d11= n;
  if (e00 > n ) e00 = n;
  if (e01 > n ) d01 = n;
  if (e10 > n ) e10 = n;
  if (e11 > n ) e11= n;
  d = zp_imax(zp_imax(zp_imax(d00+e00,d01+e10),zp_imax(d00+e01,d01+e11)),
          zp_imax(zp_imax(d10+e00,d11+e10),zp_imax(d10+e01,d11+e11)))- 1;
  if ( d <= 0 )
    { zp_array_zero (c00, n); zp_array_zero (c01, n); zp_array_zero(c10, n); zp_array_zero (c11, n); return; }
    
  // revert to zp_array_mul_mod_xn_ks if p is small or degree/pbits ratio is large  
  k = zp_bits(p);
  if ( k < ZP_POLY_TFT_MIN_PBITS || d / (2*k) > ZP_POLY_TFT_MAX_DEG_PBIT_RATIO )
    { zp_array_mul4_mod_xn_ks (c00, c01, c10, c11, a00, d00, a01, d01, a10, d10, a11, d11,
                           b00, e00, b01, e01, b10, e10, b11, e11, n, p); return; }


  // use Strassen if p is big enough
  if ( ZP_POLY_STRASSEN_MIN_PBITS >= 0 && zp_bits(p) >= ZP_POLY_STRASSEN_MIN_PBITS )
    {
      zp_array_mul4_mod_xn_strassen (c00, c01, c10, c11, a00, d00, a01, d01, a10, d10, a11, d11,
                     b00, e00, b01, e01, b10, e10, b11, e11, n, p);
      return;
    }

  for (k = 0, K = 1; K < d; k++, K <<= 1);
  m = choose_m (2 * zp_bits(p) + k + 3, K);

  // import and transform aij into vij
  fermat_vec_init (v00, K, m);  fermat_vec_import_tft (v00, a00, d00, k, d, 0);
  fermat_vec_init (v01, K, m);  fermat_vec_import_tft (v01, a01, d01, k, d, 0);
  fermat_vec_init (v10, K, m);  fermat_vec_import_tft (v10, a10, d10, k, d, 0);
  fermat_vec_init (v11, K, m);  fermat_vec_import_tft (v11, a11, d11, k, d, 0);

  // import and transform b00 and b10 into t0 and t1, init t2 and t3
  fermat_vec_init (t0, K, m);  fermat_vec_import_tft (t0, b00, e00, k, d, 0);
  fermat_vec_init (t1, K, m);  fermat_vec_import_tft (t1, b10, e10, k, d, 0);
  fermat_vec_init (t2, K, m);
  fermat_vec_init (t3, K, m);

  // compute c00 and c10 using t2 and t3 as temporaries
  temp = mp_limb_vec_alloc (2 * m);
  fermat_vec_mul (t2, v00, t0, d, k, m, temp);
  fermat_vec_mul (t3, v01, t1, d, k, m, temp);
  fermat_vec_addto (t2, t3, d, m);
  fermat_vec_sqrt2exp_itft (t2, d, k, m, 0);
  fermat_vec_export_xn_modp (c00, t2, d, n, p);
  fermat_vec_mul (t2, v10, t0, d, k, m, temp);
  fermat_vec_mul (t3, v11, t1, d, k, m, temp);
  fermat_vec_addto (t2, t3, d, m);
  fermat_vec_sqrt2exp_itft (t2, d, k, m, 0);
  fermat_vec_export_xn_modp (c10, t2, d, n, p);

  // import and transform b01 and b11 into t0 and t1
  fermat_vec_import_tft (t0, b01, e01, k, d, 0);
  fermat_vec_import_tft (t1, b11, e11, k, d, 0);

  // compute c01 and c11 using t2 and t3 as temporaries
  fermat_vec_mul (t2, v00, t0, d, k, m, temp);
  fermat_vec_mul (t3, v01, t1, d, k, m, temp);
  fermat_vec_addto (t2, t3, d, m);
  fermat_vec_sqrt2exp_itft (t2, d, k, m, 0);
  fermat_vec_export_xn_modp (c01, t2, d, n, p);
  fermat_vec_mul (t2, v10, t0, d, k, m, temp);
  fermat_vec_mul (t3, v11, t1, d, k, m, temp);
  fermat_vec_addto (t2, t3, d, m);
  fermat_vec_sqrt2exp_itft (t2, d, k, m, 0);
  fermat_vec_export_xn_modp (c11, t2, d, n, p);

  // clean up
  mp_limb_vec_free (temp, 2 * m);
  fermat_vec_clear (v00); fermat_vec_clear (v01); fermat_vec_clear (v10); fermat_vec_clear (v11);
  fermat_vec_clear (t0); fermat_vec_clear (t1); fermat_vec_clear (t2); fermat_vec_clear (t3);
}

// same interface as zp_array_mul4_mod_xn
// uses Strassen to trade 1 mult for 11 adds
void zp_array_mul4_mod_xn_strassen (mpz_t c00[], mpz_t c01[], mpz_t c10[], mpz_t c11[],
  mpz_t a00[], int d00, mpz_t a01[], int d01, mpz_t a10[], int d10, mpz_t a11[], int d11,
  mpz_t b00[], int e00, mpz_t b01[], int e01, mpz_t b10[], int e10, mpz_t b11[], int e11,
  int n, mpz_t p)
{
  fermat_vec_t t1, t2, t3, t4, t5, t6, t7, t8, t9;
  mp_limb_t *temp;
  size_t m;
  int d, k, K;

  // truncate input lengths and compute the output size d
  if (d00 > n ) d00 = n;
  if (d01 > n ) d01 = n;
  if (d10 > n ) d10 = n;
  if (d11 > n ) d11= n;
  if (e00 > n ) e00 = n;
  if (e01 > n ) d01 = n;
  if (e10 > n ) e10 = n;
  if (e11 > n ) e11= n;
  d = zp_imax(zp_imax(zp_imax(d00+e00,d01+e10),zp_imax(d00+e01,d01+e11)),
          zp_imax(zp_imax(d10+e00,d11+e10),zp_imax(d10+e01,d11+e11)))- 1;
  if ( d <= 0 )
    { zp_array_zero (c00, n); zp_array_zero (c01, n); zp_array_zero(c10, n); zp_array_zero (c11, n); return; }

  for (k = 0, K = 1; K < d; k++, K <<= 1);
  m = choose_m (2 * zp_bits(p) + k + 3, K);

  // initialize
  fermat_vec_init (t1, K, m); fermat_vec_init (t2, K, m); fermat_vec_init (t3, K, m);
  fermat_vec_init (t4, K, m); fermat_vec_init (t5, K, m); fermat_vec_init (t6, K, m);
  fermat_vec_init (t7, K, m); fermat_vec_init (t8, K, m); fermat_vec_init (t9, K, m);
  temp = mp_limb_vec_alloc (2 * m);
    
  // the following is a modified version of the code in Alg 12.1 in [GG] (uses fewer temps)
  // labels in parens refer to the intermediate symbols used in [GG] (see p. 328)
  // note that we index here from 0, so a01 is A_12 in [GG].
  fermat_vec_import_tft (t5, a00, d00, k, d, 0);    // t5 = a00
  fermat_vec_import_tft (t6, b00, e00, k, d, 0);    // t6 = b00
  fermat_vec_mul (t1, t5, t6, d, k, m, temp);       // t1 = t5 * t6 (P1)
  fermat_vec_import_tft (t2, a01, d01, k, d, 0);    // t2 = a01
  fermat_vec_import_tft (t3, b10, e10, k, d, 0);    // t3 = b10
  fermat_vec_mul (t4, t2, t3, d, k, m, temp);       // t4 = t2 * t3 (P2)
  fermat_vec_addto (t4, t1, d, m);          // t4 = t4 + t1 (U1 = P1 + P2)
  fermat_vec_sqrt2exp_itft (t4, d, k, m, 0);
  fermat_vec_export_xn_modp (c00, t4, d, n, p);     // c00 = t4 (overwriting a00 or b00 is ok)
  fermat_vec_import_tft (t4, a10, d10, k, d, 0);    // t4 = a10
  fermat_vec_import_tft (t7, a11, d11, k, d, 0);    // t7 = a11
  fermat_vec_addto (t4, t7, d, m);          // t4 = t4 + t7 (S1)
  fermat_vec_sub (t5, t5, t4, d, m);            // t5 = t5 - t4 (-S2)
  fermat_vec_import_tft (t8, b01, e01, k, d, 0);    // t8 = b01
  fermat_vec_sub (t8, t8, t6, d, m);            // t8 = t8 - t6 (T1)
  fermat_vec_mul (t4, t4, t8, d, k, m, temp);       // t4 = t4 * t8 (P5)
  fermat_vec_addto (t1, t4, d, m);          // t1 = t1 + t4 (P1 + P5)
  fermat_vec_addto (t2, t5, d, m);          // t2 = t2 + t5 (S4)
  fermat_vec_import_tft (t9, b11, e11, k, d, 0);    // t9 = b11
  fermat_vec_mul (t2, t2, t9, d, k, m, temp);       // t2 = t2 * t9 (P3)
  fermat_vec_sub (t9, t9, t8, d, m);            // t9 = t9 - t8 (T2)
  fermat_vec_add (t8, t7, t5, d, m);            // t8 = t7 + t5 (S3)
  fermat_vec_mul (t5, t5, t9, d, k, m, temp);       // t5 = t5 * t9 (-P6)
  fermat_vec_sub (t1, t1, t5, d, m);            // t1 = t1 - t5 (U4 = P1 + P5 + P6) 
  fermat_vec_addto (t2, t1, d, m);          // t2 = t2 + t1 (U5 = P1 + P3 + P5 + P6)
  fermat_vec_sqrt2exp_itft (t2, d, k, m, 0);
  fermat_vec_exports_xn_modp (c01, t2, d, n, p);    // c01 = t2 (overwriting is ok)
  fermat_vec_sub (t2, t9, t6, d, m);            // t2 = t9 - t6 (T3)
  fermat_vec_sub (t9, t9, t3, d, m);            // t9 = t9 - t3 (T4)
  fermat_vec_mul (t9, t9, t7, d, k, m, temp);       // t9 = t9 * t7 (P4)
  fermat_vec_addto (t9, t4, d, m);          // t9 = t9 + t4 (P4 + P5)
  fermat_vec_mul (t2, t2, t8, d, k, m, temp);       // t2 = t2 * t8 (P7)
  fermat_vec_addto (t2, t1, d, m);          // t2 = t2 + t1 (U7 = P1 + P5 + P6 + P7)
  fermat_vec_sub (t1, t2, t9, d, m);            // t1 = t2 - t3 (U6 = P1 - P4 + P6 + P7)
  fermat_vec_sqrt2exp_itft (t1, d, k, m, 0);
  fermat_vec_exports_xn_modp (c10, t1, d, n, p);    // c10 = t1 (overwriting is ok)
  fermat_vec_sqrt2exp_itft (t2, d, k, m, 0);
  fermat_vec_exports_xn_modp (c11, t2, d, n, p);    // c11 = t2 (overwriting is ok)
  
  // total cost above is 8 import/tft, 7 mul, 15 add, 4 itft/export/modp
    
  // clean up
  mp_limb_vec_free (temp, 2 * m);
  fermat_vec_clear (t1); fermat_vec_clear (t2); fermat_vec_clear (t3);
  fermat_vec_clear (t4); fermat_vec_clear (t5); fermat_vec_clear (t6);
  fermat_vec_clear (t7); fermat_vec_clear (t8); fermat_vec_clear (t9);
}


// ==========================================================================
// "stupid" implementation of zp_poly_mod_t, all PRIVATE


// input is g of degree d with g[d] = 1.
// computes h of degree d-1 with h[d-1] = 1 such that g*h = x^(2d-1) + r
// where deg(r) < d.
void recip_stupid (mpz_t h[], mpz_t g[], int d, mpz_t p)
{
  int i, j;

  assert (d >= 1);

  mpz_set_ui (h[d - 1], 1);
  for (i = d - 2; i >= 0; i--)
    mpz_set_ui (h[i], 0);

  for (i = d - 2; i >= 0; i--)
    {
      for (j = 0; j <= i; j++)
    mpz_submul (h[j], h[i + 1], g[j - i + d - 1]);
      mpz_mod (h[i], h[i], p);
    }
}


// assumes algo, d, p, g are already set
void zp_poly_mod_init_stupid (zp_poly_mod_t mod)
{
  mod->ginv = mpz_vec_alloc (mod->d);
  recip_stupid (mod->ginv, mod->g, mod->d, mod->p);
}


void zp_poly_mod_clear_stupid (zp_poly_mod_t mod)
{
  mpz_vec_free (mod->ginv, mod->d);
}


void zp_poly_mod_import_stupid (mpz_t h[], mpz_t f[], zp_poly_mod_t mod)
{
  int i;

  for (i = 0; i < mod->d; i++)
    mpz_set (h[i], f[i]);
}


void zp_poly_mod_import_one_stupid (mpz_t h[], zp_poly_mod_t mod)
{
  int i;

  mpz_set_ui (h[0], 1);
  for (i = 1; i < mod->d; i++)
    mpz_set_ui (h[i], 0);
}


void zp_poly_mod_export_stupid (mpz_t h[], mpz_t f[], zp_poly_mod_t mod)
{
  int i;

  for (i = 0; i < mod->d; i++)
    mpz_set (h[i], f[i]);
}


void zp_poly_mod_mul_stupid (mpz_t h[], mpz_t f1[], mpz_t f2[], int e,
                 zp_poly_mod_t mod)
{
  mpz_t* v, * q;
  int i, j;

  v = mpz_vec_alloc (2 * mod->d);
  q = mpz_vec_alloc (mod->d);

  // multiply out
  for (i = 0; i < mod->d; i++)
    for (j = 0; j < mod->d; j++)
      mpz_addmul (v[i+j+e], f1[i], f2[j]);
  for (i = 0; i < 2*mod->d; i++)
    mpz_mod (v[i], v[i], mod->p);

  // multiply top half of v by ginv to get quotient q
  for (i = 0; i < mod->d; i++)
    for (j = mod->d - i - 1; j < mod->d; j++)
      mpz_addmul (q[i + j - mod->d + 1], v[i + mod->d], mod->ginv[j]);
  for (i = 0; i < mod->d; i++)
    mpz_mod (q[i], q[i], mod->p);

  // multiply quotient by g and subtract from bottom half of v
  for (i = 0; i < mod->d; i++)
    for (j = 0; j <= mod->d - i; j++)
      mpz_submul (v[i+j], q[i], mod->g[j]);

  for (i = 0; i < mod->d; i++)
    mpz_mod (h[i], v[i], mod->p);

  mpz_vec_free (v, 2 * mod->d);
  mpz_vec_free (q, mod->d);
}



// ==========================================================================
// "KRONECKER" implementation of zp_poly_mod_t, all PRIVATE

// input is g of degree d with g[d] = 1.
// computes h of degree d-1 with h[d-1] = 1 such that g*h = x^(2d-1) + r
// where deg(r) < d.
void recip_kronecker (mpz_t h[], mpz_t g[], int d, mpz_t p)
{
  zp_poly_div_xn (h, g, d, 2*d-1, p);
}


// assumes algo, d, p, g are already set
void zp_poly_mod_init_kronecker (zp_poly_mod_t mod)
{
  int i;

  mod->ginv = mpz_vec_alloc (mod->d);
  recip_kronecker (mod->ginv, mod->g, mod->d, mod->p);

  for ( i = 1 ; (1<<i ) <= mod->d ; i++ );
  mod->cbits = i + 2*mpz_sizeinbase(mod->p,2) + 1;
  if ( mod->cbits < GMP_NUMB_BITS ) mod->cbits = GMP_NUMB_BITS;
//while ( mod->cbits&3 ) mod->cbits++;
  mpz_init2 (mod->g_ks, (mod->d+1)*mod->cbits);
  ks_import (mod->g_ks, mod->g, mod->d, mod->cbits);    // we don't need the top coeff of g
  mpz_init2 (mod->ginv_ks, mod->d*mod->cbits);
  ks_import (mod->ginv_ks, mod->ginv, mod->d, mod->cbits);
  mpz_init2 (mod->N1, mod->d*mod->cbits);
  mpz_init2 (mod->N2, mod->d*mod->cbits);
  mpz_init2 (mod->N3, 2*mod->d*mod->cbits);
  mpz_init2 (mod->n1, mod->cbits);
  mpz_init2 (mod->n2, mod->cbits);
}


void zp_poly_mod_clear_kronecker (zp_poly_mod_t mod)
{
  mpz_vec_free (mod->ginv, mod->d);
  mpz_clear (mod->g_ks);  mpz_clear (mod->ginv_ks);
  mpz_clear (mod->N1); mpz_clear (mod->N2); mpz_clear (mod->N3); 
  mpz_clear (mod->n1); mpz_clear (mod->n2);
}


void zp_poly_mod_import_kronecker (mpz_t h[], mpz_t f[], zp_poly_mod_t mod)
{
  int i;

  for (i = 0; i < mod->d; i++)
    mpz_set (h[i], f[i]);
}


void zp_poly_mod_import_one_kronecker (mpz_t h[], zp_poly_mod_t mod)
{
  int i;

  mpz_set_ui (h[0], 1);
  for (i = 1; i < mod->d; i++)
    mpz_set_ui (h[i], 0);
}


void zp_poly_mod_export_kronecker (mpz_t h[], mpz_t f[], zp_poly_mod_t mod)
{
  int i;

  for (i = 0; i < mod->d; i++)
    mpz_set (h[i], f[i]);
}


void zp_poly_mod_mul_kronecker (mpz_t h[], mpz_t f1[], mpz_t f2[], int e,
                 zp_poly_mod_t mod)
{
  int i, off;

  assert ( e == 0 || e == 1 );

  // multiply out
  ks_import (mod->N1, f1, mod->d, mod->cbits);
  ks_import (mod->N2, f2, mod->d, mod->cbits);
  mpz_mul (mod->N3, mod->N1, mod->N2);
  ks_reduce (mod->N3, 2*mod->d-1, mod->cbits, mod->p);
  mpz_fdiv_q_2exp (mod->N1, mod->N3, (mod->d-e) * mod->cbits);
  mpz_fdiv_r_2exp (mod->N2, mod->N3, (mod->d-e) * mod->cbits);
  if ( e ) mpz_mul_2exp (mod->N2, mod->N2, mod->cbits);

  // multiply top half of v by ginv to get quotient q
  mpz_mul (mod->N3, mod->N1, mod->ginv_ks);
  mpz_fdiv_q_2exp (mod->N1, mod->N3, (mod->d-1)*mod->cbits);
  ks_reduce (mod->N1, mod->d, mod->cbits, mod->p);

  // multiply quotient by g and subtract from bottom half of v
  mpz_mul (mod->N3, mod->N1, mod->g_ks);
  for ( i = off = 0 ; i < mod->d ; i++, off += mod->cbits )
    {
      ks_unpack (mod->n1, mod->N3, off, mod->cbits);  mpz_mod (mod->n1, mod->n1, mod->p);
      ks_unpack (mod->n2, mod->N2, off, mod->cbits);  if ( mpz_cmp (mod->n1, mod->n2) > 0 ) mpz_add (mod->n2, mod->n2, mod->p);
      mpz_sub (h[i], mod->n2, mod->n1);
    }
}


// ==========================================================================
// "barrett" implementation of zp_poly_mod_t, all PRIVATE


// same interface as recip_stupid()
void recip_barrett (mpz_t h[], mpz_t g[], int d, mpz_t p)
{
    zp_poly_div_xn (h, g, d, 2*d-1, p);
}


// assumes algo, d, p, g are already set
void zp_poly_mod_init_barrett (zp_poly_mod_t mod)
{
  int k, K, i;
  size_t bits, m, s;
  mpz_t* ginv;

  // choose K = 2^k >= 2*d
  k = 0;
  K = 1;
  while (K < 2 * mod->d)
    k++, K <<= 1;
  mod->k = k;
  mod->K = K;

  // choose coefficient size to accommodate products
  // TODO: think more carefully about the 2*k term here...
  bits = 2 * mpz_sizeinbase (mod->p, 2) + 2 * k + 3;
  m = choose_m (bits, K);
  mod->m = m;

  // select butterfly orders so that P(x) is defined over Z
  mod->order = bit_reverse (mod->d, k);

  // precompute transform of g, including scale factor
  mod->g_tft = (fermat_vec_struct*) zp_malloc (sizeof (fermat_vec_struct));
  fermat_vec_init (mod->g_tft, K, m);
  fermat_vec_import (mod->g_tft, mod->g, mod->d + 1);
  fermat_vec_tft (mod->g_tft, k, 0, mod->d + 1, mod->d, mod->order);

  s = 4 * GMP_NUMB_BITS * mod->m - 2 * mod->k;
  for (i = 0; i < mod->d; i++)
    fermat_mul_sqrt2exp (&mod->g_tft->x[i], &mod->g_tft->t, s, m);

  // precompute transform of ginv, including scale factor
  ginv = mpz_vec_alloc (mod->d);
  recip_barrett (ginv, mod->g, mod->d, mod->p);

  mod->ginv_tft = (fermat_vec_struct*) zp_malloc (sizeof (fermat_vec_struct));
  fermat_vec_init (mod->ginv_tft, K, m);
  fermat_vec_import (mod->ginv_tft, ginv, mod->d);
  fermat_vec_tft (mod->ginv_tft, k, 0, mod->d, 2 * mod->d - 1, mod->order);

  s = 4 * GMP_NUMB_BITS * mod->m - 2 * mod->k;
  for (i = 0; i < 2 * mod->d - 1; i++)
    fermat_mul_sqrt2exp (&mod->ginv_tft->x[i], &mod->ginv_tft->t, s, m);

  mpz_vec_free (ginv, mod->d);
}


void zp_poly_mod_clear_barrett (zp_poly_mod_t mod)
{
  fermat_vec_clear (mod->ginv_tft);
  fermat_vec_clear (mod->g_tft);
  zp_free (mod->ginv_tft, sizeof (fermat_vec_struct));
  zp_free (mod->g_tft, sizeof (fermat_vec_struct));
}


void zp_poly_mod_import_barrett (mpz_t h[], mpz_t f[], zp_poly_mod_t mod)
{
  int i;

  for (i = 0; i < mod->d; i++)
    mpz_set (h[i], f[i]);
}


void zp_poly_mod_import_one_barrett (mpz_t h[], zp_poly_mod_t mod)
{
  int i;

  mpz_set_ui (h[0], 1);
  for (i = 1; i < mod->d; i++)
    mpz_set_ui (h[i], 0);
}


void zp_poly_mod_export_barrett (mpz_t h[], mpz_t f[], zp_poly_mod_t mod)
{
  int i;

  for (i = 0; i < mod->d; i++)
    mpz_set (h[i], f[i]);
}


void zp_poly_mod_mul_barrett (mpz_t h[], mpz_t f1[], mpz_t f2[], int e,
                  zp_poly_mod_t mod)
{
  mpz_t* v, * q;
  int i, sqr;
  ptrdiff_t s;
  size_t w;
  fermat_vec_t vec1, vec2;
  mp_limb_t* temp;

  sqr = (f1 == f2);

  v = mpz_vec_alloc (2 * mod->d);
  q = mpz_vec_alloc (mod->d);
  temp = mp_limb_vec_alloc (2 * mod->m);
  fermat_vec_init (vec1, mod->K, mod->m);
  fermat_vec_init (vec2, mod->K, mod->m);

  // ================   v = f1 * f2

  // convert polys to fermat polys and transform
  fermat_vec_import (vec1, f1, mod->d);
  fermat_vec_tft (vec1, mod->k, 0, mod->d, 2 * mod->d, 0);

  if (!sqr)
    {
      fermat_vec_import (vec2, f2, mod->d);
      fermat_vec_tft (vec2, mod->k, 0, mod->d, 2 * mod->d, 0);
    }

  // K-th root of unity
  w = (4 * GMP_NUMB_BITS * mod->m) >> mod->k;

  // multiply fourier coefficients, scale, and twist if e == 1
  for (i = 0; i < 2 * mod->d; i++)
    {
      fermat_mul (vec1->x[i], vec1->x[i],
          (sqr ? vec1 : vec2)->x[i], temp, mod->m, vec1->best_k);

      // TODO: we *should* be able to push the twist directly into the
      // twist parameter of the TFT, need to think about this some more...
      s = e ? (bit_reverse (i, mod->k) * w) : 0;
      s -= 2 * mod->k;
      if (s < 0)
    s += 4 * GMP_NUMB_BITS * mod->m;
      fermat_mul_sqrt2exp (&vec1->x[i], &vec1->t, s, mod->m);
    }

  // inverse transform
  fermat_vec_itft (vec1, mod->k, 0, 2 * mod->d, 2 * mod->d, 0, 0);

  // convert back to zp_poly
  fermat_vec_export_modp (v, vec1, 0, 2 * mod->d, 0, mod->p);

  // ================   q = ((v // x^d) * ginv) // x^(d-1)

  // convert (v // x^d) to fermat poly and transform
  // (also retain transform of (v // x^d) mod P(x) for later use)
  fermat_vec_import (vec1, v + mod->d, mod->d);
  fermat_vec_tft (vec1, mod->k, 0, mod->d, 2 * mod->d - 1, mod->order);

  // multiply fourier coefficients (including scale factor)
  for (i = 0; i < 2 * mod->d - 1; i++)
    fermat_mul (vec2->x[i], vec1->x[i], mod->ginv_tft->x[i], temp,
        mod->m, vec2->best_k);

  // inverse transform
  fermat_vec_itft (vec2, mod->k, 0, 2 * mod->d - 1, 2 * mod->d - 1, 0,
           mod->order);

  // convert back to zp_poly
  fermat_vec_export_modp (q, vec2, mod->d - 1, mod->d, 0, mod->p);

  // ================   ((q * g) - x^d * (v // x^d))  mod P(x)

  // convert q to fermat poly and transform
  fermat_vec_import (vec2, q, mod->d);
  fermat_vec_tft (vec2, mod->k, 0, mod->d, mod->d, mod->order);

  // multiply fourier coefficients and scale
  for (i = 0; i < mod->d; i++)
    {
      fermat_mul (vec2->x[i], vec2->x[i], mod->g_tft->x[i], temp,
          mod->m, vec2->best_k);

      s = mod->d * w * (bit_reverse (i, mod->k) ^ mod->order);
      s += 4 * GMP_NUMB_BITS * mod->m - 2 * mod->k;
      s %= 4 * GMP_NUMB_BITS * mod->m;
      fermat_mul_sqrt2exp (&vec1->x[i], &vec1->t, s, mod->m);
      fermat_sub (vec2->x[i], vec2->x[i], vec1->x[i], mod->m);
    }

  // inverse transform
  fermat_vec_itft (vec2, mod->k, 0, mod->d, mod->d, 0, mod->order);

  // convert back to zp_poly and subtract from v
  fermat_vec_export_modp (q, vec2, 0, mod->d, 1, mod->p);

  for (i = 0; i < mod->d; i++)
    {
      mpz_sub (v[i], v[i], q[i]);
      mpz_mod (h[i], v[i], mod->p);
    }

  // ================

  fermat_vec_clear (vec1);
  fermat_vec_clear (vec2);
  mpz_vec_free (v, 2 * mod->d);
  mpz_vec_free (q, mod->d);
  mp_limb_vec_free (temp, 2 * mod->m);
}



// ==========================================================================
// "montgomery" implementation of zp_poly_mod_t, all PRIVATE


/*
  This algorithm is a variant of Montgomery multiplication, similar to the
  McLaughlin/Mihailescu variants, adapted for the TFT.

  Let P(x) be the polynomial defined by the roots of the TFT of length d,
  where the butterfly ordering is selected so that P(x) is in Z[x].

  Let Q(x) = P(x/2) * 2^d (which is monic).

  We require that (P, Q) = (g, Q) = 1.
  If this is not satisfied, we fall back on Barrett multiplication.

  Let R = 1/g mod Q, and S = 1/Q mod P.

  Now given input polynomials f1, f2 of degree < d, and 0 <= e <= 1,
  we compute

  (1) t = (x^e * f1 * f2 * R)  mod Q,
  (2) u = (x^e * f1 * f2 - g * t) * S  mod P.

  Step (2) is done by lifting the problem to Z[x], descending to Z/FZ[x]
  where F = fermat modulus, and using the TFT/ITFT to compute everything
  modulo P.

  Step (1) is similar, after first making the substitution x -> 2x
  (i.e. to convert arithmetic mod Q to arithmetic mod P).

  Claim: u = x^e * f1 * f2 / Q mod g.

  Proof: from (1) we have g * t = x^e * f1 * f2 mod Q, so

     g * t = x^e * f1 * f2 - v * Q   for some v, len(v) <= d,      (*)

  and therefore

     v = x^e * f1 * f2 / Q  mod g.

  But reading (*) mod P we get

     x^e * f1 * f2 - g * t = v * Q   mod P,

  so

     v = (x^e * f1 * f2 - g * t) * S   mod P,

  i.e. v == u.
*/


void zp_poly_mod_clear_montgomery (zp_poly_mod_t mod)
{
  fermat_vec_clear (mod->S_tft);
  fermat_vec_clear (mod->R_tft);
  fermat_vec_clear (mod->g_tft);
  zp_free (mod->S_tft, sizeof (fermat_vec_struct));
  zp_free (mod->R_tft, sizeof (fermat_vec_struct));
  zp_free (mod->g_tft, sizeof (fermat_vec_struct));
  mpz_vec_free (mod->Q_mod_g, mod->d);
}


// assumes algo, d, p, g are already set
void zp_poly_mod_init_montgomery (zp_poly_mod_t mod)
{
  int k, K, i;
  size_t bits, m, s;
  mpz_t* P, * Q, * R, * S, * T;
  mpz_t M;

  // choose K = 2^k >= d
  k = 0;
  K = 1;
  while (K < mod->d + 1)
    k++, K <<= 1;
  mod->k = k;
  mod->K = K;

  // choose coefficient size to accommodate products
  // TODO: think more carefully about the 4*k term here...
  bits = 3 * mpz_sizeinbase (mod->p, 2) + 4 * k + 10;
  m = choose_m (bits, K);
  mod->m = m;

  // select butterfly orders so that P(x) is defined over Z
  mod->order = bit_reverse (mod->d, k);

  mpz_init (M);
  P = mpz_vec_alloc (mod->d + 1);
  Q = mpz_vec_alloc (mod->d + 1);
  R = mpz_vec_alloc (mod->d);
  S = mpz_vec_alloc (mod->d);
  T = mpz_vec_alloc (mod->d);

  mod->Q_mod_g = mpz_vec_alloc (mod->d);
  mod->R_tft = (fermat_vec_struct*) zp_malloc (sizeof (fermat_vec_struct));
  mod->S_tft = (fermat_vec_struct*) zp_malloc (sizeof (fermat_vec_struct));
  mod->g_tft = (fermat_vec_struct*) zp_malloc (sizeof (fermat_vec_struct));
  fermat_vec_init (mod->R_tft, K, m);
  fermat_vec_init (mod->S_tft, K, m);
  fermat_vec_init (mod->g_tft, K, m);

  // compute P(x) = \prod (x^(2^j) + 1) where d has bit #j set
  mpz_set_ui (P[0], 1);
  for (i = 1; i <= mod->d; i++)
    if (!(i & ~mod->d))
      mpz_set_ui (P[i], 1);

  // compute Q(x) = P(x/2) * 2^d
  for (i = 0; i <= mod->d; i++)
    {
      mpz_mul_2exp (Q[i], P[i], mod->d - i);
      mpz_mod (Q[i], Q[i], mod->p);
    }

  // compute Q mod g
  for (i = 0; i < mod->d; i++)
    {
      mpz_sub (mod->Q_mod_g[i], Q[i], mod->g[i]);
      mpz_mod (mod->Q_mod_g[i], mod->Q_mod_g[i], mod->p);
    }

  // compute R = 1/g mod Q
  // (first do T = g mod Q)
  for (i = 0; i < mod->d; i++)
    {
      mpz_sub (T[i], mod->g[i], Q[i]);
      mpz_mod (T[i], T[i], mod->p);
    }
  if (zp_poly_inv_mod (R, T, Q, mod->d, mod->p) == 0)
    goto init_failed;

  // compute S = 1/Q mod P
  // (first do T = Q mod P)
  for (i = 0; i < mod->d; i++)
    {
      mpz_sub (T[i], Q[i], P[i]);
      mpz_mod (T[i], T[i], mod->p);
    }
  if (zp_poly_inv_mod (S, T, P, mod->d, mod->p) == 0)
    goto init_failed;

  // compute transform of R
  // make substitution x -> 2x, and scale by 2^(-d+1) mod p
  mpz_set_ui (M, 1);
  mpz_mul_2exp (M, M, mod->d - 1);
  mpz_invert (M, M, mod->p);    // 2^(-d+1) mod p
  for (i = 0; i < mod->d; i++)
    {
      mpz_mul (R[i], R[i], M);
      mpz_mul_2exp (R[i], R[i], i);
      mpz_mod (R[i], R[i], mod->p);
    }
  fermat_vec_import (mod->R_tft, R, mod->d);
  fermat_vec_tft (mod->R_tft, k, 0, mod->d, mod->d, mod->order);
  // include scale factor 2^(-k)
  s = 4 * GMP_NUMB_BITS * mod->m - 2 * mod->k;
  for (i = 0; i < mod->d; i++)
    fermat_mul_sqrt2exp (&mod->R_tft->x[i], &mod->R_tft->t, s, m);

  // compute transform of S
  fermat_vec_import (mod->S_tft, S, mod->d);
  fermat_vec_tft (mod->S_tft, k, 0, mod->d, mod->d, mod->order);
  // include scale factor 2^(-k)
  for (i = 0; i < mod->d; i++)
    fermat_mul_sqrt2exp (&mod->S_tft->x[i], &mod->S_tft->t, s, m);

  // compute transform of g
  fermat_vec_import (mod->g_tft, mod->g, mod->d + 1);
  fermat_vec_tft (mod->g_tft, k, 0, mod->d + 1, mod->d, mod->order);

 cleanup:
  mpz_clear (M);
  mpz_vec_free (P, mod->d + 1);
  mpz_vec_free (Q, mod->d + 1);
  mpz_vec_free (R, mod->d);
  mpz_vec_free (S, mod->d);
  mpz_vec_free (T, mod->d);
  return;

 init_failed:
puts ("Montgomery init failed");
  // fall back on barrett
  mod->algo = ZP_POLY_MOD_ALGO_BARRETT;
  zp_poly_mod_clear_montgomery (mod);
  zp_poly_mod_init_barrett (mod);
  goto cleanup;
}


void zp_poly_mod_mul_montgomery (mpz_t h[], mpz_t f1[], mpz_t f2[], int e,
                 zp_poly_mod_t mod)
{
  mpz_t* t, * v;
  mpz_t z;
  int sqr, i;
  fermat_vec_t vec1, vec2;
  mp_limb_t* temp;
  size_t w, s;

  sqr = (f1 == f2);

  mpz_init2 (z, mod->m * GMP_NUMB_BITS);
  t = mpz_vec_alloc (mod->d);
  v = mpz_vec_alloc (mod->d);
  temp = mp_limb_vec_alloc (2 * mod->m);
  fermat_vec_init (vec1, mod->K, mod->m);
  fermat_vec_init (vec2, mod->K, mod->m);

  // K-th root of unity
  w = (4 * GMP_NUMB_BITS * mod->m) >> mod->k;

  // ================   t = (x^e * f1 * f2 * R)   mod Q

  // substitute x -> 2x in f1, import and transform
  for (i = 0; i < mod->d; i++)
    {
      mpz_mul_2exp (z, f1[i], i);
      mpz_mod (v[i], z, mod->p);
    }
  fermat_vec_import (vec1, v, mod->d);
  fermat_vec_tft (vec1, mod->k, 0, mod->d, mod->d, mod->order);

  if (!sqr)
    {
      // substitute x -> 2x in f2, import and transform
      for (i = 0; i < mod->d; i++)
    {
      mpz_mul_2exp (z, f2[i], i);
      mpz_mod (v[i], z, mod->p);
    }
      fermat_vec_import (vec2, v, mod->d);
      fermat_vec_tft (vec2, mod->k, 0, mod->d, mod->d, mod->order);
    }

  // multiply transforms for f1 * f2 * R, and twist if e == 1
  // (R includes factor of 2^(-d+1) mod p)
  // (R_tft includes 2^(-k) scale factor)
  for (i = 0; i < mod->d; i++)
    {
      fermat_mul (vec1->x[i], vec1->x[i],
          (sqr ? vec1 : vec2)->x[i], temp, mod->m, vec1->best_k);
      fermat_mul (vec1->x[i], vec1->x[i], mod->R_tft->x[i], temp,
          mod->m, vec1->best_k);

      if (e)
    {
      s = w * (bit_reverse (i, mod->k) ^ mod->order) + 2;
      if (s >= 4 * GMP_NUMB_BITS * mod->m)
        s -= 4 * GMP_NUMB_BITS * mod->m;
      fermat_mul_sqrt2exp (&vec1->x[i], &vec1->t, s, mod->m);
    }
    }

  // inverse transform, substitute back x -> x/2 to obtain t
  fermat_vec_itft (vec1, mod->k, 0, mod->d, mod->d, 0, mod->order);
  fermat_vec_export_modp (t, vec1, 0, mod->d, 1, mod->p);
  for (i = 0; i < mod->d; i++)
    {
      mpz_mul_2exp (z, t[i], mod->d - 1 - i);
      mpz_mod (t[i], z, mod->p);
    }

  // ================   (x^e * f1 * f2 - g * t) * S  mod P

  // import f1 and transform
  fermat_vec_import (vec1, f1, mod->d);
  fermat_vec_tft (vec1, mod->k, 0, mod->d, mod->d, mod->order);

  if (!sqr)
    {
      // import f2 and transform
      fermat_vec_import (vec2, f2, mod->d);
      fermat_vec_tft (vec2, mod->k, 0, mod->d, mod->d, mod->order);
    }

  // multiply fourier coefficients and twist if e == 1
  for (i = 0; i < mod->d; i++)
    {
      fermat_mul (vec1->x[i], vec1->x[i],
          (sqr ? vec1 : vec2)->x[i], temp, mod->m, vec1->best_k);
      if (e)
    fermat_mul_sqrt2exp (&vec1->x[i], &vec1->t,
                 w * (bit_reverse (i, mod->k) ^ mod->order),
                 mod->m);
    }

  // import t and transform
  fermat_vec_import (vec2, t, mod->d);
  fermat_vec_tft (vec2, mod->k, 0, mod->d, mod->d, mod->order);

  // multiply by g, subtract, and multiply by S
  // (S_tft includes 2^(-k) scale factor)
  for (i = 0; i < mod->d; i++)
    {
      fermat_mul (vec2->x[i], vec2->x[i], mod->g_tft->x[i], temp,
          mod->m, vec2->best_k);
      fermat_sub (vec1->x[i], vec1->x[i], vec2->x[i], mod->m);
      fermat_mul (vec1->x[i], vec1->x[i], mod->S_tft->x[i], temp,
          mod->m, vec1->best_k);
    }

  // inverse transform and convert back to zp_poly
  fermat_vec_itft (vec1, mod->k, 0, mod->d, mod->d, 0, mod->order);
  fermat_vec_export_modp (h, vec1, 0, mod->d, 1, mod->p);

  // ================

  mp_limb_vec_free (temp, 2 * mod->m);
  mpz_vec_free (v, mod->d);
  mpz_vec_free (t, mod->d);
  mpz_clear (z);
  fermat_vec_clear (vec1);
  fermat_vec_clear (vec2);
}



void zp_poly_mod_import_montgomery (mpz_t h[], mpz_t f[],
                    zp_poly_mod_t mod)
{
  // TODO: optimise this
  zp_poly_mod_t B;
  zp_poly_mod_init (B, mod->g, mod->d, mod->p, ZP_POLY_MOD_ALGO_BARRETT);
  zp_poly_mod_mulx (h, f, mod->Q_mod_g, 0, B);
  zp_poly_mod_clear (B);
}


void zp_poly_mod_import_one_montgomery (mpz_t h[], zp_poly_mod_t mod)
{
  // TODO: optimise this
  mpz_t* one;
  one = mpz_vec_alloc (mod->d);
  mpz_set_ui (one[0], 1);
  zp_poly_mod_import_montgomery (h, one, mod);
  mpz_vec_free (one, mod->d);
}


void zp_poly_mod_export_montgomery (mpz_t h[], mpz_t f[],
                    zp_poly_mod_t mod)
{
  // TODO: optimise this
  mpz_t* one;
  one = mpz_vec_alloc (mod->d);
  mpz_set_ui (one[0], 1);
  zp_poly_mod_mul_montgomery (h, f, one, 0, mod);
  mpz_vec_free (one, mod->d);
}


// ==========================================================================
// dispatchers for zp_poly_mod_t, all PUBLIC


void zp_poly_mod_init (zp_poly_mod_t mod, mpz_t g[], int d, zp_p _p, int algo)
{
  zp_p_to_mpz p=_p;
  int i;

  // copy basic parameters to mod
  mod->d = d;
  mpz_init_set (mod->p, p);
  mod->g = mpz_vec_alloc (d + 1);
  for (i = 0; i <= d; i++)
    mpz_set (mod->g[i], g[i]);
  assert (mpz_cmp_ui (g[d], 1) == 0);

  i = zp_bits(p);
  // select algorithm
  if (algo == ZP_POLY_MOD_ALGO_AUTO)
  {
    if ( i < ZP_POLY_TFT_MIN_PBITS || d/i > ZP_POLY_TFT_MAX_DEG_PBIT_RATIO )
      algo = ZP_POLY_MOD_ALGO_KRONECKER;
    else
      algo = ZP_POLY_MOD_ALGO_BARRETT;
  }
  if ( i < 32 ) algo = ZP_POLY_MOD_ALGO_STUPID;     // force slow mult for small p

  mod->algo = algo;

  // dispatch to appropriate constructor
  switch (algo)
    {
    case ZP_POLY_MOD_ALGO_STUPID:
      zp_poly_mod_init_stupid (mod);
      break;

    case ZP_POLY_MOD_ALGO_BARRETT:
      zp_poly_mod_init_barrett (mod);
      break;

    case ZP_POLY_MOD_ALGO_MONTGOMERY:
      zp_poly_mod_init_montgomery (mod);
      break;

    case ZP_POLY_MOD_ALGO_KRONECKER:
      zp_poly_mod_init_kronecker (mod);
      break;

    default: assert (0);
    }
}


void zp_poly_mod_clear (zp_poly_mod_t mod)
{
  switch (mod->algo)
    {
    case ZP_POLY_MOD_ALGO_STUPID:
      zp_poly_mod_clear_stupid (mod);
      break;

    case ZP_POLY_MOD_ALGO_BARRETT:
      zp_poly_mod_clear_barrett (mod);
      break;

    case ZP_POLY_MOD_ALGO_MONTGOMERY:
      zp_poly_mod_clear_montgomery (mod);
      break;

    case ZP_POLY_MOD_ALGO_KRONECKER:
      zp_poly_mod_clear_kronecker (mod);
      break;

    default: assert (0);
    }

  mpz_vec_free (mod->g, mod->d + 1);
  mpz_clear (mod->p);
}


void zp_poly_mod_import (mpz_t h[], mpz_t f[], zp_poly_mod_t mod)
{
  switch (mod->algo)
    {
    case ZP_POLY_MOD_ALGO_STUPID:
      zp_poly_mod_import_stupid (h, f, mod);
      break;

    case ZP_POLY_MOD_ALGO_BARRETT:
      zp_poly_mod_import_barrett (h, f, mod);
      break;

    case ZP_POLY_MOD_ALGO_MONTGOMERY:
      zp_poly_mod_import_montgomery (h, f, mod);
      break;

    case ZP_POLY_MOD_ALGO_KRONECKER:
      zp_poly_mod_import_kronecker (h, f, mod);
      break;

    default: assert (0);
    }
}


void zp_poly_mod_import_one (mpz_t h[], zp_poly_mod_t mod)
{
  switch (mod->algo)
    {
    case ZP_POLY_MOD_ALGO_STUPID:
      zp_poly_mod_import_one_stupid (h, mod);
      break;

    case ZP_POLY_MOD_ALGO_BARRETT:
      zp_poly_mod_import_one_barrett (h, mod);
      break;

    case ZP_POLY_MOD_ALGO_MONTGOMERY:
      zp_poly_mod_import_one_montgomery (h, mod);
      break;

    case ZP_POLY_MOD_ALGO_KRONECKER:
      zp_poly_mod_import_one_kronecker (h, mod);
      break;

    default: assert (0);
    }
}


void zp_poly_mod_export (mpz_t h[], mpz_t f[], zp_poly_mod_t mod)
{
  switch (mod->algo)
    {
    case ZP_POLY_MOD_ALGO_STUPID:
      zp_poly_mod_export_stupid (h, f, mod);
      break;

    case ZP_POLY_MOD_ALGO_BARRETT:
      zp_poly_mod_export_barrett (h, f, mod);
      break;

    case ZP_POLY_MOD_ALGO_MONTGOMERY:
      zp_poly_mod_export_montgomery (h, f, mod);
      break;

    case ZP_POLY_MOD_ALGO_KRONECKER:
      zp_poly_mod_export_kronecker (h, f, mod);
      break;

    default: assert (0);
    }
}


void zp_poly_mod_mulx (mpz_t h[], mpz_t f1[], mpz_t f2[], int e,
              zp_poly_mod_t mod)
{
  switch (mod->algo)
    {
    case ZP_POLY_MOD_ALGO_STUPID:
      zp_poly_mod_mul_stupid (h, f1, f2, e, mod);
      break;

    case ZP_POLY_MOD_ALGO_BARRETT:
      zp_poly_mod_mul_barrett (h, f1, f2, e, mod);
      break;

    case ZP_POLY_MOD_ALGO_MONTGOMERY:
      zp_poly_mod_mul_montgomery (h, f1, f2, e, mod);
      break;

    case ZP_POLY_MOD_ALGO_KRONECKER:
      zp_poly_mod_mul_kronecker (h, f1, f2, e, mod);
      break;

    default: assert (0);
    }
}

#ifdef ZP_POLY_MUL_TEST


// ==========================================================================
// test code


gmp_randstate_t randstate;


void test_zp_poly_mod_mul_stupid ()
{
  zp_poly_mod_t mod;
  mpz_t p;
  mpz_t t;
  mpz_t* g, * f1, * f2, * f, * h1, * h2, * h;

  printf ("test_zp_poly_mod_mul_stupid...\n");

  mpz_init (t);

  mpz_init (p);
  mpz_set_str (p, "100000007", 10);

  g = mpz_vec_alloc (5);
  f1 = mpz_vec_alloc (4);
  f2 = mpz_vec_alloc (4);
  f = mpz_vec_alloc (4);
  h1 = mpz_vec_alloc (4);
  h2 = mpz_vec_alloc (4);
  h = mpz_vec_alloc (4);

  mpz_set_str (g[0], "1234567", 10);
  mpz_set_str (g[1], "9876543", 10);
  mpz_set_str (g[2], "3141592", 10);
  mpz_set_str (g[3], "2718282", 10);
  mpz_set_ui (g[4], 1);

  zp_poly_mod_init (mod, g, 4, p, ZP_POLY_MOD_ALGO_STUPID);

  mpz_set_str (f1[0], "99999", 10);
  mpz_set_str (f1[1], "35691623", 10);
  mpz_set_str (f1[2], "81786861", 10);
  mpz_set_str (f1[3], "23567234", 10);

  mpz_set_str (f2[0], "78643", 10);
  mpz_set_str (f2[1], "91723", 10);
  mpz_set_str (f2[2], "51762873", 10);
  mpz_set_str (f2[3], "82374873", 10);

  zp_poly_mod_import (h1, f1, mod);
  zp_poly_mod_import (h2, f2, mod);

  zp_poly_mod_mulx (h, h1, h2, 0, mod);
  zp_poly_mod_export (f, h, mod);
  mpz_set_str (t, "72166028", 10);
  assert (mpz_cmp (t, f[0]) == 0);
  mpz_set_str (t, "50314117", 10);
  assert (mpz_cmp (t, f[1]) == 0);
  mpz_set_str (t, "29720232", 10);
  assert (mpz_cmp (t, f[2]) == 0);
  mpz_set_str (t, "63393999", 10);
  assert (mpz_cmp (t, f[3]) == 0);

  zp_poly_mod_mulx (h, h1, h2, 1, mod);
  zp_poly_mod_export (f, h, mod);
  mpz_set_str (t, "66315061", 10);
  assert (mpz_cmp (t, f[0]) == 0);
  mpz_set_str (t, "58928516", 10);
  assert (mpz_cmp (t, f[1]) == 0);
  mpz_set_str (t, "84148776", 10);
  assert (mpz_cmp (t, f[2]) == 0);
  mpz_set_str (t, "75393110", 10);
  assert (mpz_cmp (t, f[3]) == 0);
  
  zp_poly_mod_clear (mod);

  mpz_clear (p);
  mpz_clear (t);
  mpz_vec_free (g, 5);
  mpz_vec_free (f1, 4);
  mpz_vec_free (f2, 4);
  mpz_vec_free (f, 4);
  mpz_vec_free (h1, 4);
  mpz_vec_free (h2, 4);
  mpz_vec_free (h, 4);
}



void test_recip_stupid ()
{
  int i, j, k, d;
  mpz_t p;
  mpz_t* g, * f, * h;

  printf ("test_recip_stupid...\n");

  d = 4;

  mpz_init (p);
  mpz_set_str (p, "100000007", 10);

  g = mpz_vec_alloc (d + 1);
  h = mpz_vec_alloc (d);
  f = mpz_vec_alloc (2*d);

  for (i = 0; i < 10; i++)
    {
      mpz_set_ui (g[d], 1);
      for (j = 0; j < d; j++)
    mpz_urandomm (g[j], randstate, p);

      recip_stupid (h, g, d, p);

      // multiply out
      for (j = 0; j < 2*d; j++) mpz_set_ui (f[j], 0);
      for (j = 0; j <= d; j++) for (k = 0; k < d; k++) mpz_addmul (f[j+k], g[j], h[k]);
      for (j = 0; j < 2*d; j++) mpz_mod (f[j], f[j], p);

      // check answer
      assert (mpz_cmp_ui (f[2*d - 1], 1) == 0);
      for (j = d; j < 2*d - 1; j++)
    assert (mpz_cmp_ui (f[j], 0) == 0);
    }

  mpz_vec_free (g, d + 1);
  mpz_vec_free (h, d);
  mpz_vec_free (f, 2*d);
  mpz_clear (p);
}



void testcase_fft (int k, int K, size_t m, int order)
{
  int i;
  mpz_t p;
  zp_poly_mod_t mod;
  mpz_t* g, * f1, * f2, * h, * r;
  fermat_vec_t vec1, vec2;

  // p = B^m + 1
  mpz_init_set_ui (p, 1);
  mpz_mul_2exp (p, p, m * GMP_NUMB_BITS);
  mpz_add_ui (p, p, 1);

  // set up for multiplication mod x^K - 1 mod p
  g = mpz_vec_alloc (K + 1);
  mpz_set (g[0], p);
  mpz_sub_ui (g[0], g[0], 1);
  mpz_set_ui (g[K], 1);
  zp_poly_mod_init (mod, g, K, p, ZP_POLY_MOD_ALGO_STUPID);

  // make up random polys, convert to fermat polys
  f1 = mpz_vec_alloc (K);
  f2 = mpz_vec_alloc (K);
  fermat_vec_init (vec1, K, m);
  fermat_vec_init (vec2, K, m);
  for (i = 0; i < K; i++)
    {
      mpz_rrandomb (f1[i], randstate, (m + 1) * GMP_NUMB_BITS);
      mpz_mod (f1[i], f1[i], p);
      fermat_import (vec1->x[i], f1[i], m);

      mpz_rrandomb (f2[i], randstate, (m + 1) * GMP_NUMB_BITS);
      mpz_mod (f2[i], f2[i], p);
      fermat_import (vec2->x[i], f2[i], m);
    }

  // do convolution using zp_poly_mod
  h = mpz_vec_alloc (K);
  r = mpz_vec_alloc (K);
  zp_poly_mod_mulx (h, f1, f2, 0, mod);

  // do convolution using fermat_vec
  mp_limb_t* temp = mp_limb_vec_alloc (2 * m);
  fermat_vec_fft (vec1, k, order);
  fermat_vec_fft (vec2, k, order);

  for (i = 0; i < K; i++)
    fermat_mul (vec1->x[i], vec1->x[i], vec2->x[i], temp, m, vec1->best_k);

  fermat_vec_ifft (vec1, k, order);

  // scale by 1/K and normalise
  for (i = 0; i < K; i++)
    {
      if (k > 0)
    fermat_mul_sqrt2exp (&vec1->x[i], &vec1->t,
                 4 * m * GMP_NUMB_BITS - 2 * k, m);
      fermat_normalise (vec1->x[i], m);
    }

  for (i = 0; i < K; i++)
    fermat_export (r[i], vec1->x[i], vec1->t, m);

  // check results agree
  for (i = 0; i < K; i++)
    assert (mpz_cmp (r[i], h[i]) == 0);

  mp_limb_vec_free (temp, 2 * m);
  fermat_vec_clear (vec2);
  fermat_vec_clear (vec1);
  zp_poly_mod_clear (mod);
  mpz_vec_free (g, K + 1);
  mpz_vec_free (f1, K);
  mpz_vec_free (f2, K);
  mpz_vec_free (h, K);
  mpz_vec_free (r, K);
  mpz_clear (p);
}


void test_fft ()
{
  // perform some convolutions, check we get the right result
  int k, K, order;
  size_t m, mdiv;

  printf ("test_fft...\n");

  for (k = 0; k <= 11; k++)
    {
      K = 1 << k;
      mdiv = get_mdiv (K);
      for (m = mdiv; m < 5*mdiv; m += mdiv)
    {
      order = random () % K;
      testcase_fft (k, K, m, order);
    }
    }
}



void testcase_tft (int k, int K, size_t m, size_t y, int z, int n, int order)
{
  int i;
  mpz_t p, x, x2;
  fermat_vec_t vec1, vec2;

  mpz_init (x);
  mpz_init (x2);

  // p = B^m + 1
  mpz_init_set_ui (p, 1);
  mpz_mul_2exp (p, p, m * GMP_NUMB_BITS);
  mpz_add_ui (p, p, 1);

  // make up random poly
  fermat_vec_init (vec1, K, m);
  for (i = 0; i < K; i++)
    {
      mpz_rrandomb (x, randstate, (m + 1) * GMP_NUMB_BITS);
      mpz_mod (x, x, p);
      fermat_import (vec1->x[i], x, m);
    }

  // copy first n coeffs to another poly and zero-pad
  fermat_vec_init (vec2, K, m);
  for (i = 0; i < z; i++)
    fermat_set (vec2->x[i], vec1->x[i], m);
  for (; i < K; i++)
    fermat_zero (vec2->x[i], m);

  // run both FFT and TFT and compare results
  fermat_vec_fft (vec2, k, order);
  fermat_vec_tft (vec1, k, y, z, n, order);

  for (i = 0; i < n; i++)
    {
      fermat_mul_sqrt2exp (&vec2->x[i], &vec2->t,
               y * (bit_reverse (i, k) ^ order), m);
      fermat_export (x, vec1->x[i], vec1->t, m);
      fermat_export (x2, vec2->x[i], vec2->t, m);
      assert (mpz_cmp (x, x2) == 0);
    }

  mpz_clear (x);
  mpz_clear (x2);
  mpz_clear (p);
  fermat_vec_clear (vec1);
  fermat_vec_clear (vec2);
}



void test_tft ()
{
  // compare against non-truncated fft

  int i, k, K, z, n, order;
  size_t m, mdiv, y;

  printf ("test_tft...\n");

  for (i = 0; i < 100000; i++)
    {
      k = random () % 10;
      K = 1 << k;
      mdiv = get_mdiv (K);

      m = ((random () % 7) + 1) * mdiv;
      y = random () % (4 * GMP_NUMB_BITS * m / K);
      z = (random () % (K + 1));
      n = (random () % K) + 1;
      order = random () % K;

      testcase_tft (k, K, m, y, z, n, order);
    }
}



void testcase_itft (int k, int K, size_t m, size_t y, int z, int n, int f,
            int order)
{
  int i;
  mpz_t p, x, x2;
  fermat_vec_t vec1, vec2, vec3;

  mpz_init (x);
  mpz_init (x2);

  // p = B^m + 1
  mpz_init_set_ui (p, 1);
  mpz_mul_2exp (p, p, m * GMP_NUMB_BITS);
  mpz_add_ui (p, p, 1);

  // make up random poly of length z
  fermat_vec_init (vec1, K, m);
  for (i = 0; i < z; i++)
    {
      mpz_rrandomb (x, randstate, (m + 1) * GMP_NUMB_BITS);
      mpz_mod (x, x, p);
      fermat_import (vec1->x[i], x, m);
    }
  for (; i < K; i++)
    fermat_zero (vec1->x[i], m);

  // compute its full FFT
  fermat_vec_init (vec2, K, m);
  for (i = 0; i < K; i++)
    fermat_set (vec2->x[i], vec1->x[i], m);
  fermat_vec_fft (vec2, k, order);
  for (i = 0; i < K; i++)
    fermat_mul_sqrt2exp (&vec2->x[i], &vec2->t,
             y * (bit_reverse (i, k) ^ order), m);

  for (i = 0; i < K; i++)
    fermat_mul_sqrt2exp (&vec1->x[i], &vec1->t, 2 * k, m);

  // run itft
  fermat_vec_init (vec3, K, m);
  for (i = 0; i < n; i++)
    fermat_set (vec3->x[i], vec2->x[i], m);
  for (; i < z; i++)
    fermat_set (vec3->x[i], vec1->x[i], m);
  for (; i < K; i++)
    {
      mpz_rrandomb (x, randstate, (m + 1) * GMP_NUMB_BITS);
      mpz_mod (x, x, p);
      fermat_import (vec3->x[i], x, m);
    }
  fermat_vec_itft (vec3, k, y, z, n, f, order);

  // check output
  for (i = 0; i < n; i++)
    {
      fermat_export (x, vec1->x[i], vec1->t, m);
      fermat_export (x2, vec3->x[i], vec3->t, m);
      assert (mpz_cmp (x, x2) == 0);
    }
  if (f)
    {
      fermat_export (x, vec2->x[n], vec2->t, m);
      fermat_export (x2, vec3->x[n], vec3->t, m);
      assert (mpz_cmp (x, x2) == 0);
    }

  mpz_clear (x);
  mpz_clear (x2);
  mpz_clear (p);
  fermat_vec_clear (vec1);
  fermat_vec_clear (vec2);
  fermat_vec_clear (vec3);
}



void test_itft ()
{
  // compare against truncated fft

  int i, k, K, z, n, f, order;
  size_t m, mdiv, y;

  printf ("test_itft...\n");

  for (i = 0; i < 100000; i++)
    {
      int good;

      do
    {
      k = random () % 10;
      K = 1 << k;
      f = random () % 2;
      z = (random () % K) + 1;
      n = (random () % K) + 1 - f;
      good = (z >= n);
    }
      while (!good);

      mdiv = get_mdiv (K);
      m = ((random () % 7) + 1) * mdiv;
      y = random () % (4 * GMP_NUMB_BITS * m / K);
      order = random () % K;

      testcase_itft (k, K, m, y, z, n, f, order);
    }
}


void testcase_zp_poly_mod_mul_barrett (int d, int e, int sqr, mpz_t p)
{
  int i;
  mpz_t* f1, * f2, * g, * h1, * h2;
  zp_poly_mod_t mod1, mod2;

  f1 = mpz_vec_alloc (d);
  f2 = mpz_vec_alloc (d);
  g = mpz_vec_alloc (d + 1);
  h1 = mpz_vec_alloc (d);
  h2 = mpz_vec_alloc (d);

  mpz_set_ui (g[d], 1);
  for (i = 0; i < d; i++)
    {
      mpz_urandomm (f1[i], randstate, p);
      mpz_urandomm (f2[i], randstate, p);
      mpz_urandomm (g[i], randstate, p);
    }

  zp_poly_mod_init (mod1, g, d, p, ZP_POLY_MOD_ALGO_STUPID);
  zp_poly_mod_init (mod2, g, d, p, ZP_POLY_MOD_ALGO_BARRETT);

  if (sqr)
    {
      zp_poly_mod_mulx (h1, f1, f1, e, mod1);
      zp_poly_mod_mulx (h2, f1, f1, e, mod2);
    }
  else
    {
      zp_poly_mod_mulx (h1, f1, f2, e, mod1);
      zp_poly_mod_mulx (h2, f1, f2, e, mod2);
    }

  for (i = 0; i < d; i++)
    assert (mpz_cmp (h1[i], h2[i]) == 0);

  zp_poly_mod_clear (mod1);
  zp_poly_mod_clear (mod2);
  mpz_vec_free (f1, d);
  mpz_vec_free (f2, d);
  mpz_vec_free (g, d + 1);
  mpz_vec_free (h1, d);
  mpz_vec_free (h2, d);
}


void test_zp_poly_mod_mul_barrett ()
{
  int i, d, e, bits, sqr;
  mpz_t p;

  printf ("test_zp_poly_mod_mul_barrett...\n");

  mpz_init (p);

  for (i = 0; i < 1000; i++)
    {
      bits = random () % 1000 + 64;
      e = random () % 2;
      sqr = random () % 2;
      d = random () % 100 + 1;

      mpz_urandomb (p, randstate, bits);
      mpz_nextprime (p, p);

      testcase_zp_poly_mod_mul_barrett (d, e, sqr, p);
    }

  mpz_clear (p);
}



void testcase_zp_poly_mod_mul_kronecker (int d, int e, int sqr, mpz_t p)
{
  int i;
  mpz_t* f1, * f2, * g, * h1, * h2;
  zp_poly_mod_t mod1, mod2;
    
  f1 = mpz_vec_alloc (d);
  f2 = mpz_vec_alloc (d);
  g = mpz_vec_alloc (d + 1);
  h1 = mpz_vec_alloc (d);
  h2 = mpz_vec_alloc (d);

  mpz_set_ui (g[d], 1);
  for (i = 0; i < d; i++)
    {
      mpz_urandomm (f1[i], randstate, p);
      mpz_urandomm (f2[i], randstate, p);
      mpz_urandomm (g[i], randstate, p);
    }
  zp_poly_mod_init (mod1, g, d, p, ZP_POLY_MOD_ALGO_STUPID);
  zp_poly_mod_init (mod2, g, d, p, ZP_POLY_MOD_ALGO_KRONECKER);

  if (sqr)
    {
      zp_poly_mod_mulx (h1, f1, f1, e, mod1);
      zp_poly_mod_mulx (h2, f1, f1, e, mod2);
    }
  else
    {
      zp_poly_mod_mulx (h1, f1, f2, e, mod1);
      zp_poly_mod_mulx (h2, f1, f2, e, mod2);
    }

  for (i = 0; i < d; i++)
    assert (mpz_cmp (h1[i], h2[i]) == 0);

  zp_poly_mod_clear (mod1);
  zp_poly_mod_clear (mod2);
  mpz_vec_free (f1, d);
  mpz_vec_free (f2, d);
  mpz_vec_free (g, d + 1);
  mpz_vec_free (h1, d);
  mpz_vec_free (h2, d);
}


void test_zp_poly_mod_mul_kronecker ()
{
  int i, d, e, bits, sqr;
  mpz_t p;

  printf ("test_zp_poly_mod_mul_kronecker...\n");

  mpz_init (p);

  for (i = 0; i < 1000; i++)
    {
      bits = random () % 500 + 64;
      e = random () % 2;
      sqr = random () % 2;
      d = random () % 500 + 1;

      mpz_urandomb (p, randstate, bits);
      mpz_nextprime (p, p);

      testcase_zp_poly_mod_mul_kronecker (d, e, sqr, p);
    }

  mpz_clear (p);
}


void testcase_zp_poly_mod_mul_montgomery (int d, int e, int sqr, mpz_t p)
{
  int i;
  mpz_t* f1, * f2, * g, * h1, * h2, * h3, * u1 ,* u2;
  zp_poly_mod_t mod1, mod2;
    
  f1 = mpz_vec_alloc (d);
  f2 = mpz_vec_alloc (d);
  u1 = mpz_vec_alloc (d);
  u2 = mpz_vec_alloc (d);
  g = mpz_vec_alloc (d + 1);
  h1 = mpz_vec_alloc (d);
  h2 = mpz_vec_alloc (d);
  h3 = mpz_vec_alloc (d);

  mpz_set_ui (g[d], 1);
  for (i = 0; i < d; i++)
    {
      mpz_urandomm (f1[i], randstate, p);
      mpz_urandomm (f2[i], randstate, p);
      mpz_urandomm (g[i], randstate, p);
    }

  zp_poly_mod_init (mod1, g, d, p, ZP_POLY_MOD_ALGO_STUPID);
  zp_poly_mod_init (mod2, g, d, p, ZP_POLY_MOD_ALGO_MONTGOMERY);
  zp_poly_mod_import (u1, f1, mod2);
  zp_poly_mod_import (u2, f2, mod2);
  if (sqr)
    {
      zp_poly_mod_mulx (h1, f1, f1, e, mod1);
      zp_poly_mod_mulx (h2, u1, u1, e, mod2);
    }
  else
    {
      zp_poly_mod_mulx (h1, f1, f2, e, mod1);
      zp_poly_mod_mulx (h2, u1, u2, e, mod2);
    }

  zp_poly_mod_export (h3, h2, mod2);

  for (i = 0; i < d; i++)
    assert (mpz_cmp (h1[i], h3[i]) == 0);

  zp_poly_mod_clear (mod2);
  zp_poly_mod_clear (mod1);
  mpz_vec_free (f1, d);
  mpz_vec_free (f2, d);
  mpz_vec_free (u1, d);
  mpz_vec_free (u2, d);
  mpz_vec_free (g, d + 1);
  mpz_vec_free (h1, d);
  mpz_vec_free (h2, d);
  mpz_vec_free (h3, d);
}

// TODO: there is a bug in the montgomery code that shows up in small cases
void test_zp_poly_mod_mul_montgomery ()
{
  int i, d, e, bits, sqr;
  mpz_t p;

  printf ("test_zp_poly_mod_mul_montgomery...\n");

  mpz_init (p);

  for (i = 0; i < 1000; i++)
    {
      bits = random () % 1000 + 64; // changed 1000 to 100 to expose bug
      e = random () % 2;
      sqr = random () % 2;
      d = random () % 100 + 1;      // changed 100 to 10 to expose bug

      mpz_urandomb (p, randstate, bits);
      mpz_nextprime (p, p);

      testcase_zp_poly_mod_mul_montgomery (d, e, sqr, p);
    }

  mpz_clear (p);
}


void testcase_zp_array_mul (int d1, int d2, int sqr, mpz_t p)
{
  int i, j, d3;
  mpz_t* f1, * f2, * h1, * h2;

  f1 = mpz_vec_alloc (d1);
  for (i = 0; i < d1; i++)
    mpz_urandomm (f1[i], randstate, p);

  if (sqr)
    {
      f2 = f1;
      d2 = d1;
    }
  else
    {
      f2 = mpz_vec_alloc (d2);
      for (i = 0; i < d2; i++)
    mpz_urandomm (f2[i], randstate, p);
    }

  d3 = d1 + d2 - 1;
  h1 = mpz_vec_alloc (d3);
  h2 = mpz_vec_alloc (d3);

  // multiply out stupidly
  for (i = 0; i < d1; i++)
    for (j = 0; j < d2; j++)
      mpz_addmul (h1[i+j], f1[i], f2[j]);
  for (i = 0; i < d3; i++)
    mpz_mod (h1[i], h1[i], p);

  // multiply with zp_array_mul
  zp_array_mul (h2, f1, d1, f2, d2, p);

  // compare results
  for (i = 0; i < d3; i++)
    assert (mpz_cmp (h1[i], h2[i]) == 0);

  mpz_vec_free (f1, d1);
  if (!sqr)
    mpz_vec_free (f2, d2);
  mpz_vec_free (h1, d3);
  mpz_vec_free (h2, d3);
}


void test_zp_array_mul ()
{
  mpz_t p;
  int i, d1, d2, sqr;
  size_t bits;

  printf ("test_zp_array_mul...\n");

  mpz_init (p);

  for (i = 0; i < 1000; i++)
    {
      bits = random () % 500 + 64;
      sqr = random () % 2;

      mpz_urandomb (p, randstate, bits);
      mpz_nextprime (p, p);

      d1 = random () % 500 + 1;
      d2 = random () % 500 + 1;

      testcase_zp_array_mul (d1, d2, sqr, p);
    }

  mpz_clear (p);
}

void testcase_zp_array_mul2n (int d1, int d2, int d3, int d4, int n,
                    int sqr1, int sqr2, mpz_t p)
{
  int d, i, j;
  mpz_t *f1, *f2, *f3, *f4, *h1, *h2, *r, *s;

  // random polys f1, f2, f3, f4
  f1 = mpz_vec_alloc (d1);
  for (i = 0; i < d1; i++)
    mpz_urandomm (f1[i], randstate, p);
  if (sqr1)
    {
      f2 = f1;
      d2 = d1;
    }
  else
    {
      f2 = mpz_vec_alloc (d2);
      for (i = 0; i < d2; i++)
    mpz_urandomm (f2[i], randstate, p);
    }

  f3 = mpz_vec_alloc (d3);
  for (i = 0; i < d3; i++)
    mpz_urandomm (f3[i], randstate, p);
  if (sqr2)
    {
      f4 = f3;
      d4 = d3;
    }
  else
    {
      f4 = mpz_vec_alloc (d4);
      for (i = 0; i < d4; i++)
    mpz_urandomm (f4[i], randstate, p);
    }


  // multiply out stupidly
  d = zp_imax (d1 + d2, d3 + d4) - 1;
  h1 = mpz_vec_alloc (d);
  for (i = 0; i < d1; i++)
    for (j = 0; j < d2; j++)
      mpz_addmul (h1[i+j], f1[i], f2[j]);
  for (i = 0; i < d; i++)
    mpz_mod (h1[i], h1[i], p);
  h2 = mpz_vec_alloc (d);
  for (i = 0; i < d3; i++)
    for (j = 0; j < d4; j++)
      mpz_addmul (h2[i+j], f3[i], f4[j]);
  for (i = 0; i < d; i++)
    mpz_mod (h2[i], h2[i], p);
  r = mpz_vec_alloc(n);
  for (i = 0; i < d && i < n; i++)
    { mpz_add(r[i],h1[i],h2[i]);  mpz_mod(r[i],r[i],p); }

  s = mpz_vec_alloc (n);
  zp_array_mul2_mod_xn (s, f1, d1, f2, d2, f3, d3, f4, d4, n, p);
  
  // compare results
  for (i = 0; i < n; i++)
    assert (mpz_cmp (r[i], s[i]) == 0);
  
  mpz_vec_free (r, n);
  mpz_vec_free (s, n);
  mpz_vec_free (h1, d);
  mpz_vec_free (h2, d);
  mpz_vec_free (f1, d1);
  if (!sqr1)
    mpz_vec_free (f2, d2);
  mpz_vec_free (f3, d3);
  if (!sqr2)
    mpz_vec_free (f4, d4);
}


void test_zp_array_mul2_mod_xn ()
{
  mpz_t p;
  int i, d1, d2, d3, d4, n, sqr1, sqr2;
  size_t bits;

  printf ("test_zp_array_mul2_mod_xn...\n");

  mpz_init (p);

  for (i = 0; i < 1000; i++)
    {
      bits = random () % 1000 + 64;
      sqr1 = random () % 2;
      sqr2 = random () % 2;

      d1 = random () % 100 + 1;
      d2 = random () % 100 + 1;
      d3 = random () % 100 + 1;
      d4 = random () % 100 + 1;
      n = random () % 100 + 1;
      mpz_urandomb (p, randstate, bits);     
      testcase_zp_array_mul2n (d1, d2, d3, d4, n, sqr1, sqr2, p);
    }

  mpz_clear (p);
}



int main (int argc, char* argv[])
{
  gmp_randinit_default (randstate);
  mp_set_memory_functions (zp_malloc, zp_realloc, zp_free);

  test_zp_array_mul ();  zp_mem_report (0);
  test_zp_array_mul2_mod_xn (); zp_mem_report(0);
  test_zp_poly_mod_mul_stupid ();  zp_mem_report (0);
  test_recip_stupid ();  zp_mem_report (0);
  test_fft ();  zp_mem_report (0);
  test_tft ();  zp_mem_report (0);
  test_itft ();  zp_mem_report (0);
  test_zp_poly_mod_mul_barrett ();  zp_mem_report (0);
  test_zp_poly_mod_mul_kronecker ();  zp_mem_report (0);
  test_zp_poly_mod_mul_montgomery ();  zp_mem_report (0);

  gmp_randclear (randstate);
  return 0;
}

#endif
