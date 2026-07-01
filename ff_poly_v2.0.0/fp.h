#ifndef _FP_INCLUDE_
#define _FP_INCLUDE_

/*
    Copyright 2018 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <gmp.h>
#include "polyparse.h"
#include "ff_poly/mm.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef mm_t fp_t;

typedef struct fp_ctx_struct {
    mm_t p;         // odd prime less than or equal to FF_MAX_P
    mm_t pinv;      // -1/p mod B, where B=2^MM_BITS
    mm_t R;         // B mod p (1 in Montgomery rep)
    mm_t R2;        // B^2 mod p
    mm_t R3;        // B^3 mod p
    mm_t half;      // 1/2 mod p
} fp_ctx_t[1];


static inline struct fp_ctx_struct *fp_init (fp_ctx_t ctx, long p)
{
    ctx->p = p; ctx->pinv = mm_pinv(ctx->p);
    ctx->R = mm_R (ctx->p);  ctx->R2 = mm_R2(ctx->R,ctx->p,ctx->pinv); ctx->R3 = mm_R3 (ctx->R2,ctx->p,ctx->pinv);
    ctx->half = mm_div2 (ctx->R, ctx->p);
    return ctx;
}

/*
    This is the recommended interface to the ff_poly library.
*/
/*
typedef struct {
    mm_t p;         // odd prime less than or equal to FF_MAX_P
    mm_t pinv;      // -1/p mod B, where B=2^MM_BITS
    mm_t R;         // B mod p (1 in Montgomery rep)
    mm_t R2;        // B^2 mod p
    mm_t R3;        // B^3 mod p
    mm_t g2;        // generator for 2-sylow subgroup of Fp* (also a non-residue)
    mm_t g2inv;     // 1/g2 mod p
    mm_t half;      // B/2 mod p (1/2 in Montgomery rep)
    mm_t third;     // B/3 mod p (1/3 in Montgomery rep)
    mm_t fourth;    // B/4 mod p (1/4 in Montgomery rep)
    mm_t fifth;     // B/5 mod p (1/5 in Montgomery rep)
    mm_t sixth;     // B/6 mod p (1/6 in Montgomery rep)
    mm_t seventh;   // B/7 mod p (1/7 in Montgomery rep)
    mm_t z3;        // primitive cube root of unity or 1 if p=2mod3
    int v2;         // 2-adic valuation of p-1
    int v3;         // 3-adic valuation of p-1
    mm_t m2;        // (p-1)/2^v2   // integer
    mm_t m3;        // (p-1)/3^v3   // integer
} fffctx_t[1];
*/
static inline long fp_char (fp_ctx_t ctx) { return ctx->p; }
static inline fp_t fp_zero (fp_ctx_t ctx) { return 0; }
static inline fp_t fp_one (fp_ctx_t ctx) { return ctx->R; }
static inline fp_t fp_negone (fp_ctx_t ctx) { return ctx->p-ctx->R; }
static inline fp_t fp_one_half (fp_ctx_t ctx) { return ctx->half; }

static inline fp_t fp_from_int (long n, fp_ctx_t ctx) { return mm_from_si(n,ctx->R2,ctx->p,ctx->pinv); }
static inline long fp_to_int (fp_t a, fp_ctx_t ctx) { return mm_to_ui(a,ctx->p,ctx->pinv); }
static inline fp_t fp_from_mpz (mpz_t n, fp_ctx_t ctx) { return mm_from_si(mpz_get_si(n),ctx->R2,ctx->p,ctx->pinv); }
static inline void fp_to_mpz (mpz_t A, fp_t a, fp_ctx_t ctx) { mpz_set_ui (A, fp_to_int(a, ctx)); }

static inline fp_t fp_is_zero (fp_t a, fp_ctx_t ctx) { return !a; }
static inline fp_t fp_is_nonzero (fp_t a, fp_ctx_t ctx) { return a; }
static inline fp_t fp_is_one (fp_t a, fp_ctx_t ctx) { return a == ctx->R; }

static inline fp_t fp_neg (fp_t a, fp_ctx_t ctx) { return mm_neg (a, ctx->p); }
static inline fp_t fp_inc (fp_t a, fp_ctx_t ctx) { return mm_add (a, ctx->R, ctx->p); }
static inline fp_t fp_dec (fp_t a, fp_ctx_t ctx) { return mm_sub (a, ctx->R, ctx->p); }
static inline fp_t fp_add (fp_t a, fp_t b, fp_ctx_t ctx) { return mm_add (a, b, ctx->p); }
static inline fp_t fp_sub (fp_t a, fp_t b, fp_ctx_t ctx) { return mm_sub (a, b, ctx->p); }
static inline fp_t fp_mul2 (fp_t a, fp_ctx_t ctx) { return mm_add (a, a, ctx->p); }
static inline fp_t fp_div2 (fp_t a, fp_ctx_t ctx) { return mm_div2 (a, ctx->p); }
static inline fp_t fp_mul (fp_t a, fp_t b, fp_ctx_t ctx) { return mm_mul (a, b, ctx->p, ctx->pinv); }
static inline fp_t fp_mul_ui (fp_t a, unsigned long b, fp_ctx_t ctx) { return mm_mul_ui (a, b, ctx->R, ctx->p, ctx->pinv); }
static inline fp_t fp_addmul (fp_t a, fp_t b, fp_t c, fp_ctx_t ctx) { return mm_addmul (a, b, c, ctx->R, ctx->p, ctx->pinv); }
static inline fp_t fp_exp_ui (fp_t a, unsigned long e, fp_ctx_t ctx) { return mm_exp_ui (a, e, ctx->R, ctx->p, ctx->pinv); }

static inline fp_t fp_sqr (fp_t a, fp_ctx_t ctx) { return mm_sqr (a, ctx->p, ctx->pinv); }
static inline fp_t fp_inv (fp_t a, fp_ctx_t ctx) { return mm_inv (a, ctx->R2, ctx->R3, ctx->p, ctx->pinv); }
static inline int fp_legendre (fp_t a, fp_ctx_t ctx) { return mm_legendre (a, ctx->R, ctx->p, ctx->pinv); }
static inline int fp_is_square (fp_t a, fp_ctx_t ctx) { return fp_legendre (a, ctx) >= 0; }

static inline fp_t fp_dot_2r (fp_t a1, fp_t a2, fp_t b1, fp_t b2, fp_ctx_t ctx)
    { return mm_dot_2r (a1, a2, b1, b2, ctx->p, ctx->pinv); }
static inline fp_t fp_conv_2r (fp_t a1, fp_t a2, fp_t b1, fp_t b2, fp_ctx_t ctx)
    { return mm_conv_2r (a1, a2, b1, b2, ctx->p, ctx->pinv); }
    
static inline fp_t fp_dot (fp_t a[], fp_t b[], int n, fp_ctx_t ctx) { return mm_dot (a, b, n, ctx->p, ctx->pinv); }
static inline fp_t fp_dot_ui (fp_t a[], unsigned long b[], int n, fp_ctx_t ctx) { return mm_dot_ui (a, b, n, ctx->R, ctx->p, ctx->pinv); }
static inline fp_t fp_conv (fp_t a[], fp_t b[], int n, fp_ctx_t ctx) { return mm_conv (a, b, n, ctx->p, ctx->pinv); }
      
static inline int fp_poly_degree (fp_t f[], int d, fp_ctx_t ctx)
    { while ( d >= 0 && fp_is_zero(f[d],ctx) ) d--; return d; }
static inline fp_t *fp_poly_from_int (fp_t f[], long F[], int d, fp_ctx_t ctx) { for ( int i = 0 ; i <= d ; i++ ) f[i] = fp_from_int (F[i], ctx); return f; }
static inline fp_t *fp_poly_from_mpz (fp_t f[], mpz_t F[], int d, fp_ctx_t ctx) { for ( int i = 0 ; i <= d ; i++ ) f[i] = fp_from_mpz (F[i], ctx); return f; }
static inline long *fp_poly_to_int (long F[], fp_t f[], int d, fp_ctx_t ctx) { for ( int i = 0 ; i <= d ; i++ ) F[i] = fp_to_int (f[i], ctx); return F; }
static inline mpz_t *fp_poly_to_mpz (mpz_t F[], fp_t f[], int d, fp_ctx_t ctx) { for ( int i = 0 ; i <= d ; i++ ) fp_to_mpz (F[i], f[i], ctx);  return F; }
static inline void fp_poly_print (fp_t f[], int d, fp_ctx_t ctx)
    { int64_t o[d+1]; i_poly_print (fp_poly_to_int (o, f, d, ctx), d); }

static inline fp_t fp_poly_eval (fp_t f[], int d, fp_t x, fp_ctx_t ctx)
    { return mm_poly_eval (f, d, x, ctx->R, ctx->p, ctx->pinv); }
static inline fp_t *fp_poly_translate (fp_t g[], fp_t f[], int d, fp_t a, fp_ctx_t ctx)
    { return mm_poly_translate (g, f, d, a, ctx->R, ctx->p, ctx->pinv); }
static inline fp_t *fp_poly_make_monic (fp_t o[], fp_t f[], int d, fp_ctx_t ctx)
    { return mm_poly_make_monic (o, f, d, ctx->R, ctx->R2, ctx->R3, ctx->p, ctx->pinv); }
static inline fp_t *fp_poly_depress_monic (fp_t o[], fp_t a[1], fp_t f[], int d, fp_ctx_t ctx)
    { return mm_poly_depress_monic (o, a, f, d, ctx->R, ctx->R2, ctx->R3, ctx->p, ctx->pinv); }
static inline fp_t *fp_poly_depress_and_make_monic (fp_t o[], fp_t a[1], fp_t f[], int d, fp_ctx_t ctx)
    { return mm_poly_depress_and_make_monic (o, a, f, d, ctx->R, ctx->R2, ctx->R3, ctx->p, ctx->pinv); }

static inline fp_t fp_poly_disc (fp_t f[], int d, fp_ctx_t ctx)
    { return mm_poly_disc (f, d, ctx->R, ctx->R2, ctx->R3, ctx->p, ctx->pinv); }    

static inline int fp_poly_is_squarefree (fp_t f[], int d, fp_ctx_t ctx)
    { return mm_poly_squarefree (f, d, ctx->R, ctx->R2, ctx->R3, ctx->p, ctx->pinv); }
    
#ifdef __cplusplus
}
#endif

#endif
