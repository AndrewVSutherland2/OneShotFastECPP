#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "zn_poly/zn_poly.h"
#include "zn_poly/zn_poly_internal.h"
#include "ff.h"
#include "ffpoly.h"
#include "ffpolysmall.h"
#include "ffpolybig.h"
#include "cstd.h"

/*
    Copyright 2008-2016 Andrew V. Sutherland
    See LICENSE file for license details.
*/

static zn_mod_t _ff_poly_zn_mod_ctx;
static unsigned long _ff_poly_zn_mod_p;
static int _ff_poly_zn_mod_init;
static ff_t *_ff_poly_big_a;
static ff_t *_ff_poly_big_b;
static int _ff_poly_alloc_len;

void ff_poly_zn_mod_setup (unsigned long p, int d)
{
    if ( p != _ff_poly_zn_mod_p ) {
        if ( _ff_poly_zn_mod_init ) zn_mod_clear(_ff_poly_zn_mod_ctx);
        zn_mod_init(_ff_poly_zn_mod_ctx, _ff_p);
        _ff_poly_zn_mod_p = _ff_p;
        _ff_poly_zn_mod_init = 1;
    }
    if ( d >= _ff_poly_alloc_len ) {
        _ff_poly_alloc_len = (1<<ui_len(d+1));
        _ff_poly_big_a = realloc(_ff_poly_big_a,_ff_poly_alloc_len*sizeof(ff_t));
        _ff_poly_big_b = realloc(_ff_poly_big_b,_ff_poly_alloc_len*sizeof(ff_t));
    }
}

void ff_poly_zn_mod_clear ()
{
    if ( _ff_poly_zn_mod_init ) zn_mod_clear (_ff_poly_zn_mod_ctx);
    _ff_poly_zn_mod_init = 0;
    _ff_poly_zn_mod_p = 0;
    if ( _ff_poly_big_a ) free (_ff_poly_big_a);
    if ( _ff_poly_big_b ) free (_ff_poly_big_b);
    _ff_poly_alloc_len = 0;
    _ff_poly_big_a = 0;
    _ff_poly_big_b = 0;
}

void ff_poly_square_big (ff_t h[], ff_t f[], int d)
{
    register int i;

    ff_poly_zn_mod_setup(_ff_p, d);
    for ( i = 0 ; i <= d ; i++ ) _ff_poly_big_a[i] = _ff_get_ui(f[i]);
    zn_array_mul (h, _ff_poly_big_a, d+1, _ff_poly_big_a, d+1, _ff_poly_zn_mod_ctx);
    for ( i = 0 ; i <= 2*d ; i++ ) _ff_set_ui(h[i],h[i]);
}

void *ff_poly_middle_product_setup1 (ff_t a[], int d_a, int d_b)
{
    zn_array_mulmid_precomp1_t *ctx;

    register int i;
    
    assert (d_a >= d_b);
    ff_poly_zn_mod_setup(_ff_p, d_a);
    ctx = malloc(sizeof(*ctx));
    for ( i = 0 ; i <= d_a ; i++ ) _ff_poly_big_a[i] = _ff_get_ui(a[i]);
    zn_array_mulmid_precomp1_init (*ctx, _ff_poly_big_a, d_a+1, d_b+1, _ff_poly_zn_mod_ctx);
    return ctx;
}

void ff_poly_middle_product_execute1 (ff_t c[], ff_t b[], void *ctx)
{
    zn_array_mulmid_precomp1_t *mpctx = ctx;
    register int i;
    
    for ( i = 0 ; i < (*mpctx)->n2 ; i++ ) _ff_poly_big_b[i] = _ff_get_ui(b[i]);
    zn_array_mulmid_precomp1_execute (c, _ff_poly_big_b, *mpctx);
    for ( i = 0 ; i <= (*mpctx)->n1-(*mpctx)->n2 ; i++ ) _ff_set_ui(c[i],c[i]);
}

void ff_poly_middle_product_clear1 (void *ctx)
{
    zn_array_mulmid_precomp1_t *mpctx = ctx;
    zn_array_mulmid_precomp1_clear (*mpctx);
    free (ctx);
}

// computes c[i] = sum_{j+k=d_b+i} a[j]*b[k] for i in [0,d_a-d_b]
void ff_poly_middle_product_big (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b)
{
    register int i;

    assert (d_a >= d_b);
    ff_poly_zn_mod_setup(_ff_p, d_a);
    // this conversion is wasteful
    for ( i = 0 ; i <= d_a ; i++ ) _ff_poly_big_a[i] = _ff_get_ui(a[i]);
    for ( i = 0 ; i <= d_b ; i++ ) _ff_poly_big_b[i] = _ff_get_ui(b[i]);
    zn_array_mulmid (c, _ff_poly_big_a, d_a+1, _ff_poly_big_b, d_b+1, _ff_poly_zn_mod_ctx);
    for ( i = 0 ; i <= d_a-d_b ; i++ ) _ff_set_ui(c[i],c[i]);
}

void ff_poly_mult_big (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b)
{
    register int i;

    ff_poly_zn_mod_setup(_ff_p, _max(d_a,d_b));
    // this conversion is wasteful
    for ( i = 0 ; i <= d_a ; i++ ) _ff_poly_big_a[i] = _ff_get_ui(a[i]);
    for ( i = 0 ; i <= d_b ; i++ ) _ff_poly_big_b[i] = _ff_get_ui(b[i]);
    if ( d_a < d_b ) {
        zn_array_mul (c, _ff_poly_big_b, d_b+1, _ff_poly_big_a, d_a+1, _ff_poly_zn_mod_ctx);
    } else {
        zn_array_mul (c, _ff_poly_big_a, d_a+1, _ff_poly_big_b, d_b+1, _ff_poly_zn_mod_ctx);
    }
    for ( i = 0 ; i <= d_a+d_b ; i++ ) _ff_set_ui(c[i],c[i]);
}

// no conversion out of montgomery rep
void ff_poly_mult_bigx (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b)
{
    ff_poly_zn_mod_setup(_ff_p,0);
    if ( d_a < d_b ) {
        zn_array_mul (c, b, d_b+1, a, d_a+1, _ff_poly_zn_mod_ctx);
    } else {
        zn_array_mul (c, a, d_a+1, b, d_b+1, _ff_poly_zn_mod_ctx);
    }
}

// no conversion out of montgomery rep
void ff_poly_square_bigx (ff_t h[], ff_t f[], int d)
{
    ff_poly_zn_mod_setup(_ff_p, 0);
    zn_array_mul (h, f, d+1, f, d+1, _ff_poly_zn_mod_ctx);
}

void ff_poly_print_raw(ff_t f[], int d)
{
    register int i;
    
    if ( d == 0 ) { printf ("%ld;\n", f[0]); return; }
    if ( d == 1 ) { printf ("%ld*x + %ld;\n", f[1], f[0]); return; }
    printf ("%ld*x^%d", f[d], d);
    for ( i = d-1 ; i > 1 ; i-- ) if ( f[i] ) printf (" +%ld*x^%d", f[i], i);
    printf (" + %ld*x + %ld;\n", f[1], f[0]);
}
