#include <assert.h>
#include <stdio.h>
#include "mm.h"
#include "ffpoly.h"
#include "mmpoly.h"


// computes o(x)=f(x)*g(x) where deg(f)=deg(g)=d+1 with f and g monic (coeff of x^d implicit, not examined in input and not set in output)
static inline void mm_poly_mul_m_d_d (mm_t o[], mm_t f[], mm_t g[], int d, mm_t R, mm_t p, mm_t pinv)
{
    o[2*d-1] = mm_add(f[d-1],g[d-1],p);
    for ( register int i = d-1 ; i > 0 ; i-- ) o[d+i-1] = mm_redc (mm2_conv(f+i, g+i, d-i) + (mm2_t)R*(f[i-1]+g[i-1]), p, pinv);
    for ( register int i = d-1 ; i >= 0 ; i-- ) o[i] = mm_conv(f, g, i+1, p, pinv);
}
 
// computes o(x)=f(x)*g(x) where deg(f)=d+1 and deg(g)=d with f and g monic (coeff of x^d implicit, not examined in input and not set in output)
static inline void mm_poly_mul_m_d_do (mm_t o[], mm_t f[], mm_t g[], int d, mm_t R, mm_t p, mm_t pinv)
{
    o[2*d-2] = mm_add(f[d-1],g[d-2],p);
    for ( register int i = d-1 ; i > 1 ; i-- ) o[d+i-2] = mm_redc (mm2_conv(f+i, g+i-1, d-i) + (mm2_t)R*(f[i-1]+g[i-2]), p, pinv);
    o[d-1] = mm_redc (mm2_conv(f+1,g,d-1)+(mm2_t)R*f[0],p,pinv);
    for ( register int i = d-2 ; i >= 0 ; i-- ) o[i] = mm_conv(f, g, i+1, p, pinv);
}


// computes f(x) = \prod_{i=0}^{n-1} (x-r[i]) without setting the (implicit) leading monic coefficient
void mm_poly_from_roots (mm_t o[], mm_t r[], int d, mm_t R, mm_t p, mm_t pinv)
{
    switch (d) {
    case 0: return;
    case 1: mm_poly_from_roots_1 (o, r, R, p, pinv); return;
    case 2: mm_poly_from_roots_2 (o, r, R, p, pinv); return;
    case 3: mm_poly_from_roots_3 (o, r, R, p, pinv); return;
    case 4: mm_poly_from_roots_4 (o, r, R, p, pinv); return;
    case 5: mm_poly_from_roots_5 (o, r, R, p, pinv); return;
    case 6: mm_poly_from_roots_6 (o, r, R, p, pinv); return;
    case 7: mm_poly_from_roots_7 (o, r, R, p, pinv); return;
    case 8: mm_poly_from_roots_8 (o, r, R, p, pinv); return;
    case 9: mm_poly_from_roots_9 (o, r, R, p, pinv); return;
    case 10: mm_poly_from_roots_10 (o, r, R, p, pinv); return;
    case 11: mm_poly_from_roots_11 (o, r, R, p, pinv); return;
    case 12: mm_poly_from_roots_12 (o, r, R, p, pinv); return;
    case 13: mm_poly_from_roots_13 (o, r, R, p, pinv); return;
    case 14: mm_poly_from_roots_14 (o, r, R, p, pinv); return;
    case 15: mm_poly_from_roots_15 (o, r, R, p, pinv); return;
    case 16: mm_poly_from_roots_16 (o, r, R, p, pinv); return;
    case 17: mm_poly_from_roots_17 (o, r, R, p, pinv); return;
    case 18: mm_poly_from_roots_18 (o, r, R, p, pinv); return;
    case 19: mm_poly_from_roots_19 (o, r, R, p, pinv); return;
    case 20: mm_poly_from_roots_20 (o, r, R, p, pinv); return;
    case 21: mm_poly_from_roots_21 (o, r, R, p, pinv); return;
    case 22: mm_poly_from_roots_22 (o, r, R, p, pinv); return;
    case 23: mm_poly_from_roots_23 (o, r, R, p, pinv); return;
    case 24: mm_poly_from_roots_24 (o, r, R, p, pinv); return;
    case 25: mm_poly_from_roots_25 (o, r, R, p, pinv); return;
    case 26: mm_poly_from_roots_26 (o, r, R, p, pinv); return;
    case 27: mm_poly_from_roots_27 (o, r, R, p, pinv); return;
    case 28: mm_poly_from_roots_28 (o, r, R, p, pinv); return;
    case 29: mm_poly_from_roots_29 (o, r, R, p, pinv); return;
    case 30: mm_poly_from_roots_30 (o, r, R, p, pinv); return;
    case 31: mm_poly_from_roots_31 (o, r, R, p, pinv); return;
    case 32: mm_poly_from_roots_32 (o, r, R, p, pinv); return;
    }
    int m = (d+1)/2;
    int n = d/2;
    mm_t f[m], g[n];
    
    mm_poly_from_roots (f,r,m,R,p,pinv); mm_poly_from_roots(g,r+m,n,R,p,pinv);
    if ( m == n ) mm_poly_mul_m_d_d (o,f,g,m,R,p,pinv);
    else mm_poly_mul_m_d_do (o,f,g,m,R,p,pinv);
}


// given f(x) and g(x) with deg(f) >= deg(g) > 0 replaces f(x) with a linear combination of f and g with degree less than deg(f)
// if deg(g) < deg(f) the result will have degree at most deg(f)-2 -- typical case is deg(f)=deg(g)+1 with output of degree deg(g)-1.
static inline int mm_poly_red (mm_t f[], int df, mm_t g[], int dg, mm_t p, mm_t pinv)
{
    if ( df == dg ) {
        for ( register int i = df-1 ; i >= 0 ; i-- ) f[i] = mm_det_2r (g[dg],f[df],g[i],f[i],p,pinv);
        return mm_poly_deg(f,df-1);
    } else {
        mm_t s, r0, r1;
        int j = df-dg-1;
        s = mm_mul(g[dg],g[dg],p,pinv);
        r0 = mm_det_2r (f[df],f[df-1],g[dg],g[dg-1],p,pinv);
        r1 = mm_neg(mm_mul(f[df],g[dg],p,pinv),p);
        for ( register int i = df-2 ; i > j ; i-- ) f[i] = mm_dot_3r (s,r0,r1,f[i],g[i-j],g[i-j-1],p,pinv);
        f[j] = mm_dot_2r (s,r0,f[j],g[0],p,pinv);
        for ( register int i=j-1 ; i >= 0 ; i-- ) f[i] = mm_mul (s,f[i],p,pinv);
        return mm_poly_deg(f,df-2);
    }
}

// deg(f)=d+1, deg(g)=d > 0, deg(f) < d on return, zero padded to degree d
static inline void mm_poly_red_normal (mm_t f[], mm_t g[], int d, mm_t p, mm_t pinv)
{
    register mm_t s, r0, r1;
    s = mm_mul(g[d],g[d],p,pinv);
    r0 = mm_det_2r (f[d+1],f[d],g[d],g[d-1],p,pinv);
    r1 = mm_neg(mm_mul(f[d+1],g[d],p,pinv),p);
    f[d] = 0;
    for ( register int i = d-1 ; i ; i-- ) f[i] = mm_dot_3r (s,r0,r1,f[i],g[i],g[i-1],p,pinv);
    f[0] = mm_dot_2r (s,r0,f[0],g[0],p,pinv);
}

// given deg(f)=d+1 and deg(g)=d > 1, replaces f and g with lower degree generators of the ideal (f,g)
// returns m=max(deg(f),deg(g)) and ensures that f and g are both zero-padded to this degree (unless m=0 in which case f and g should be ignored)
// in the typical case one of f and g is reduced to zero while the other is their GCD
static inline int mm_poly_gcd_red_normal (mm_t f[], mm_t g[], int d, mm_t p, mm_t pinv)
{
    do {
        mm_poly_red_normal (f,g,d--,p,pinv);  if ( !f[d] ) return d+1;
        mm_poly_red_normal (g,f,d--,p,pinv);  if ( !g[d] ) return d+1;
    } while ( d > 1 );
    if ( !d ) return 0;
    // if f mod g is zero then g holds linear gcd and o.w. gcd is 1; in either case, setting f[1]=f[0]=f mod g does what we want
    f[1] = f[0] = mm_dot_2r (mm_mul(g[1],g[1],p,pinv),mm_det_2r(f[2],f[1],g[1],g[0],p,pinv),f[0],g[0],p,pinv);
    return f[0] ? 0 : 1;
}

static inline int mm_poly_gcd_tiny (mm_t o[], mm_t f[], int df, mm_t g[], int dg, mm_t p, mm_t pinv)
{
    switch (df) {
    case 0: o[0] = 1; return 0;
    case 1: switch (dg) { case 0: o[0]=1; return 0; case 1: return mm_poly_gcd_1_1 (o,f,g,p,pinv); case 2: return mm_poly_gcd_2_1 (o,g,f,p,pinv); case 3: return mm_poly_gcd_3_1 (o,g,f,p,pinv); case 4: return mm_poly_gcd_4_1 (o,g,f,p,pinv);
                          case 5: return mm_poly_gcd_5_1 (o,g,f,p,pinv); case 6: return mm_poly_gcd_6_1 (o,g,f,p,pinv); }
    case 2: switch (dg) { case 0: o[0]=1; return 0; case 1: return mm_poly_gcd_2_1 (o,f,g,p,pinv); case 2: return mm_poly_gcd_2_2 (o,f,g,p,pinv); case 3: return mm_poly_gcd_3_2 (o,g,f,p,pinv); case 4: return mm_poly_gcd_4_2 (o,g,f,p,pinv);
                          case 5: return mm_poly_gcd_5_2 (o,g,f,p,pinv); case 6: return mm_poly_gcd_6_2 (o,g,f,p,pinv); }
    case 3: switch (dg) { case 0: o[0]=1; return 0; case 1: return mm_poly_gcd_3_1 (o,f,g,p,pinv); case 2: return mm_poly_gcd_3_2 (o,f,g,p,pinv); case 3: return mm_poly_gcd_3_3 (o,f,g,p,pinv); case 4: return mm_poly_gcd_4_3 (o,g,f,p,pinv);
                          case 5: return mm_poly_gcd_5_3 (o,g,f,p,pinv); case 6: return mm_poly_gcd_6_3 (o,g,f,p,pinv); }
    case 4: switch (dg) { case 0: o[0]=1; return 0; case 1: return mm_poly_gcd_4_1 (o,f,g,p,pinv); case 2: return mm_poly_gcd_4_2 (o,f,g,p,pinv); case 3: return mm_poly_gcd_4_3 (o,f,g,p,pinv); case 4: return mm_poly_gcd_4_4 (o,f,g,p,pinv);
                          case 5: return mm_poly_gcd_5_4 (o,g,f,p,pinv); case 6: return mm_poly_gcd_6_4 (o,g,f,p,pinv); }
    case 5: switch (dg) { case 0: o[0]=1; return 0; case 1: return mm_poly_gcd_5_1 (o,f,g,p,pinv); case 2: return mm_poly_gcd_5_2 (o,f,g,p,pinv); case 3: return mm_poly_gcd_5_3 (o,f,g,p,pinv); case 4: return mm_poly_gcd_5_4 (o,f,g,p,pinv);
                          case 5: return mm_poly_gcd_5_5 (o,f,g,p,pinv); case 6: return mm_poly_gcd_6_5 (o,g,f,p,pinv); }
    case 6: switch (dg) { case 0: o[0]=1; return 0; case 1: return mm_poly_gcd_6_1 (o,f,g,p,pinv); case 2: return mm_poly_gcd_6_2 (o,f,g,p,pinv); case 3: return mm_poly_gcd_6_3 (o,f,g,p,pinv); case 4: return mm_poly_gcd_6_4 (o,f,g,p,pinv);
                          case 5: return mm_poly_gcd_6_5 (o,f,g,p,pinv); case 6: return mm_poly_gcd_6_6 (o,f,g,p,pinv); }
    }
    return -1;
}


// given f,g with deg(f) >= deg(g), destructively computes f = gcd(f,g) destroying g in the process.
// returns deg(f) (returns -1 if df and dg are both 0)
static inline int mm_poly_gcd_red (mm_t f[], int df, mm_t g[], int dg, mm_t p, mm_t pinv)
{
    while ( dg > 0 ) {
        while ( df >= dg ) df = mm_poly_red (f,df,g,dg,p,pinv);
        if ( df <= 0 ) break;
        while ( dg >= df ) dg = mm_poly_red (g,dg,f,df,p,pinv);
    }
    if ( !df || !dg ) { f[0]=1; return 0; }
    if ( df < 0 ) { memcpy(f,g,(dg+1)*sizeof(g[0])); df = dg; }
    return df;
}


// computes o=gcd(f,g), one of f and g must be nonzero -- optimized for the case df=dg+1
int mm_poly_gcd (mm_t o[], mm_t f[], int df, mm_t g[], int dg, mm_t R, mm_t R2, mm_t R3, mm_t p, mm_t pinv)
{
    if ( df < 0 ) { memmove(o,g,(dg+1)*sizeof(*g)); return dg; }
    if ( dg < 0 ) { memmove(o,f,(df+1)*sizeof(*f)); return df; }
    if ( df <= 6 && dg <= 6 ) return mm_poly_gcd_tiny (o,f,df,g,dg,p,pinv);
    
    if ( !df || !dg ) { o[0] = 1; return 0; }

    mm_t S[df+1], T[dg+1];
    mm_t *s=S, *t=T;
    memcpy(s,f,(df+1)*sizeof(*f)); memcpy (t,g,(dg+1)*sizeof(*g));
    if ( df < dg ) { SWAP(s,t); SWAP(df,dg); }
    
    if ( df > 15 && df > dg+1 ) {
        mm_poly_make_monic(t,t,dg,R,R2,R3,p,pinv);
        for ( int i = 0 ; i < dg ; i ++ ) t[i] = mm_neg(t[i],p);
        mm_poly_mod_m (s,s,df,t,dg,R,p,pinv);
        t[dg] = mm_neg(t[dg],p);
        return mm_poly_gcd(o,t,dg,s,mm_poly_deg(s,dg-1),R,R2,R3,p,pinv);
    }
    if ( df == dg+1 ) {
        df = mm_poly_gcd_red_normal(s,t,dg,p,pinv);
        if ( !df ) { o[0]=1; return 0; }
        if ( ! s[df] ) SWAP(s,t);
        dg = mm_poly_deg(t,df-2);
        if ( dg < 0 ) { memmove (o,s,(df+1)*sizeof(s[0]));  return df; }
        if ( !dg ) { o[0]=1; return 0; }
    }
    df = mm_poly_gcd_red (s, df, t, dg, p, pinv);
    memmove (o,s,(df+1)*sizeof(s[0]));
    return df;
}


// 2M + 1A +1R
static inline void mm_poly_res_1_1 (mm_t o[2],  mm_t f[2], mm_t g[2],  mm_t R, mm_t p, mm_t pinv)
    { o[0] = mm_det_2r(f[1],f[0],g[1],g[0],p,pinv); o[1] = R; }

// 5M + 2A + 3R
static inline void mm_poly_res_2_1 (mm_t o[2],  mm_t f[3], mm_t g[2], mm_t R, mm_t p, mm_t pinv)
    { o[0] = mm_dot_2r (mm_mul(g[1],g[1],p,pinv),mm_det_2r(f[2],f[1],g[1],g[0],p,pinv),f[0],g[0],p,pinv); o[1] = R; }

// 10M + 4A + 6R
static inline void mm_poly_res_2_2 (mm_t o[2],  mm_t f[3], mm_t g[3], mm_t R, mm_t p, mm_t pinv) {
    mm_t t[2];
    t[1] = mm_det_2r (g[2],f[2],g[1],f[1],p,pinv);
    t[0] = mm_det_2r (g[2],f[2],g[0],f[0],p,pinv);
    if ( t[1] ) { mm_poly_res_2_1 (o,g,t,R,p,pinv);  o[1] = mm_mul(o[1],g[2],p,pinv); return; }
    o[1] = R;
    o[0] = t[0] ? mm_sqr(t[0],p,pinv) : 0;
}

// 8M + 4A + 5R
static inline void mm_poly_res_3_1 (mm_t o[2], mm_t f[4], mm_t g[2], mm_t R, mm_t p, mm_t pinv)
{
    mm_t t[2];
    register mm_t  s = mm_mul(g[1],g[1],p,pinv);
    register mm_t r0 = mm_det_2r (f[3],f[2],g[1],g[0],p,pinv);
    t[1] = mm_dot_2r (s,r0,f[1],g[0],p,pinv);
    t[0] = mm_mul (s,f[0],p,pinv);
    if ( t[1] ) { mm_poly_res_1_1 (o,g,t,R,p,pinv); o[0] = mm_neg(o[0],p); return; }
    o[0] = t[0] ? mm_neg(mm_mul(t[0],g[1],p,pinv),p) : 0;
    o[1] = R;
}

// 15M + 7A + 9R
static inline void mm_poly_res_3_2 (mm_t o[2], mm_t f[4], mm_t g[3], mm_t R, mm_t p, mm_t pinv)
{
    mm_t t[2];
    register mm_t  s = mm_mul(g[2],g[2],p,pinv);
    register mm_t r0 = mm_det_2r (f[3],f[2],g[2],g[1],p,pinv);
    register mm_t r1 = mm_neg(mm_mul(f[3],g[2],p,pinv),p);
    t[1] = mm_dot_3r (s,r0,r1,f[1],g[1],g[0],p,pinv);
    t[0] = mm_dot_2r (s,r0,f[0],g[0],p,pinv);
    if ( t[1] ) { mm_poly_res_2_1 (o,g,t,R,p,pinv); o[1] = mm_sqr(g[2],p,pinv); return; }
    if ( t[0] ) { o[0] = mm_sqr(t[0],p,pinv); o[1] = g[2]; return; }
    o[0] = 0; o[1] = R;
}

// 23M + 11A + 14R
static inline void mm_poly_res_3_3 (mm_t o[2], mm_t f[4], mm_t g[4], mm_t R, mm_t p, mm_t pinv)
{
    mm_t t[3];
    t[2] = mm_det_2r (g[3],f[3],g[2],f[2],p,pinv);
    t[1] = mm_det_2r (g[3],f[3],g[1],f[1],p,pinv);
    t[0] = mm_det_2r (g[3],f[3],g[0],f[0],p,pinv);
    if ( t[2] ) { mm_poly_res_3_2 (o,g,t,R,p,pinv); o[0] = mm_neg(o[0],p); o[1] = mm_mul(o[1],mm_sqr(g[3],p,pinv),p,pinv); return; }
    if ( t[1] ) { mm_poly_res_3_1 (o,g,t,R,p,pinv); o[0] = mm_neg(o[0],p); o[1] = mm_mul(o[1],g[3],p,pinv); return; }
    o[0] = t[0] ? mm_neg(mm_mul(t[0],mm_sqr(t[0],p,pinv),p,pinv),p) : 0;
    o[1] = R;
}

// 12M + 4A + 8R
static inline void mm_poly_res_4_1 (mm_t o[2], mm_t f[5], mm_t g[2], mm_t R, mm_t p, mm_t pinv)
{
    mm_t t[3];
    register mm_t s  = mm_mul(g[1],g[1],p,pinv);
    register mm_t r0 = mm_det_2r (f[4],f[3],g[1],g[0],p,pinv);
    t[2] = mm_dot_2r (s,r0,f[2],g[0],p,pinv);
    t[1] = mm_mul (s,f[1],p,pinv);
    t[0] = mm_mul (s,f[0],p,pinv);
    if ( t[2] ) { mm_poly_res_2_1 (o,t,g,R,p,pinv); return; }
    if ( t[1] ) { mm_poly_res_1_1 (o,g,t,R,p,pinv); o[0] = mm_mul(o[0],g[1],p,pinv); return; }
    o[0] = t[0] ? mm_mul(t[0],mm_sqr(g[1],p,pinv),p,pinv) : 0;
    o[1] = R;
}

// 22M + 9A + 14R
static inline void mm_poly_res_4_2 (mm_t o[2], mm_t f[5], mm_t g[3], mm_t R, mm_t p, mm_t pinv)
{
    mm_t t[3];
    register mm_t  s = mm_mul(g[2],g[2],p,pinv);
    register mm_t r0 = mm_det_2r (f[4],f[3],g[2],g[1],p,pinv);
    register mm_t r1 = mm_neg(mm_mul(f[4],g[2],p,pinv),p);
    t[2] = mm_dot_3r (s,r0,r1,f[2],g[1],g[0],p,pinv);
    t[1] = mm_dot_2r (s,r0,f[1],g[0],p,pinv);
    t[0] = mm_mul (s,f[0],p,pinv);
    if ( t[2] ) { mm_poly_res_2_2 (o,g,t,R,p,pinv); o[1] = mm_mul (o[1],mm_sqr(g[2],p,pinv),p,pinv); return; }
    if ( t[1] ) { mm_poly_res_2_1 (o,g,t,R,p,pinv); o[1] = mm_mul (o[1],g[2],p,pinv); return; }
    o[0] = t[0] ? mm_sqr(t[0],p,pinv) : 0;
    o[1] = R;
}

// 30M + 14A + 18R
static inline void mm_poly_res_4_3 (mm_t o[2], mm_t f[5], mm_t g[4], mm_t R, mm_t p, mm_t pinv)
{
    mm_t t[3];
    register mm_t  s = mm_mul(g[3],g[3],p,pinv);
    register mm_t r0 = mm_det_2r (f[4],f[3],g[3],g[2],p,pinv);
    register mm_t r1 = mm_neg(mm_mul(f[4],g[3],p,pinv),p);
    t[2] = mm_dot_3r (s,r0,r1,f[2],g[2],g[1],p,pinv);
    t[1] = mm_dot_3r (s,r0,r1,f[1],g[1],g[0],p,pinv);
    t[0] = mm_dot_2r (s,r0,f[0],g[0],p,pinv);
    if ( t[2] ) { mm_poly_res_3_2 (o,g,t,R,p,pinv); o[1] = mm_mul (o[1],mm_sqr(mm_sqr(g[3],p,pinv),p,pinv),p,pinv); return; }
    if ( t[1] ) { mm_poly_res_3_1 (o,g,t,R,p,pinv); o[1] = mm_mul (mm_mul(o[1],g[3],p,pinv),mm_sqr(g[3],p,pinv),p,pinv); return; }
    if ( t[0] ) { o[0] = mm_mul(t[0],mm_sqr(t[0],p,pinv),p,pinv); o[1] = mm_sqr(g[3],p,pinv); return; }
    o[0] = 0; o[1] = R;
}

// 41M + 22A + 29R
static inline void mm_poly_res_4_4 (mm_t o[2], mm_t f[5], mm_t g[5], mm_t R, mm_t p, mm_t pinv)
{
    mm_t t[4];
    t[3] = mm_det_2r (g[4],f[4],g[3],f[3],p,pinv);
    t[2] = mm_det_2r (g[4],f[4],g[2],f[2],p,pinv);
    t[1] = mm_det_2r (g[4],f[4],g[1],f[1],p,pinv);
    t[0] = mm_det_2r (g[4],f[4],g[0],f[0],p,pinv);
    if ( t[3] ) { mm_poly_res_4_3 (o,g,t,R,p,pinv); o[1] = mm_mul (mm_mul(o[1],g[4],p,pinv),mm_sqr(g[4],p,pinv),p,pinv); return; }
    if ( t[2] ) { mm_poly_res_4_2 (o,g,t,R,p,pinv); o[1] = mm_mul (o[1],mm_sqr(g[4],p,pinv),p,pinv); return; }
    if ( t[1] ) { mm_poly_res_4_1 (o,g,t,R,p,pinv); o[1] = mm_mul (o[1],g[4],p,pinv); return; }
    o[0] = t[0] ? mm_sqr(mm_sqr(t[0],p,pinv),p,pinv) : 0;
    o[1] = R;
}

static inline void mm_poly_res_tiny (mm_t o[2], mm_t f[], int df, mm_t g[], int dg, mm_t R, mm_t p, mm_t pinv)
{
    switch (df) {
    case 0: o[1] = R; switch (dg) { case 0: o[0] = R; return;
                                    case 1: o[0] = f[0]; return; case 2: o[0] = mm_sqr(f[0],p,pinv); return; case 3: o[0] = mm_mul(f[0],mm_sqr(f[0],p,pinv),p,pinv); return; case 4: o[0] = mm_sqr(mm_sqr(f[0],p,pinv),p,pinv); return; }
    case 1: switch (dg) { case 0: o[0] = g[0]; o[1] = R; return;
                          case 1: mm_poly_res_1_1 (o,f,g,R,p,pinv); return; case 2: mm_poly_res_2_1 (o,g,f,R,p,pinv); return; case 3: mm_poly_res_3_1 (o,g,f,R,p,pinv); o[0] = mm_neg(o[0],p); return; case 4: mm_poly_res_4_1 (o,g,f,R,p,pinv); return; }
    case 2: switch (dg) { case 0: o[0] = mm_sqr(g[0],p,pinv); o[1] = R; return;
                          case 1: mm_poly_res_2_1 (o,f,g,R,p,pinv); return; case 2: mm_poly_res_2_2 (o,f,g,R,p,pinv); return; case 3: mm_poly_res_3_2 (o,g,f,R,p,pinv); return; case 4: mm_poly_res_4_2 (o,g,f,R,p,pinv); return; }
    case 3: switch (dg) { case 0: o[0] = mm_mul(g[0],mm_sqr(g[0],p,pinv),p,pinv); o[1] = R; return;
                          case 1: mm_poly_res_3_1 (o,f,g,R,p,pinv); return; case 2: mm_poly_res_3_2 (o,f,g,R,p,pinv); return; case 3: mm_poly_res_3_3 (o,f,g,R,p,pinv); return; case 4: mm_poly_res_4_3 (o,g,f,R,p,pinv); return; }
    case 4: switch (dg) { case 0: o[0] = mm_sqr(mm_sqr(g[0],p,pinv),p,pinv); o[1] = R; return;
                          case 1: mm_poly_res_4_1 (o,f,g,R,p,pinv); return; case 2: mm_poly_res_4_2 (o,f,g,R,p,pinv); return; case 3: mm_poly_res_4_3 (o,f,g,R,p,pinv); return; case 4: mm_poly_res_4_4 (o,f,g,R,p,pinv); return; }
    }
}


// given deg(f)=d+1 and deg(g)=d > 1, destructively computes res(f,g) = o[0]/o[1] provided that the degree sequence is normal down to the gcd
// retuerns zero if this is not the case, in which case the caller will need to revert to the general method (note that f and g will be trashed)
static inline int mm_poly_res_normal (mm_t o[], mm_t f[], mm_t g[], int d, mm_t R, mm_t p, mm_t pinv)
{
    mm_t t, u;

    // in the normal degree sequence case either the resultant is zero or the degree of f or g drops by 2 in each call to mm_poly_red_normal
    // in particular, the degrees of f and g always have opposite parity and the degree drop of 2 always has opposite parity from the smaller degree
    // this means that we don't have to keep track of signs (there is never a change) and after each call to mm_poly_red_normal we want to multiply
    // the denominator of the resultant by lc(g)^(2deg(g)-2) or lc(f)^(2deg(f)-2) (depending on which of g or f has smaller degree).

    t = 0;
    do {
        mm_poly_red_normal (f,g,d--,p,pinv);  if ( !f[d] ) { for ( d-- ; d >=0 && !f[d] ; d-- );  o[0] = 0; o[1] = R; return d < 0 ? 1 : 0; }
        if ( t ) { t = mm_mul (t, g[d+1],p,pinv); u = mm_mul (u,t,p,pinv); } else { u = t = g[d+1]; }
        mm_poly_red_normal (g,f,d--,p,pinv);  if ( !g[d] ) { for ( d-- ; d >=0 && !g[d] ; d-- );  o[0] = 0; o[1] = R; return d < 0 ? 1 : 0; }
        if ( d ) { t = mm_mul (t, f[d+1],p,pinv); u = mm_mul (u,t,p,pinv); }
    } while ( d > 1 );
    if ( d ) { mm_poly_res_2_1 (o, f, g, R, p, pinv); } else { o[0] = g[0]; }
    o[1] = mm_sqr (u,p,pinv); // note that o[1] = R in either case above so this assignment is equivalent to a multiplication
    return 1;
}

// computes the resultant res(f,g) as o[0]/o[1] by computing r = lc(g)^a*f - h(x)*g with a=1 or 2 and deg(r) <= deg(f)-a and applying the identity
//     res(f,g) = (-1)^(deg(f)*deg(g)) * lc(g)^(deg(f)-deg(r)-a*deg(g)) * res(g,r)
void _mm_poly_res (mm_t o[2], mm_t f[], int df, mm_t g[], int dg, mm_t R, mm_t p, mm_t pinv)
{
    int d,e,s;

    if ( df < 0 || dg < 0 ) { o[0] = 0; o[1] = R; return; }
    if ( df <= 4 && dg <= 4 ) { mm_poly_res_tiny (o,f,df,g,dg,R,p,pinv); return; }
    if ( !df ) { o[0] = mm_exp_ui(f[0],dg,R,p,pinv); o[1] = R; return; }
    if ( !dg ) { o[0] = mm_exp_ui(g[0],df,R,p,pinv); o[1] = R; return; }
    
    mm_t F[df+1], G[dg+1], *sf, *sg;
    memcpy (F, f, (df+1)*sizeof(mm_t)); sf = f; f = F;
    memcpy (G, g, (dg+1)*sizeof(mm_t)); sg = g; g = G;
    if ( df == dg+1 ) { if ( mm_poly_res_normal (o,f,g,dg,R,p,pinv) ) return; else { memcpy (f,sf,(df+1)*sizeof(mm_t)); memcpy (g,sg,(dg+1)*sizeof(mm_t)); } }
    if ( dg == df+1 ) { if ( mm_poly_res_normal (o,g,f,df,R,p,pinv) ) return; else { memcpy (f,sf,(df+1)*sizeof(mm_t)); memcpy (g,sg,(dg+1)*sizeof(mm_t)); } }
    
    o[0] = o[1] = R;
    s = 1;
    do { // df and dg are both positive here
        if ( df < dg ) { SWAP(f,g); SWAP(df,dg); s = ((df&dg)&1) ? -s : s; } // R(f,g) = (-1)^(df*dg)R(g,f)
        d = mm_poly_red (f, df, g, dg, p, pinv); // reduces degree of f by a=1,2 and multiplies by lc(g)^a in the process
        if ( d < 0 ) { o[0] = 0; o[1] = R; return; }
        e = df - d - (df == dg ? dg : 2*dg);     // a=1 if df=dg and a=2 if df > dg
        // R(f,g) = (-1)^((df+d)*dg) * lc(g)^(df-d+a*dg)*R(lc(g)^a*f mod g, g)
        if ( e < 0 ) o[1] = mm_mul(o[1],mm_exp_ui(g[dg],-e,R,p,pinv),p,pinv);
        if ( e > 0 ) o[0] = mm_mul(o[0],mm_exp_ui(g[dg],e,R,p,pinv),p,pinv);
        s = (((df+d)&dg)&1) ? -s : s;
        df = d;
    } while ( df );
    o[0] = mm_mul(o[0],mm_exp_ui(f[0],dg,R,p,pinv),p,pinv);
    if ( s < 0 ) o[0] = mm_neg(o[0],p);
}


// given deg(f)=d+1 and deg(g)=d > 1, destructively computes res(f,g) = o[0]/o[1] provided that the degree sequence is normal down to the gcd
// retuerns zero if this is not the case, in which case the caller will need to revert to the general method (note that f and g will be trashed)
static inline int mm_poly_res_legendre_normal (mm_t f[], mm_t g[], int d, mm_t R, mm_t p, mm_t pinv)
{
    // in the normal degree sequence case either the resultant is zero or the degree of f or g drops by 2 in each call to mm_poly_red_normal
    // in particular, the degrees of f and g always have opposite parity and the degree drop of 2 always has opposite parity from the smaller degree
    // this means that we don't have to keep track of signs (there is never a change) and after each call to mm_poly_red_normal we want to multiply
    // the denominator of the resultant by a square, which we can safely ignore

    do {
        mm_poly_red_normal (f,g,d--,p,pinv);  if ( !f[d] ) { for ( d-- ; d >=0 && !f[d] ; d-- );  return d < 0 ? 0 : -2; }
        mm_poly_red_normal (g,f,d--,p,pinv);  if ( !g[d] ) { for ( d-- ; d >=0 && !g[d] ; d-- );  return d < 0 ? 0 : -2; }
    } while ( d > 1 );
    mm_t o[2];
    if ( d ) { mm_poly_res_2_1 (o, f, g, R, p, pinv); } else { o[0] = g[0]; }
    return mm_legendre (o[0],R,p,pinv);
}

// computes the resultant res(f,g) as o[0]/o[1] by computing r = lc(g)^a*f - h(x)*g with a=1 or 2 and deg(r) <= deg(f)-a and applying the identity
//     res(f,g) = (-1)^(deg(f)*deg(g)) * lc(g)^(deg(f)-deg(r)-a*deg(g)) * res(g,r)
int mm_poly_res_legendre (mm_t f[], int df, mm_t g[], int dg, mm_t R, mm_t p, mm_t pinv)
{
    mm_t t,o[2];
    int d,e,s;

    if ( df < 0 || dg < 0 ) return 0;
    if ( !df ) return ( (dg&1) ? 1 : mm_legendre (f[0],R,p,pinv) );
    if ( !dg ) return ( (df&1) ? 1 : mm_legendre (g[0],R,p,pinv) );
    if ( df == 1 && dg == 1 ) { mm_poly_res_1_1 (o,f,g,R,p,pinv); return mm_legendre (o[0],R,p,pinv); }
    if ( df == 2 && dg == 1 ) { mm_poly_res_2_1 (o,f,g,R,p,pinv); return mm_legendre (o[0],R,p,pinv); }
    if ( df == 1 && dg == 2 ) { mm_poly_res_2_1 (o,g,f,R,p,pinv); return mm_legendre (o[0],R,p,pinv); }
    
    mm_t F[df+1], G[dg+1], *sf, *sg;
    memcpy (F, f, (df+1)*sizeof(mm_t)); sf = f; f = F;
    memcpy (G, g, (dg+1)*sizeof(mm_t)); sg = g; g = G;
    if ( df == dg+1 ) { d = mm_poly_res_legendre_normal (f,g,dg,R,p,pinv); if ( d >= -1 ) return d; else { memcpy (f,sf,(df+1)*sizeof(mm_t)); memcpy (g,sg,(dg+1)*sizeof(mm_t)); } }
    if ( dg == df+1 ) { d = mm_poly_res_legendre_normal (g,f,df,R,p,pinv); if ( d >= -1 ) return d; else { memcpy (f,sf,(df+1)*sizeof(mm_t)); memcpy (g,sg,(dg+1)*sizeof(mm_t)); } }

    t = R; s = 1;
    do { // df and dg are both positive here
        if ( df < dg ) { SWAP(f,g); SWAP(df,dg); s = ((df&dg)&1) ? -s : s; } // R(f,g) = (-1)^(df*dg)R(g,f)
        d = mm_poly_red (f, df, g, dg, p, pinv); // reduces degree of f by a=1,2 and multiplies by lc(g)^a in the process
        if ( d < 0 ) return 0;
        e = df + d + (df == dg ? dg : 2*dg);     // a=1 if df=dg and a=2 if df > dg, flip signs because we only care about the parity of e
        if ( (e&1) ) t = mm_mul(t,g[dg],p,pinv);
        s = (((df+d)&dg)&1) ? -s : s;
        df = d;
    } while ( df );
    if ( dg&1 ) t = mm_mul(t,f[0],p,pinv);
    d = mm_legendre(t,R,p,pinv);
    return ( s < 0 && (p&3) == 3 ) ? -d : d;    // account for sign
}

// destructively computes determinant of n x n matrix n without using inversions
// numerator of the determinant is stored in m[0] and the denominator is returned
mm_t _mm_det (mm_t m[], int n, mm_t R, mm_t p, mm_t pinv)
{
    register int i,k;
    register mm_t x,*s,*t,*e,*me=m+n*n;

    switch (n) {
    case 1: return R;
    case 2: m[0] = mm_det_2 (m,p,pinv); return R;
    case 3: m[0] = mm_det_3 (m,p,pinv); return R;
    case 4: m[0] = mm_det_4 (m,p,pinv); return R;
    }
    if ( n==1 ) return R;
    x = R; t = m;
    for ( k = 0 ; k < n-2 ; k++, t+=(n+1) ) {  // *t = m[k*n+k] is the kth diagonal entry
        if ( ! *t ) {
            for ( i = k+1 ; i < n ; i++ ) if ( m[i*n+k] ) break;
            if ( i == n ) { m[0] = 0; return R; }
            for ( register int j = k ; j < n ; j++ ) SWAP(m[i*n+j],m[k*n+j]);
            x = mm_neg(x,p);
        }
        register mm_t a = *t;                  // a = m[k*n+k]
        x = mm_mul(x,mm_exp_ui(a,n-k-2,R,p,pinv),p,pinv);
        for ( s = t+n, e = s+n-k ; s < me ; s+=k, e+=n ) { // *s = m[i*n+k], *e = m[i*n+n] with i ranging from k+1 to n-1
            register mm_t b = *s++;
            for ( register mm_t *u=t+1 ; s < e ; s++, u++ ) *s = mm_det_2r (a,*u,b,*s,p,pinv);
        }

    }
    m[0] = mm_det_2r (t[0],t[1],t[n],t[n+1],p,pinv);
    return x;
}

static inline int mm_bdlog_1 (mm_t a[1], mm_t ainv[1], mm_t b, mm_t R, mm_t p, mm_t pinv)
    { return ( b == R ? 0 : 1 ); }

static inline int mm_bdlog_1x (mm_t a[1], mm_t ainv[1], mm_t b, mm_t R, mm_t p, mm_t pinv)
    { return ( b == R ? 0 : ( b == p-R ? 1 : -1 ) ); }

static inline int mm_bdlog_2 (mm_t a[2], mm_t ainv[2], mm_t b, mm_t R, mm_t p, mm_t pinv)
    { return ( b == R ? 0 : ( b == a[0] ? 1 : ( b == a[1] ? 2 : 3 ) ) ); }
    
static inline int mm_bdlog_2x (mm_t a[2], mm_t ainv[2], mm_t b, mm_t R, mm_t p, mm_t pinv)
    { return ( b == R ? 0 : ( b == a[0] ? 1 : ( b == a[1] ? 2 : (b == p-a[0] ? 3 : -1 ) ) ) ); }

static inline int mm_bdlog_3 (mm_t a[3], mm_t ainv[3], mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    if ( b == R ) return 0;
    if ( b == a[0] ) return 1;
    if ( b == a[1] ) return 2;
    if ( b == a[2] ) return 4;
    if ( b == ainv[0] ) return 7;
    if ( b == ainv[1] ) return 6;
    return ( b == p-a[0] ? 5 : 3 );
}
static inline int mm_bdlog_3x (mm_t a[3], mm_t ainv[3], mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    if ( b == R ) return 0;
    if ( b == a[0] ) return 1;
    if ( b == a[1] ) return 2;
    if ( b == a[2] ) return 4;
    if ( b == ainv[0] ) return 7;
    if ( b == ainv[1] ) return 6;
    return ( b == p-a[0] ? 5 : ( b == p-ainv[0] ? 3 : -1 ) );
}

static inline int mm_bdlog_4 (mm_t a[4], mm_t ainv[4], mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    int k0 = mm_bdlog_2 (a+2, ainv+2, mm_sqr(mm_sqr(b,p,pinv),p,pinv), R, p, pinv);
    return (mm_bdlog_2 (a+2, ainv+2, mm_mul (b, mm_exp_pow (ainv, k0, R, p, pinv), p, pinv), R, p, pinv) << 2) + k0;
}

static inline int mm_bdlog_4x (mm_t a[4], mm_t ainv[4], mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    int k0 = mm_bdlog_2x (a+2, ainv+2, mm_sqr(mm_sqr(b,p,pinv),p,pinv), R, p, pinv);
    if ( k0 < 0 ) return -1;
    return (mm_bdlog_2 (a+2, ainv+2, mm_mul (b, mm_exp_pow (ainv, k0, R, p, pinv), p, pinv), R, p, pinv) << 2) + k0;
}

static inline int mm_bdlog_5 (mm_t a[5], mm_t ainv[5], mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    int k0 = mm_bdlog_3 (a+2, ainv+2, mm_sqr(mm_sqr(b,p,pinv),p,pinv), R, p, pinv);
    return (mm_bdlog_2 (a+3, ainv+3, mm_mul (b, mm_exp_pow (ainv, k0, R, p, pinv), p, pinv), R, p, pinv) << 3) + k0;
}

static inline int mm_bdlog_5x (mm_t a[5], mm_t ainv[5], mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    int k0 = mm_bdlog_3x (a+2, ainv+2, mm_sqr(mm_sqr(b,p,pinv),p,pinv), R, p, pinv);
    if ( k0 < 0 ) return -1;
    return (mm_bdlog_2 (a+3, ainv+3, mm_mul (b, mm_exp_pow (ainv, k0, R, p, pinv), p, pinv), R, p, pinv) << 3) + k0;
}

static inline int mm_bdlog_6 (mm_t a[6], mm_t ainv[6], mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    int k0 = mm_bdlog_3 (a+3, ainv+3, mm_sqr(mm_sqr(mm_sqr(b,p,pinv),p,pinv),p,pinv), R, p, pinv);
    return (mm_bdlog_3 (a+3, ainv+3, mm_mul (b, mm_exp_pow (ainv, k0, R, p, pinv), p, pinv), R, p, pinv) << 3) + k0;
}

static inline int mm_bdlog_6x (mm_t a[6], mm_t ainv[6], mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    int k0 = mm_bdlog_3x (a+3, ainv+3, mm_sqr(mm_sqr(mm_sqr(b,p,pinv),p,pinv),p,pinv), R, p, pinv);
    if ( k0 < 0 ) return -1;
    return (mm_bdlog_3 (a+3, ainv+3, mm_mul (b, mm_exp_pow (ainv, k0, R, p, pinv), p, pinv), R, p, pinv) << 3) + k0;
}

static inline int mm_bdlog_7 (mm_t a[7], mm_t ainv[7], mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    int k0 = mm_bdlog_4 (a+3, ainv+3, mm_sqr(mm_sqr(mm_sqr(b,p,pinv),p,pinv),p,pinv), R, p, pinv);
    return (mm_bdlog_3 (a+4, ainv+4, mm_mul (b, mm_exp_pow (ainv, k0, R, p, pinv), p, pinv), R, p, pinv) << 4) + k0;
}

static inline int mm_bdlog_7x (mm_t a[7], mm_t ainv[7], mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    int k0 = mm_bdlog_4x (a+3, ainv+3, mm_sqr(mm_sqr(mm_sqr(b,p,pinv),p,pinv),p,pinv), R, p, pinv);
    if ( k0 < 0 ) return -1;
    return (mm_bdlog_3 (a+4, ainv+4, mm_mul (b, mm_exp_pow (ainv, k0, R, p, pinv), p, pinv), R, p, pinv) << 4) + k0;
}

/*
static inline int mm_bdlog_8 (mm_t a[8], mm_t ainv[8], mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    int k0 = mm_bdlog_4 (a+4, ainv+4, mm_sqr(mm_sqr(mm_sqr(mm_sqr(b,p,pinv),p,pinv),p,pinv),p,pinv), R, p, pinv);
    return (mm_bdlog_4 (a+4, ainv+4, mm_mul (b, mm_exp_pow (ainv, k0, R, p, pinv), p, pinv), R, p, pinv) << 4) + k0;
}
*/
static inline int mm_bdlog_8 (mm_t a[8], mm_t ainv[8], mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    int k0 = mm_bdlog_5 (a+3, ainv+3, mm_sqr(mm_sqr(mm_sqr(b,p,pinv),p,pinv),p,pinv), R, p, pinv);
    return (mm_bdlog_3 (a+5, ainv+5, mm_mul (b, mm_exp_pow (ainv, k0, R, p, pinv), p, pinv), R, p, pinv) << 5) + k0;
}

static inline int mm_bdlog_8x (mm_t a[8], mm_t ainv[8], mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    int k0 = mm_bdlog_5x (a+3, ainv+3, mm_sqr(mm_sqr(mm_sqr(b,p,pinv),p,pinv),p,pinv), R, p, pinv);
    if ( k0 < 0 ) return -1;
    return (mm_bdlog_3 (a+5, ainv+5, mm_mul (b, mm_exp_pow (ainv, k0, R, p, pinv), p, pinv), R, p, pinv) << 5) + k0;
}

// Recursive binary discrete log algorithm: given a of order 2^e and b in <a>, returnsk in [0,2^3) such that a^k = b
// Takes precomputed powers a[i]=a^(2^i) and ainv=a^{-2^i)
uint64_t mm_bdlog (mm_t a[], mm_t ainv[], int e, mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    switch (e) {
    case 1: return mm_bdlog_1 (a,ainv,b,R,p,pinv);
    case 2: return mm_bdlog_2 (a,ainv,b,R,p,pinv);
    case 3: return mm_bdlog_3 (a,ainv,b,R,p,pinv);
    case 4: return mm_bdlog_4 (a,ainv,b,R,p,pinv);
    case 5: return mm_bdlog_5 (a,ainv,b,R,p,pinv);
    case 6: return mm_bdlog_6 (a,ainv,b,R,p,pinv);
    case 7: return mm_bdlog_7 (a,ainv,b,R,p,pinv);
    case 8: return mm_bdlog_8 (a,ainv,b,R,p,pinv);
    }
    int d = e/2; e -= d;
    uint64_t k0 = mm_bdlog (a+d, ainv+d, e, mm_exp_2k (b, d, p, pinv), R, p, pinv);
    return ((uint64_t)mm_bdlog (a+e, ainv+e, d, mm_mul (b,mm_exp_pow (ainv, k0, R, p, pinv), p, pinv), R, p, pinv) << e) + k0;
}

uint64_t mm_bdlogx (mm_t a[], mm_t ainv[], int e, mm_t b, mm_t R, mm_t p, mm_t pinv)
{
    switch (e) {
    case 0: return b == R ? 0 : -1;
    case 1: return mm_bdlog_1x (a,ainv,b,R,p,pinv);
    case 2: return mm_bdlog_2x (a,ainv,b,R,p,pinv);
    case 3: return mm_bdlog_3x (a,ainv,b,R,p,pinv);
    case 4: return mm_bdlog_4x (a,ainv,b,R,p,pinv);
    case 5: return mm_bdlog_5x (a,ainv,b,R,p,pinv);
    case 6: return mm_bdlog_6x (a,ainv,b,R,p,pinv);
    case 7: return mm_bdlog_7x (a,ainv,b,R,p,pinv);
    case 8: return mm_bdlog_8x (a,ainv,b,R,p,pinv);
    }
    int d = e/2; e -= d;
    int64_t k0 = mm_bdlogx (a+d, ainv+d, e, mm_exp_2k (b, d, p, pinv), R, p, pinv);
    if ( k0 < 0 ) return -1;
    return ((uint64_t)mm_bdlog (a+e, ainv+e, d, mm_mul (b,mm_exp_pow (ainv, k0, R, p, pinv), p, pinv), R, p, pinv) << e) + k0;
}

// computes the square root of 1/x modulo p=2^e*m+1 (with m odd), given lists of binary powers of generator a for the 2-Sylow and ainv=1/a
// returns 0 if x does not have a square-root in Fp (if ext is set then sqrt(x) = z*sqrt(a[0]) in Fp2)
int mm_invsqrt (mm_t *z, mm_t x, mm_t a[], mm_t ainv[], int e, mm_t R, mm_t p, mm_t pinv, int ext)
{
    int64_t k, m;
    
    if ( ! x ) { *z = 0; return 1; }
    m = p >> e;                                                                     // p = 2^e*m+1 (we assume m is odd)
    mm_t y = mm_exp_ui (x, m>>1, R, p, pinv);                                       // y = x^((m-1)/2)
    mm_t b = mm_mul(mm_sqr(y,p,pinv),x,p,pinv);                                     // b = y^2*x = x^m is in the 2-Sylow
    if ( ext ) {
        k = mm_bdlog (a, ainv, e, b, R, p, pinv);                                   // k = dlog(a,b)
        *z = k ? mm_mul (y, mm_exp_pow (ainv, k>>1, R, p, pinv), p, pinv) : y;      // z = y*b^(-1/2) = x^((m-1)/2)*x^(-m/2) = x^(-1/2)
        if ( (k&1) ) *z = mm_mul(*z,ainv[0],p,pinv);
        return !(k&1);
    } else {
        k = mm_bdlogx (a+1, ainv+1, e-1, b, R, p, pinv);
        if ( k < 0 ) return 0;
        *z = k ? mm_mul (y, mm_exp_pow (ainv, k, R, p, pinv), p, pinv) : y;         // z = y*b^(-1/2) = x^((m-1)/2)*x^(-m/2) = x^(-1/2)
        return 1;
    }
}

static unsigned char _nr35_tab[35] =    { 0, 0, 5, 5, 0, 7, 7, 5, 5, 0, 7, 0, 5, 5, 0, 0, 0, 5, 5, 7, 7, 0, 5, 5, 7, 0, 7, 5, 5, 0, 0, 7, 5, 5, 7};
static unsigned char _p256_tab[55] = {  2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 0 };

// conditional sqrt -- returns r such that r^2=a conditional on a being square  (caller should check r^2=a if a is not known to be a residue)
// R2 may be zero (if so, it will computed as/when needed)
mm_t mm_sqrtc (mm_t a, mm_t R, mm_t R2, mm_t p, mm_t pinv)
{
    register mm_t r,b,x,z;
    uint64_t m;
    int d, e;

    if ( !a ) return 0;
    if ( (p&3) == 3 ) return mm_exp_ui (a,(p+1)>>2,R,p,pinv);
    if ( (p&7) == 5 ) {
        r = mm_exp_ui (a,(p+3)>>3,R,p,pinv);
        if ( mm_sqr (r,p,pinv) == a ) return r;
        z = mm_exp_ui(mm_add(R,R,p),(p-1)>>2,R,p,pinv); // z = 2^((p-1)/4) = sqrt(-1)
        // if a is square and r^2 != a then we must have r^2 = -a, use 2^((p-1)/4) to get sqrt(-1)
        return mm_mul(z,r,p,pinv);
    }
    e = __builtin_ctzl(p-1);  m = (p-1)>>e;         // p =2^e*m+1
    r = mm_exp_ui (a,(m-1)>>1,R,p,pinv);            // r = a^((m-1)/2), r^2 = a^m * a^-1
    b = mm_mul(a,mm_sqr(r,p,pinv),p,pinv);          // b = a^m is the image of a in the 2-Sylow
    r = mm_mul(r,a,p,pinv);                         // Now r^2 = b*a
    if ( b == R ) return r;                         // Good news, no need to find a nonresidue
    // Get a nonresidue mod p < 2^64, which we expect to be < 256 (not proved but surely true)
    z = ( (p%3) == 2 ? 3 : _nr35_tab[p%35] );
    if ( ! z ) for ( unsigned char *q = _p256_tab+5 ; ui_legendre(z=*q,p) > 0 ; q++ );  // *q is small so these should be quick
    assert (z); // check for counterexample to claim above (note _p256_tab is 0-terminated)
    z = mm_exp_ui (mm_from_ui(z,R2,p,pinv),m,R,p,pinv); // z generates 2-Sylow
    // From here use a standard O(e^2) Tonnelli-Shanks approach, which amounts to computing the
    // discrete log of b wrt z bit-by-bit (low to high)
    // This can be improved to O((e log e)/ loglog(e)) via https://www.ams.org/journals/mcom/2011-80-273/S0025-5718-10-02356-2/
    // But this would be useful only when v_2(p-1) is largish, say greater than 8, which is fairly rare, so we don't bother
    do {
        // z has order 2^e and generates subgroup of 2-Sylow containg b
        for ( d = 0, x = b ; x != R && d < e-1 ; d++ ) x = mm_sqr(x,p,pinv);
        if ( x != R ) return 0; // no square root
        while ( e > d+1 ) { z = mm_sqr(z,p,pinv); e--; } // power z to have order 2^(d+1)
        r = mm_mul(r,z,p,pinv); // adjust sqrt(b) in r (flip bottom bit of dlog)
        if ( d == 1 ) break;    // quit now if know that was the last bit
        z = mm_sqr(z,p,pinv); e--;
        b = mm_mul(b,z,p,pinv); // adjust b (flip bottom bit of dlog)
    } while ( b != R );
    return r;
}
