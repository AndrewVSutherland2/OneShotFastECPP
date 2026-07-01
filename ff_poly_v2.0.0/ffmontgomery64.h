#ifndef _FFMONTGOMERY64_INCLUDE_
#define _FFMONTGOMERY64_INCLUDE_

#include <stdint.h>

/*
    Copyright 2009-2017 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#ifdef __cplusplus
extern "C" {
#endif

#define _redc(a1,a0,x1,x0)  { a0=x0*_ff_mont_pni;_asm_mult_1_1 (a1,a0,a0,_ff_p); _asm_addto_2_2 (a1,a0,x1,x0);  if ( a1>=_ff_p ) a1 -= _ff_p;  }

// Note the functions below all require their inputs to be recuded mod p !!
// Additionally, p must be small enough for the unreduced sum of the reduced products to fit in 64 bits (this means less than 59 bits for the largest cases)

static inline uint64_t ff_montgomery1_mult (uint64_t x, uint64_t y)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x,y);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0x1, assumes _ff_p is less than 63 bits
// I have no idea why this is faster than multiplying and then doubling, but it appears to be (on an AMD Athlon, YMMV)
static inline uint64_t ff_montgomery1_sum_2_mults_s2 (uint64_t x0, uint64_t x1)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x1);
    _asm_addto_2_2 (b1,b0,b1,b0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0y1+y0x1 assumes _ff_p is less than 63 bits
static inline uint64_t ff_montgomery1_sum_2_mults (uint64_t x0, uint64_t x1, uint64_t y0, uint64_t y1)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y1);
    _asm_mult_1_1 (a1,a0,x1,y0);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0y1-y0x1 assumes _ff_p is less than 63 bits
static inline uint64_t ff_montgomery1_diff_2_mults (uint64_t x0, uint64_t x1, uint64_t y0, uint64_t y1)
{
    register uint64_t b0, b1,a0, a1,t;

    t = _ff_p-y0;
    _asm_mult_1_1 (b1,b0,x0,y1);
    _asm_mult_1_1 (a1,a0,x1,t);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y1+y0x1 assumes _ff_p is less than 62 bits
static inline uint64_t ff_montgomery1_sum_2_mults_d1 (uint64_t x0, uint64_t x1, uint64_t y0, uint64_t y1)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y1);
    _asm_addto_2_2 (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y0);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0y2+ x1y1 + x2y0  assumes _ff_p is less than 62 bits
static inline uint64_t ff_montgomery1_sum_3_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t y0, uint64_t y1, uint64_t y2)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y2);
    _asm_mult_1_1 (a1,a0,x1,y1);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y0);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y2+ x1y1 + x2y0  assumes _ff_p is less than 62 bits
static inline uint64_t ff_montgomery1_sum_3_mults_d1 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t y0, uint64_t y1, uint64_t y2)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y2);
    _asm_addto_2_2 (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y1);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y0);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y2+ 2x1y1 + x2y0  assumes _ff_p is less than 62 bits
static inline uint64_t ff_montgomery1_sum_3_mults_d2 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t y0, uint64_t y1, uint64_t y2)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y2);
    _asm_addto_2_2 (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y1);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y0);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0x2+ x1x1 + x2x0 = 2x0x2+x1^2,  assumes _ff_p is less than 62 bits
static inline uint64_t ff_montgomery1_sum_3_mults_s3 (uint64_t x0, uint64_t x1, uint64_t x2)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x2);
    _asm_addto_2_2 (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x1);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0y3+x1y2+x2y1+x3y0,  assumes _ff_p is less than 62 bits
static inline uint64_t ff_montgomery1_sum_4_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y3);
    _asm_mult_1_1 (a1,a0,x1,y2);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y1);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y0);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y3+x1y2+x2y1+x3y0,  assumes _ff_p is less than 62 bits
static inline uint64_t ff_montgomery1_sum_4_mults_d1 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y3);
    _asm_addto_2_2 (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y2);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y1);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y0);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


// computes x0x3+x1x2+x2x1+x3x0 = 2x0x3+2x1x2,  assumes _ff_p is less than 62 bits
static inline uint64_t ff_montgomery1_sum_4_mults_s4 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x3);
    _asm_addto_2_2 (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x2);
    _asm_addto_2_2 (a1,a0,a1,a0);
    _asm_addto_2_2 (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// for p < 2^61, we only need this if we are exceeding 56*p^2
#define _asm_addto_2_2_redp(b1,b0,a1,a0)        { _asm_addto_2_2(b1,b0,a1,a0);  while ( b1 >= _ff_p ) b1 -= _ff_p;}

// computes 2x0y3+2x1y2+x2y1+x3y0
static inline uint64_t ff_montgomery1_sum_4_mults_d2 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y3);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y2);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0y4+x1y3+x2y2+x3y1+x4y0
static inline uint64_t ff_montgomery1_sum_5_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4,
                                                  uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y4);
    _asm_mult_1_1 (a1,a0,x1,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y4+x1y3+x2y2+x3y1+x4y0
static inline uint64_t ff_montgomery1_sum_5_mults_d1 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4,
                                                       uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y4);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0x4+x1x3+x2x2+x3x1+x4x0 = 2x0x4+2x1x3+x2^2
static inline uint64_t ff_montgomery1_sum_5_mults_s5 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x4);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x3);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0y5+x1y4+x2y3+x3y2+x4y1+x5y0
static inline uint64_t ff_montgomery1_sum_6_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                  uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y5);
    _asm_mult_1_1 (a1,a0,x1,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y5+x1y4+x2y3+x3y2+x4y1+x5y0
static inline uint64_t ff_montgomery1_sum_6_mults_d1 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                       uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y5);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y5+2x1y4+x2y3+x3y2+x4y1+x5y0
static inline uint64_t ff_montgomery1_sum_6_mults_d2 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                  uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y5);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y4);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0x5+x1x4+x2x3+x3x2+x4x1+x5x0 = 2x0x5+2x1x4+2x2x3
static inline uint64_t ff_montgomery1_sum_6_mults_s6 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x5);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x4);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x3);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0y6+x1y5+x2y4+x3y3+x4y2+x5y1+x6y0
static inline uint64_t ff_montgomery1_sum_7_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6,
                                                  uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y6);
    _asm_mult_1_1 (a1,a0,x1,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y6+2x1y5+x2y4+x3y3+x4y2+x5y1+x6y0
static inline uint64_t ff_montgomery1_sum_7_mults_d2 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6,
                                                  uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y6);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y5);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y6+2x1y5+2x2y4+x3y3+x4y2+x5y1+x6y0
static inline uint64_t ff_montgomery1_sum_7_mults_d3 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6,
                                                  uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y6);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y5);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y4);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0x6+x1x5+x2x4+x3x3+x4x2+x5x1+x6x0 = 2x0x6+2x1x5+2x2x4+x3^2
static inline uint64_t ff_montgomery1_sum_7_mults_s7 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x6);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x5);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x4);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0y7+x1y6+x2y5+x3y4+x4y3+x5y2+x6y1+x7y0
static inline uint64_t ff_montgomery1_sum_8_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7,
                                                  uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y7);
    _asm_mult_1_1 (a1,a0,x1,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y7+2x1y6+x2y5+x3y4+x4y3+x5y2+x6y1+x7y0
static inline uint64_t ff_montgomery1_sum_8_mults_d2 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7,
                                                  uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y7);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y6);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0x7+x1x6+x2x5+x3x4+x4x3+x5x2+x6x1+x7x0 = 2x0x7+2x1x6+2x2x5+2x3x4
static inline uint64_t ff_montgomery1_sum_8_mults_s8 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7)
{
    register uint64_t b0, b1,a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x7);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x6);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x5);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x4);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0y8+x1y7+x2y6+x3y5+x4y4+x5y3+x6y2+x7y1+x8y0
static inline uint64_t ff_montgomery1_sum_9_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                  uint64_t x6, uint64_t x7, uint64_t x8,
                                                  uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                  uint64_t y6, uint64_t y7, uint64_t y8)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y8);
    _asm_mult_1_1 (a1,a0,x1,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y8+2x1y7+x2y6+x3y5+x4y4+x5y3+x6y2+x7y1+x8y0
static inline uint64_t ff_montgomery1_sum_9_mults_d2 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                  uint64_t x6, uint64_t x7, uint64_t x8,
                                                  uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                  uint64_t y6, uint64_t y7, uint64_t y8)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y8);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y7);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y8+2x1y7+2x2y6+x3y5+x4y4+x5y3+x6y2+x7y1+x8y0
static inline uint64_t ff_montgomery1_sum_9_mults_d3 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                  uint64_t x6, uint64_t x7, uint64_t x8,
                                                  uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                  uint64_t y6, uint64_t y7, uint64_t y8)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y8);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y7);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y6);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0x8+x1x7+x2x6+x3x5+x4x4+x5x3+x6x2+x7x1+x8x0 = 2x0x8+2x1x7+2x2x6+2x3x5+x4^2
static inline uint64_t ff_montgomery1_sum_9_mults_s9 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                       uint64_t x6, uint64_t x7, uint64_t x8)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x8);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x7);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x6);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x5);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0y9+x1y8+x2y7+x3y6+x4y5+x5y4+x6y3+x7y2+x8y1+x9y0
static inline uint64_t ff_montgomery1_sum_10_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                    uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                   uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                           uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y9);
    _asm_mult_1_1 (a1,a0,x1,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y9+2x1y8+2x2y7+x3y6+x4y5+x5y4+x6y3+x7y2+x8y1+x9y0
static inline uint64_t ff_montgomery1_sum_10_mults_d3 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                    uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                   uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                           uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y9);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y8);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y7);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes 2x0y9+2x1y8+2x2y7+2x3y6+x4y5+x5y4+x6y3+x7y2+x8y1+x9y0
static inline uint64_t ff_montgomery1_sum_10_mults_d4 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                    uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                   uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                           uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y9);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y8);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y7);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y6);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// computes x0x9+x1x8+x2x7+x3x6+x4x5+x5x4+x6x3+x7x2+x8x1+x9x0 = 2x0x9+2x1x8+2x2x7+2x3x6+2x4x5
static inline uint64_t ff_montgomery1_sum_10_mults_s10 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                       uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x9);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x8);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x7);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x6);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x5);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_11_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                    uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10,
                                                   uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                           uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, uint64_t y10)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y10);
    _asm_mult_1_1 (a1,a0,x1,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_11_mults_d3 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                    uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10,
                                                   uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                           uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, uint64_t y10)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y10);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y9);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y8);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_11_mults_s11 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                       uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x10);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x9);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x8);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x7);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x6);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_12_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                    uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11,
                                                   uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                           uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y11);
    _asm_mult_1_1 (a1,a0,x1,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_12_mults_d3 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                    uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11,
                                                   uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                           uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y11);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y9);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_12_mults_d4 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                    uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11,
                                                   uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                           uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y11);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y9);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y8);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_12_mults_s12 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                       uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x11);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x9);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x8);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x7);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x6);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_13_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                    uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12,
                                                   uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                           uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y12);
    _asm_mult_1_1 (a1,a0,x1,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_13_mults_d4 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                    uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12,
                                                   uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                           uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y12);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y9);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_13_mults_d5 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                    uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12,
                                                   uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5,
                                                           uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y12);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y9);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y8);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_13_mults_s13 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5,
                                                       uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x12);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x9);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x8);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x7);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_14_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6,
                                                    uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6,
                                                    uint64_t y7, uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y13);
    _asm_mult_1_1 (a1,a0,x1,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_14_mults_d4 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6,
                                                    uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6,
                                                    uint64_t y7, uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y13);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_14_mults_s14 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6,
                                                       uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x13);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x9);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x8);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x7);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_15_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6,
                                                    uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6,
                                                    uint64_t y7, uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y14);
    _asm_mult_1_1 (a1,a0,x1,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_15_mults_d4 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6,
                                                    uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6,
                                                    uint64_t y7, uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y14);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_15_mults_d5 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6,
                                                    uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6,
                                                    uint64_t y7, uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y14);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_15_mults_s15 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6,
                                                       uint64_t x7, uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x14);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x9);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x8);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_16_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7,
                                                    uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7,
                                                    uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y15);
    _asm_mult_1_1 (a1,a0,x1,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_16_mults_d5 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7,
                                                    uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7,
                                                    uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y15);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_16_mults_d6 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7,
                                                    uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7,
                                                    uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y15);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_16_mults_s16 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7,
                                                       uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x15);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x9);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x8);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_17_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7,
                                                    uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7,
                                                    uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y16);
    _asm_mult_1_1 (a1,a0,x1,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}
static inline uint64_t ff_montgomery1_sum_17_mults_d5 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7,
                                                    uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7,
                                                    uint64_t y8, uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y16);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_17_mults_s17 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7,
                                                       uint64_t x8, uint64_t x9, uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x16);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x9);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_18_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8,
                                                    uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y17);
    _asm_mult_1_1 (a1,a0,x1,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_18_mults_d5 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8,
                                                    uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y17);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_18_mults_d6 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8,
                                                    uint64_t y9, uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y17);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_18_mults_s18 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x17);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x9);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


static inline uint64_t ff_montgomery1_sum_19_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y18);
    _asm_mult_1_1 (a1,a0,x1,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_19_mults_d6 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y18);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_19_mults_d7 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y18);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_19_mults_s19 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x18);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,x9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_20_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y19);
    _asm_mult_1_1 (a1,a0,x1,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_20_mults_d6 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y19);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_20_mults_s20 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x19);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,x10);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_21_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y20);
    _asm_mult_1_1 (a1,a0,x1,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_21_mults_d6 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y20);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_21_mults_d7 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y20);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_21_mults_s21 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x20);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,x11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,x10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_22_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y21);
    _asm_mult_1_1 (a1,a0,x1,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_22_mults_d7 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y21);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_22_mults_d8 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y21);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_22_mults_s22 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x21);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,x12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,x11);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_23_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21, uint64_t y22)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y22);
    _asm_mult_1_1 (a1,a0,x1,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_23_mults_d7 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21, uint64_t y22)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y22);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_23_mults_s23 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x22);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,x13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,x12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,x11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_24_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y23);
    _asm_mult_1_1 (a1,a0,x1,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_24_mults_d7 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y23);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_24_mults_d8 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y23);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


static inline uint64_t ff_montgomery1_sum_24_mults_s24 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x23);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,x14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,x13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,x12);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


static inline uint64_t ff_montgomery1_sum_25_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y24);
    _asm_mult_1_1 (a1,a0,x1,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_25_mults_d8 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y24);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_25_mults_d9 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y24);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_25_mults_s25 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x24);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,x14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,x13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,x12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}



static inline uint64_t ff_montgomery1_sum_26_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y25);
    _asm_mult_1_1 (a1,a0,x1,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_26_mults_d8 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y25);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_26_mults_s26 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                       uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x25);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,x14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,x13);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


static inline uint64_t ff_montgomery1_sum_27_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y26);
    _asm_mult_1_1 (a1,a0,x1,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_27_mults_d8 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y26);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

/*
    For no apparent reason, doubling (a1,a0) 8 times generates substantially faster code than doubling (b1,b0) once.
*/
static inline uint64_t ff_montgomery1_sum_27_mults_d9 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y26);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_27_mults_s27 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                       uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x26);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,x17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,x14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,x13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


static inline uint64_t ff_montgomery1_sum_28_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y27);
    _asm_mult_1_1 (a1,a0,x1,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_28_mults_d9 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y27);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_28_mults_d10 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y27);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_28_mults_s28 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                       uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x27);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,x18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,x17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,x14);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


static inline uint64_t ff_montgomery1_sum_29_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y28);
    _asm_mult_1_1 (a1,a0,x1,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_29_mults_d9 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y28);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_29_mults_s29 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                       uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x28);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,x19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,x18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,x17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,x14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_30_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y29);
    _asm_mult_1_1 (a1,a0,x1,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_30_mults_d9 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y29);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_30_mults_d10 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y29);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_30_mults_s30 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                       uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x29);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,x20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,x19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,x18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,x17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,x15);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_31_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, uint64_t x30,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, uint64_t y30)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y30);
    _asm_mult_1_1 (a1,a0,x1,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


static inline uint64_t ff_montgomery1_sum_31_mults_d10 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, uint64_t x30,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, uint64_t y30)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y30);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}
static inline uint64_t ff_montgomery1_sum_31_mults_d11 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, uint64_t x30,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, uint64_t y30)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y30);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_31_mults_s31 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                       uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                       uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, uint64_t x30)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,x30);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,x29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,x28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,x27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,x26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,x25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,x24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,x23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,x22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,x21);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,x20);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,x19);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,x18);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,x17);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,x16);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,x15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


static inline uint64_t ff_montgomery1_sum_32_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, uint64_t x30, uint64_t x31,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, uint64_t y30, uint64_t y31)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y31);
    _asm_mult_1_1 (a1,a0,x1,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


static inline uint64_t ff_montgomery1_sum_32_mults_d10 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, uint64_t x30, uint64_t x31,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, uint64_t y30, uint64_t y31)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y31);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


static inline uint64_t ff_montgomery1_sum_33_mults (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, uint64_t x30, uint64_t x31, uint64_t x32,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, uint64_t y30, uint64_t y31, uint64_t y32)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y32);
    _asm_mult_1_1 (a1,a0,x1,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


static inline uint64_t ff_montgomery1_sum_33_mults_d10 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, uint64_t x30, uint64_t x31, uint64_t x32,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, uint64_t y30, uint64_t y31, uint64_t y32)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y32);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


static inline uint64_t ff_montgomery1_sum_33_mults_d11 (uint64_t x0, uint64_t x1, uint64_t x2,uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9,
                                                    uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19,
                                                    uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, uint64_t x30, uint64_t x31, uint64_t x32,
                                                    uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9,
                                                    uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19,
                                                    uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, uint64_t y30, uint64_t y31, uint64_t y32)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y32);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_34_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y33);
    _asm_mult_1_1 (a1,a0,x1,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_34_mults_d11 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y33);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_34_mults_d12 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y33);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y22);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_35_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y34);
    _asm_mult_1_1 (a1,a0,x1,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_35_mults_d11 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y34);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_35_mults_d12 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y34);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y23);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_36_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y35);
    _asm_mult_1_1 (a1,a0,x1,y34);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_36_mults_d11 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y35);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_36_mults_d12 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y35);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}
static inline uint64_t ff_montgomery1_sum_37_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y36);
    _asm_mult_1_1 (a1,a0,x1,y35);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y34);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_37_mults_d12 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y36);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_37_mults_d13 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y36);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y24);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_38_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y37);
    _asm_mult_1_1 (a1,a0,x1,y36);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y35);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y34);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_38_mults_d12 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y37);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_38_mults_d13 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y37);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y25);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_39_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y38);
    _asm_mult_1_1 (a1,a0,x1,y37);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y36);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y35);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y34);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_39_mults_d12 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y38);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_39_mults_d13 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y38);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_40_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y39);
    _asm_mult_1_1 (a1,a0,x1,y38);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y37);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y36);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y35);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y34);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_40_mults_d13 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y39);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_40_mults_d14 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y39);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y26);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_41_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, uint64_t y40)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y40);
    _asm_mult_1_1 (a1,a0,x1,y39);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y38);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y37);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y36);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y35);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y34);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_41_mults_d13 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, uint64_t y40)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y40);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_41_mults_d14 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, uint64_t y40)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y40);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y27);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_42_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y41);
    _asm_mult_1_1 (a1,a0,x1,y40);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y39);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y38);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y37);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y36);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y35);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y34);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_42_mults_d13 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y41);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_42_mults_d14 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y41);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}


static inline uint64_t ff_montgomery1_sum_43_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y42);
    _asm_mult_1_1 (a1,a0,x1,y41);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y40);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y39);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y38);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y37);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y36);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y35);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y34);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_43_mults_d14 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y42);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y41);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_43_mults_d15 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y42);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y41);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y28);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_44_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y43);
    _asm_mult_1_1 (a1,a0,x1,y42);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y41);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y40);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y39);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y38);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y37);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y36);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y35);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y34);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_44_mults_d14 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y43);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y42);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y41);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_44_mults_d15 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y43);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y42);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y41);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y29);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_45_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t x44, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43, uint64_t y44)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y44);
    _asm_mult_1_1 (a1,a0,x1,y43);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y42);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y41);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y40);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y39);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y38);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y37);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y36);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y35);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y34);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x44,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_45_mults_d14 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t x44, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43, uint64_t y44)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y44);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y43);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y42);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y41);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x44,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_45_mults_d15 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t x44, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43, uint64_t y44)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y44);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y43);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y42);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y41);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x44,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_46_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t x44, uint64_t x45, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43, uint64_t y44, uint64_t y45)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y45);
    _asm_mult_1_1 (a1,a0,x1,y44);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y43);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y42);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y41);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y40);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y39);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y38);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y37);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y36);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y35);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y34);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x44,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x45,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_46_mults_d15 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t x44, uint64_t x45, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43, uint64_t y44, uint64_t y45)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y45);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y44);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y43);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y42);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y41);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x44,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x45,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_46_mults_d16 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t x44, uint64_t x45, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43, uint64_t y44, uint64_t y45)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y45);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y44);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y43);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y42);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y41);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y30);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x44,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x45,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_47_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t x44, uint64_t x45, uint64_t x46, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43, uint64_t y44, uint64_t y45, uint64_t y46)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y46);
    _asm_mult_1_1 (a1,a0,x1,y45);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y44);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y43);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y42);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y41);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y40);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y39);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y38);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y37);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y36);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y35);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y34);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x44,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x45,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x46,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_47_mults_d15 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t x44, uint64_t x45, uint64_t x46, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43, uint64_t y44, uint64_t y45, uint64_t y46)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y46);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y45);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y44);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y43);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y42);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y41);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x44,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x45,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x46,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_47_mults_d16 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t x44, uint64_t x45, uint64_t x46, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43, uint64_t y44, uint64_t y45, uint64_t y46)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y46);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y45);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y44);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y43);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y42);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y41);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y31);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x44,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x45,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x46,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_48_mults (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t x44, uint64_t x45, uint64_t x46, uint64_t x47, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43, uint64_t y44, uint64_t y45, uint64_t y46, uint64_t y47)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y47);
    _asm_mult_1_1 (a1,a0,x1,y46);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y45);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y44);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y43);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y42);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y41);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y40);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y39);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y38);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y37);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y36);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y35);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y34);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y33);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x44,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x45,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x46,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x47,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_48_mults_d15 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t x44, uint64_t x45, uint64_t x46, uint64_t x47, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43, uint64_t y44, uint64_t y45, uint64_t y46, uint64_t y47)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y47);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y46);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y45);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y44);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y43);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y42);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y41);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y32);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x44,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x45,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x46,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x47,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

static inline uint64_t ff_montgomery1_sum_48_mults_d16 (uint64_t x0, uint64_t x1, uint64_t x2, uint64_t x3, uint64_t x4, uint64_t x5, uint64_t x6, uint64_t x7, uint64_t x8, uint64_t x9, 
                                                         uint64_t x10, uint64_t x11, uint64_t x12, uint64_t x13, uint64_t x14, uint64_t x15, uint64_t x16, uint64_t x17, uint64_t x18, uint64_t x19, 
                                                         uint64_t x20, uint64_t x21, uint64_t x22, uint64_t x23, uint64_t x24, uint64_t x25, uint64_t x26, uint64_t x27, uint64_t x28, uint64_t x29, 
                                                         uint64_t x30, uint64_t x31, uint64_t x32, uint64_t x33, uint64_t x34, uint64_t x35, uint64_t x36, uint64_t x37, uint64_t x38, uint64_t x39, 
                                                         uint64_t x40, uint64_t x41, uint64_t x42, uint64_t x43, uint64_t x44, uint64_t x45, uint64_t x46, uint64_t x47, uint64_t y0, uint64_t y1, uint64_t y2, uint64_t y3, uint64_t y4, uint64_t y5, uint64_t y6, uint64_t y7, uint64_t y8, uint64_t y9, 
                                                         uint64_t y10, uint64_t y11, uint64_t y12, uint64_t y13, uint64_t y14, uint64_t y15, uint64_t y16, uint64_t y17, uint64_t y18, uint64_t y19, 
                                                         uint64_t y20, uint64_t y21, uint64_t y22, uint64_t y23, uint64_t y24, uint64_t y25, uint64_t y26, uint64_t y27, uint64_t y28, uint64_t y29, 
                                                         uint64_t y30, uint64_t y31, uint64_t y32, uint64_t y33, uint64_t y34, uint64_t y35, uint64_t y36, uint64_t y37, uint64_t y38, uint64_t y39, 
                                                         uint64_t y40, uint64_t y41, uint64_t y42, uint64_t y43, uint64_t y44, uint64_t y45, uint64_t y46, uint64_t y47)
{
    register uint64_t b0, b1, a0, a1;

    _asm_mult_1_1 (b1,b0,x0,y47);
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    _asm_mult_1_1 (a1,a0,x1,y46);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x2,y45);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x3,y44);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x4,y43);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x5,y42);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x6,y41);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x7,y40);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x8,y39);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x9,y38);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x10,y37);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x11,y36);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x12,y35);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x13,y34);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x14,y33);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x15,y32);
    _asm_addto_2_2_redp (a1,a0,a1,a0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x16,y31);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x17,y30);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x18,y29);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x19,y28);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x20,y27);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x21,y26);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x22,y25);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x23,y24);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x24,y23);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x25,y22);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x26,y21);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x27,y20);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x28,y19);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x29,y18);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x30,y17);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x31,y16);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x32,y15);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x33,y14);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x34,y13);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x35,y12);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x36,y11);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x37,y10);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x38,y9);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x39,y8);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x40,y7);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x41,y6);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x42,y5);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x43,y4);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x44,y3);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x45,y2);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x46,y1);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _asm_mult_1_1 (a1,a0,x47,y0);
    _asm_addto_2_2_redp (b1,b0,a1,a0);
    _redc(a1,a0,b1,b0);
    return a1;
}

// this is about 20% slower than the unwound versions above.
static inline uint64_t ff_montgomery1_sum_mults (uint64_t x[], uint64_t y[], int n)
{
    register int i;
    register uint64_t b0, b1, a0, a1, c1, c2;

    n--;
    b0=b1=0;
    for ( i = 0 ; i <= n  ; i++ ) {
        c1 = x[i]; c2=y[n-i];
        _asm_mult_1_1(a1,a0,c1,c2);
        _asm_addto_2_2_redp (b1,b0,a1,a0);
    }
    _redc(a1,a0,b1,b0);
    return a1;
}

// returns coeff of t^k in (x[0]+x[1]t+x[2]t^2+...)^2.  Assumes x has at least k+1 entries (degree at least k)
static inline uint64_t ff_montgomery1_sq_coeff (uint64_t x[], int k)
{
    register int i, j;
    register uint64_t b0, b1, a0, a1, c1, c2;

    b0=b1=0;
    j = (k+1)>>1;
    for ( i = 0 ; i < j ; i++ ) {
        c1 = x[i]; c2=x[k-i];
        _asm_mult_1_1(a1,a0,c1,c2);
        _asm_addto_2_2_redp (b1,b0,a1,a0);
    }
    _asm_addto_2_2_redp (b1,b0,b1,b0);
    if ( i==(k>>1) ) {      // true exactly when k is even
        c1 = x[i]; c2=x[i];
        _asm_mult_1_1(a1,a0,c1,c2);
        _asm_addto_2_2_redp (b1,b0,a1,a0);
    }
    _redc(a1,a0,b1,b0);
    return a1;
}

#ifdef __cplusplus
}
#endif

#endif
