/*
   tuning.c:  global tuning values

   Copyright (C) 2007, 2008, David Harvey

   This file is part of the zn_poly library (version 0.9).

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) version 3 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "zn_poly_internal.h"

size_t ZNP_mpn_smp_kara_thresh = 50;
size_t ZNP_mpn_mulmid_fallback_thresh = 50;

tuning_info_t tuning_info[] = 
{
   {  // bits = 0
       0,0,0,0,0,0,0,0,0,0,0
   },
   {  // bits = 1
       0,0,0,0,0,0,0,0,0,0,0
   },
   {  // bits = 2
        302,   // KS1 -> KS2 multiplication threshold
   SIZE_MAX,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
        377,   // KS1 -> KS2 squaring threshold
   SIZE_MAX,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
        158,   // KS1 -> KS2 middle product threshold
   SIZE_MAX,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 3
        173,   // KS1 -> KS2 multiplication threshold
         90,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
        241,   // KS1 -> KS2 squaring threshold
       7878,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
        141,   // KS1 -> KS2 middle product threshold
       1030,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 4
        138,   // KS1 -> KS2 multiplication threshold
   SIZE_MAX,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
        173,   // KS1 -> KS2 squaring threshold
   SIZE_MAX,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
        132,   // KS1 -> KS2 middle product threshold
       1053,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 5
        123,   // KS1 -> KS2 multiplication threshold
       4505,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
        158,   // KS1 -> KS2 squaring threshold
       4212,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
        132,   // KS1 -> KS2 middle product threshold
       2106,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 6
        132,   // KS1 -> KS2 multiplication threshold
       2694,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
        132,   // KS1 -> KS2 squaring threshold
   SIZE_MAX,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
        102,   // KS1 -> KS2 middle product threshold
        705,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 7
         94,   // KS1 -> KS2 multiplication threshold
       5151,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
        132,   // KS1 -> KS2 squaring threshold
       3939,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         86,   // KS1 -> KS2 middle product threshold
        322,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 8
         86,   // KS1 -> KS2 multiplication threshold
       1053,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
        112,   // KS1 -> KS2 squaring threshold
       1378,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         78,   // KS1 -> KS2 middle product threshold
        421,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 9
         72,   // KS1 -> KS2 multiplication threshold
       1884,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
        102,   // KS1 -> KS2 squaring threshold
       2154,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         75,   // KS1 -> KS2 middle product threshold
        247,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 10
         60,   // KS1 -> KS2 multiplication threshold
        985,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         86,   // KS1 -> KS2 squaring threshold
       1127,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         66,   // KS1 -> KS2 middle product threshold
        394,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 11
         55,   // KS1 -> KS2 multiplication threshold
        461,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         72,   // KS1 -> KS2 squaring threshold
        901,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         61,   // KS1 -> KS2 middle product threshold
        216,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 12
         47,   // KS1 -> KS2 multiplication threshold
        576,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         62,   // KS1 -> KS2 squaring threshold
        901,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         51,   // KS1 -> KS2 middle product threshold
        221,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 13
         47,   // KS1 -> KS2 multiplication threshold
        345,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         56,   // KS1 -> KS2 squaring threshold
        515,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         51,   // KS1 -> KS2 middle product threshold
        189,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 14
         51,   // KS1 -> KS2 multiplication threshold
        385,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         50,   // KS1 -> KS2 squaring threshold
        689,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         47,   // KS1 -> KS2 middle product threshold
        173,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
         13    // nussbaumer squaring threshold
   },
   {  // bits = 15
         39,   // KS1 -> KS2 multiplication threshold
        322,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         47,   // KS1 -> KS2 squaring threshold
        345,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         39,   // KS1 -> KS2 middle product threshold
        161,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 16
         39,   // KS1 -> KS2 multiplication threshold
        302,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         47,   // KS1 -> KS2 squaring threshold
        539,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         43,   // KS1 -> KS2 middle product threshold
        151,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
         13    // nussbaumer squaring threshold
   },
   {  // bits = 17
         33,   // KS1 -> KS2 multiplication threshold
        282,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         39,   // KS1 -> KS2 squaring threshold
        288,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         36,   // KS1 -> KS2 middle product threshold
        151,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
       1000    // nussbaumer squaring threshold
   },
   {  // bits = 18
         33,   // KS1 -> KS2 multiplication threshold
        206,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         43,   // KS1 -> KS2 squaring threshold
        345,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         39,   // KS1 -> KS2 middle product threshold
        132,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         13,   // nussbaumer multiplication threshold
         13    // nussbaumer squaring threshold
   },
   {  // bits = 19
         33,   // KS1 -> KS2 multiplication threshold
        206,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         39,   // KS1 -> KS2 squaring threshold
        264,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         36,   // KS1 -> KS2 middle product threshold
        124,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
       1000,   // nussbaumer multiplication threshold
         13    // nussbaumer squaring threshold
   },
   {  // bits = 20
         31,   // KS1 -> KS2 multiplication threshold
        158,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         43,   // KS1 -> KS2 squaring threshold
        241,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         35,   // KS1 -> KS2 middle product threshold
        124,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         13,   // nussbaumer multiplication threshold
         13    // nussbaumer squaring threshold
   },
   {  // bits = 21
         30,   // KS1 -> KS2 multiplication threshold
        206,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         38,   // KS1 -> KS2 squaring threshold
        247,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         33,   // KS1 -> KS2 middle product threshold
        116,   // KS2 -> KS4 middle product threshold
      449378,   // KS4 -> FFT middle product threshold
         13,   // nussbaumer multiplication threshold
         13    // nussbaumer squaring threshold
   },
   {  // bits = 22
         31,   // KS1 -> KS2 multiplication threshold
        121,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         35,   // KS1 -> KS2 squaring threshold
        193,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         31,   // KS1 -> KS2 middle product threshold
        116,   // KS2 -> KS4 middle product threshold
      449378,   // KS4 -> FFT middle product threshold
         13,   // nussbaumer multiplication threshold
         13    // nussbaumer squaring threshold
   },
   {  // bits = 23
         31,   // KS1 -> KS2 multiplication threshold
        149,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         35,   // KS1 -> KS2 squaring threshold
        226,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         33,   // KS1 -> KS2 middle product threshold
        108,   // KS2 -> KS4 middle product threshold
      449378,   // KS4 -> FFT middle product threshold
         13,   // nussbaumer multiplication threshold
         13    // nussbaumer squaring threshold
   },
   {  // bits = 24
         27,   // KS1 -> KS2 multiplication threshold
         86,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         33,   // KS1 -> KS2 squaring threshold
        158,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         31,   // KS1 -> KS2 middle product threshold
        102,   // KS2 -> KS4 middle product threshold
      449378,   // KS4 -> FFT middle product threshold
         13,   // nussbaumer multiplication threshold
         13    // nussbaumer squaring threshold
   },
   {  // bits = 25
         25,   // KS1 -> KS2 multiplication threshold
        107,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         33,   // KS1 -> KS2 squaring threshold
        158,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         33,   // KS1 -> KS2 middle product threshold
        102,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         13,   // nussbaumer multiplication threshold
         13    // nussbaumer squaring threshold
   },
   {  // bits = 26
         31,   // KS1 -> KS2 multiplication threshold
        102,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         31,   // KS1 -> KS2 squaring threshold
        132,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         31,   // KS1 -> KS2 middle product threshold
         94,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         13,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 27
         25,   // KS1 -> KS2 multiplication threshold
         95,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         31,   // KS1 -> KS2 squaring threshold
        149,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         30,   // KS1 -> KS2 middle product threshold
         94,   // KS2 -> KS4 middle product threshold
      368640,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 28
         31,   // KS1 -> KS2 multiplication threshold
        107,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         31,   // KS1 -> KS2 squaring threshold
        158,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         29,   // KS1 -> KS2 middle product threshold
         94,   // KS2 -> KS4 middle product threshold
      407012,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 29
         19,   // KS1 -> KS2 multiplication threshold
         85,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         26,   // KS1 -> KS2 squaring threshold
        173,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         25,   // KS1 -> KS2 middle product threshold
         86,   // KS2 -> KS4 middle product threshold
      496153,   // KS4 -> FFT middle product threshold
         13,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 30
         19,   // KS1 -> KS2 multiplication threshold
         86,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         23,   // KS1 -> KS2 squaring threshold
        102,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         19,   // KS1 -> KS2 middle product threshold
         80,   // KS2 -> KS4 middle product threshold
      547797,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 31
         21,   // KS1 -> KS2 multiplication threshold
         85,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         25,   // KS1 -> KS2 squaring threshold
        119,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         23,   // KS1 -> KS2 middle product threshold
         86,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         13,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 32
         19,   // KS1 -> KS2 multiplication threshold
         86,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         25,   // KS1 -> KS2 squaring threshold
         94,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         23,   // KS1 -> KS2 middle product threshold
         70,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 33
         19,   // KS1 -> KS2 multiplication threshold
         72,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         25,   // KS1 -> KS2 squaring threshold
         94,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         19,   // KS1 -> KS2 middle product threshold
         80,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 34
         17,   // KS1 -> KS2 multiplication threshold
         65,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         23,   // KS1 -> KS2 squaring threshold
         78,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         21,   // KS1 -> KS2 middle product threshold
         70,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 35
         16,   // KS1 -> KS2 multiplication threshold
         75,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         23,   // KS1 -> KS2 squaring threshold
         94,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         19,   // KS1 -> KS2 middle product threshold
         75,   // KS2 -> KS4 middle product threshold
      333886,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 36
         17,   // KS1 -> KS2 multiplication threshold
         62,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         21,   // KS1 -> KS2 squaring threshold
         78,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         21,   // KS1 -> KS2 middle product threshold
         66,   // KS2 -> KS4 middle product threshold
      368640,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 37
         16,   // KS1 -> KS2 multiplication threshold
         51,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         21,   // KS1 -> KS2 squaring threshold
         78,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         19,   // KS1 -> KS2 middle product threshold
         66,   // KS2 -> KS4 middle product threshold
      302409,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 38
         16,   // KS1 -> KS2 multiplication threshold
         61,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         21,   // KS1 -> KS2 squaring threshold
         65,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         23,   // KS1 -> KS2 middle product threshold
         66,   // KS2 -> KS4 middle product threshold
      273899,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 39
         16,   // KS1 -> KS2 multiplication threshold
         51,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         23,   // KS1 -> KS2 squaring threshold
         55,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         17,   // KS1 -> KS2 middle product threshold
         66,   // KS2 -> KS4 middle product threshold
      273899,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 40
         16,   // KS1 -> KS2 multiplication threshold
         51,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         21,   // KS1 -> KS2 squaring threshold
         55,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         17,   // KS1 -> KS2 middle product threshold
         66,   // KS2 -> KS4 middle product threshold
      273899,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 41
         14,   // KS1 -> KS2 multiplication threshold
         43,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         19,   // KS1 -> KS2 squaring threshold
         75,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         19,   // KS1 -> KS2 middle product threshold
         57,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 42
         14,   // KS1 -> KS2 multiplication threshold
         43,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         19,   // KS1 -> KS2 squaring threshold
         61,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         17,   // KS1 -> KS2 middle product threshold
         61,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 43
         14,   // KS1 -> KS2 multiplication threshold
         34,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         21,   // KS1 -> KS2 squaring threshold
         54,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         16,   // KS1 -> KS2 middle product threshold
         51,   // KS2 -> KS4 middle product threshold
      224689,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 44
         13,   // KS1 -> KS2 multiplication threshold
         47,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         17,   // KS1 -> KS2 squaring threshold
         57,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         17,   // KS1 -> KS2 middle product threshold
         54,   // KS2 -> KS4 middle product threshold
      224689,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 45
         14,   // KS1 -> KS2 multiplication threshold
         39,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         19,   // KS1 -> KS2 squaring threshold
         47,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         16,   // KS1 -> KS2 middle product threshold
         51,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 46
         14,   // KS1 -> KS2 multiplication threshold
         36,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         19,   // KS1 -> KS2 squaring threshold
         47,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         17,   // KS1 -> KS2 middle product threshold
         51,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 47
         13,   // KS1 -> KS2 multiplication threshold
         30,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         19,   // KS1 -> KS2 squaring threshold
         43,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         16,   // KS1 -> KS2 middle product threshold
         47,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 48
         13,   // KS1 -> KS2 multiplication threshold
         33,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         19,   // KS1 -> KS2 squaring threshold
         56,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         16,   // KS1 -> KS2 middle product threshold
         47,   // KS2 -> KS4 middle product threshold
   SIZE_MAX,   // KS4 -> FFT middle product threshold
         12,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 49
         13,   // KS1 -> KS2 multiplication threshold
         29,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         19,   // KS1 -> KS2 squaring threshold
         43,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         14,   // KS1 -> KS2 middle product threshold
         47,   // KS2 -> KS4 middle product threshold
      16563,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 50
         13,   // KS1 -> KS2 multiplication threshold
         36,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         16,   // KS1 -> KS2 squaring threshold
         39,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         14,   // KS1 -> KS2 middle product threshold
         43,   // KS2 -> KS4 middle product threshold
      92160,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         12    // nussbaumer squaring threshold
   },
   {  // bits = 51
         12,   // KS1 -> KS2 multiplication threshold
         25,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         19,   // KS1 -> KS2 squaring threshold
         39,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         14,   // KS1 -> KS2 middle product threshold
         43,   // KS2 -> KS4 middle product threshold
      25439,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         11    // nussbaumer squaring threshold
   },
   {  // bits = 52
         12,   // KS1 -> KS2 multiplication threshold
         25,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         19,   // KS1 -> KS2 squaring threshold
         39,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         16,   // KS1 -> KS2 middle product threshold
         43,   // KS2 -> KS4 middle product threshold
       5760,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         11    // nussbaumer squaring threshold
   },
   {  // bits = 53
         13,   // KS1 -> KS2 multiplication threshold
         25,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         17,   // KS1 -> KS2 squaring threshold
         36,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         13,   // KS1 -> KS2 middle product threshold
         43,   // KS2 -> KS4 middle product threshold
      36574,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         11    // nussbaumer squaring threshold
   },
   {  // bits = 54
         13,   // KS1 -> KS2 multiplication threshold
         23,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         17,   // KS1 -> KS2 squaring threshold
         39,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         13,   // KS1 -> KS2 middle product threshold
         39,   // KS2 -> KS4 middle product threshold
       6360,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         11    // nussbaumer squaring threshold
   },
   {  // bits = 55
         12,   // KS1 -> KS2 multiplication threshold
         23,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         17,   // KS1 -> KS2 squaring threshold
         33,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         16,   // KS1 -> KS2 middle product threshold
         38,   // KS2 -> KS4 middle product threshold
       7022,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         11    // nussbaumer squaring threshold
   },
   {  // bits = 56
         12,   // KS1 -> KS2 multiplication threshold
         25,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         17,   // KS1 -> KS2 squaring threshold
         35,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         13,   // KS1 -> KS2 middle product threshold
         39,   // KS2 -> KS4 middle product threshold
       6360,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         11    // nussbaumer squaring threshold
   },
   {  // bits = 57
         12,   // KS1 -> KS2 multiplication threshold
         21,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         16,   // KS1 -> KS2 squaring threshold
         35,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         12,   // KS1 -> KS2 middle product threshold
         35,   // KS2 -> KS4 middle product threshold
       7022,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         11    // nussbaumer squaring threshold
   },
   {  // bits = 58
         12,   // KS1 -> KS2 multiplication threshold
         29,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         14,   // KS1 -> KS2 squaring threshold
         33,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         13,   // KS1 -> KS2 middle product threshold
         35,   // KS2 -> KS4 middle product threshold
       4280,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         11    // nussbaumer squaring threshold
   },
   {  // bits = 59
         12,   // KS1 -> KS2 multiplication threshold
         23,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         16,   // KS1 -> KS2 squaring threshold
         33,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         14,   // KS1 -> KS2 middle product threshold
         35,   // KS2 -> KS4 middle product threshold
       4280,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         11    // nussbaumer squaring threshold
   },
   {  // bits = 60
         10,   // KS1 -> KS2 multiplication threshold
         17,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         16,   // KS1 -> KS2 squaring threshold
         19,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         13,   // KS1 -> KS2 middle product threshold
         29,   // KS2 -> KS4 middle product threshold
       3077,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         10    // nussbaumer squaring threshold
   },
   {  // bits = 61
         10,   // KS1 -> KS2 multiplication threshold
         21,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         14,   // KS1 -> KS2 squaring threshold
         33,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         13,   // KS1 -> KS2 middle product threshold
         33,   // KS2 -> KS4 middle product threshold
       6794,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         11    // nussbaumer squaring threshold
   },
   {  // bits = 62
         13,   // KS1 -> KS2 multiplication threshold
         48,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         17,   // KS1 -> KS2 squaring threshold
         72,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         17,   // KS1 -> KS2 middle product threshold
         40,   // KS2 -> KS4 middle product threshold
       6360,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         11    // nussbaumer squaring threshold
   },
   {  // bits = 63
         10,   // KS1 -> KS2 multiplication threshold
         43,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         14,   // KS1 -> KS2 squaring threshold
         78,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         13,   // KS1 -> KS2 middle product threshold
         35,   // KS2 -> KS4 middle product threshold
       6360,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         11    // nussbaumer squaring threshold
   },
   {  // bits = 64
         10,   // KS1 -> KS2 multiplication threshold
         54,   // KS2 -> KS4 multiplication threshold
   SIZE_MAX,   // KS4 -> FFT multiplication threshold
         14,   // KS1 -> KS2 squaring threshold
         76,   // KS2 -> KS4 squaring threshold
   SIZE_MAX,   // KS4 -> FFT squaring threshold
         13,   // KS1 -> KS2 middle product threshold
         40,   // KS2 -> KS4 middle product threshold
       6360,   // KS4 -> FFT middle product threshold
         11,   // nussbaumer multiplication threshold
         11    // nussbaumer squaring threshold
   },
};

// end of file ****************************************************************
