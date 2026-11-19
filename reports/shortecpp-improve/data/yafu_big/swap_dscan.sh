#!/bin/bash
# swap_dscan.sh (on a box, after dscan.c was shipped): build the scan program beside the one in use,
# check it against that one on a small scan, and put it in place with a rename
cd ~/rule/repo/ecpp || exit 1
F="-O3 -march=native -m64 -funroll-loops -std=gnu11 -Wall -Wextra -fopenmp -isystem ../classpoly_v1.0.3 -isystem ../local/include"
Z=$(ls ../classpoly_v1.0.3/zp_poly_small.o ../classpoly_v1.0.3/zp_poly_util.o ../classpoly_v1.0.3/zp_poly_mul.o ../classpoly_v1.0.3/zp_poly_div.o ../classpoly_v1.0.3/zp_poly_gcd.o ../classpoly_v1.0.3/zp_poly_exp.o ../classpoly_v1.0.3/zp_poly_comp.o ../classpoly_v1.0.3/zp2.o ../classpoly_v1.0.3/zp_roots.o)
gcc $F -c dscan.c -o dscan_new.o 2>&1 | tail -n 3
gcc $F -o dscan.new dscan_new.o cornacchia.o smooth.o $Z -lgmp -lm 2>&1 | tail -n 3
[ -x dscan.new ] || { echo "build failed"; exit 1; }
a=$(./dscan pbits=384 seed=3 B=30000000 threads=4 dump 2>/dev/null | sort | md5sum | cut -c1-12)
b=$(./dscan.new pbits=384 seed=3 B=30000000 threads=4 dump 2>/dev/null | sort | md5sum | cut -c1-12)
c=$(./dscan.new pbits=384 seed=3 B=30000000 threads=4 dump sd=50000 2>/dev/null | sort | md5sum | cut -c1-12)
n=$(./dscan.new pbits=384 seed=3 B=30000000 threads=4 dump 2>/dev/null | wc -l)
echo "old $a new $b new,sd=50000 $c lines $n"
if [ "$a" = "$b" ] && [ "$a" = "$c" ] && [ "$n" -gt 100 ]; then mv dscan.new dscan && echo "SWAPPED on $(hostname)"; else echo "MISMATCH: not swapped"; rm -f dscan.new; fi
