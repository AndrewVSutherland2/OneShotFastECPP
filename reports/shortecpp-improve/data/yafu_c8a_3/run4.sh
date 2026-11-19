#!/bin/bash
# yafu as the ECM engine on a 192-core box, one prime after the other, all cores each:
#   y768_k   the 13 primes of 768 bits of the gmp-ecm runs (rule_c8a/a768_k)
#   c<c>     nextprime(10^c), c = 160, 170, ..., 320: the entries of the ShortPrimalityProofs table
#   y1030_k  the 1030-bit primes k = 1, 2 of rule_c8a/b1030_k (k = 0 is c310)
#   then c = 10, ..., 150 (the entries of the table made with short.gp)
cd ~/rule/repo
. ./setenv.sh > /dev/null
OUT=~/rule/out4
Y=~/yb/yafu/yafu
mkdir -p $OUT
run() {  # tag prime timeout
  rm -f work/shortcamp/$1.*
  /usr/bin/time -f "RUN $1 wall=%e user=%U sys=%S" timeout -k 20 $3 python3 ecpp/shortECPP.py p=$2 threads=192 v2=1 tag=$1 out=$OUT/$1.cert.txt sched=rule engine=yafu yafu=$Y > $OUT/$1.log 2>&1
  echo "EXIT $?" >> $OUT/$1.log
  pkill -x ecm; pkill -x yafu; pkill -x dscan; pkill -x classpoly; pkill -x cm_method; pkill -x smoothtest; pkill -x gp
  echo "$(date +%T) $1 done: $(grep -o 'RUN.*' $OUT/$1.log) verified=$(grep -c 'verified twice' $OUT/$1.log)" >> $OUT/progress
  sleep 2
}
echo "$(date +%T) start on $(hostname)" >> $OUT/progress
k=0; while read P; do [ $k -le 12 ] && run y768_$k $P 1500; k=$((k+1)); done < ~/rule/p768.txt
while read C P; do [ $C -le 220 ] && run c$C $P 3600; done < ~/rule/pall.txt
echo "$(date +%T) ALL_DONE" >> $OUT/progress
