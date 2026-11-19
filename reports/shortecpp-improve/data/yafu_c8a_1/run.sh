#!/bin/bash
# yafu as the ECM engine on a 192-core box, one prime after the other, all cores each:
#   y768_k   the 13 primes of 768 bits of the gmp-ecm runs (rule_c8a/a768_k)
#   c<c>     nextprime(10^c), c = 160, 170, ..., 320: the entries of the ShortPrimalityProofs table
#   y1030_k  the 1030-bit primes k = 1, 2 of rule_c8a/b1030_k (k = 0 is c310)
cd ~/rule/repo
. ./setenv.sh > /dev/null
OUT=~/rule/out2
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
python3 reports/shortecpp-improve/data/yafu_cost2.py $Y $OUT/yafu_cost2_c8a.json sizes=532,768,1030 deep=0 > $OUT/yafu_cost2_c8a.log 2>&1
echo "$(date +%T) cost table done" >> $OUT/progress
k=0; while read P; do [ $k -le 12 ] && run y768_$k $P 1500; k=$((k+1)); done < ~/rule/p768.txt
while read C P; do run c$C $P 3600; done < ~/rule/pbench.txt
k=0; while read P; do [ $k -ge 1 ] && [ $k -le 2 ] && run y1030_$k $P 3600; k=$((k+1)); done < ~/rule/p1030.txt
echo "$(date +%T) ALL_DONE" >> $OUT/progress
