#!/bin/bash
# one c8a.48xlarge, yafu, revised format, one prime after the other, all cores each:
#   y768_k  the 13 primes of 768 bits again (out5): does keeping a costly winner in reserve remove the long builds?
#   c<c>    nextprime(10^c), c = 330, ..., 400: the rest of the ShortPrimalityProofs table, then c = 410 (out6)
cd ~/rule/repo
. ./setenv.sh > /dev/null
Y=~/yb/yafu/yafu
run() {  # tag prime timeout outdir
  OUT=$4
  mkdir -p $OUT
  rm -f work/shortcamp/$1.*
  /usr/bin/time -f "RUN $1 wall=%e user=%U sys=%S" timeout -k 20 $3 python3 ecpp/shortECPP.py p=$2 threads=192 v2=1 tag=$1 out=$OUT/$1.cert.txt sched=rule engine=yafu yafu=$Y resume=1 > $OUT/$1.log 2>&1
  echo "EXIT $?" >> $OUT/$1.log
  pkill -x ecm; pkill -x yafu; pkill -x dscan; pkill -x classpoly; pkill -x cm_method; pkill -x smoothtest; pkill -x gp
  echo "$(date +%T) $1 done: $(grep -o 'RUN.*' $OUT/$1.log) verified=$(grep -c 'verified twice' $OUT/$1.log)" >> $OUT/progress
  sleep 2
}
echo "$(date +%T) start on $(hostname)" >> ~/rule/progress5
k=0; while read P; do [ $k -le 12 ] && run y768_$k $P 1500 ~/rule/out5; k=$((k+1)); done < ~/rule/p768.txt
echo "$(date +%T) y768 done" >> ~/rule/progress5
while read C P; do run c$C $P 43200 ~/rule/out6; done < ~/rule/pbig.txt
echo "$(date +%T) ALL_DONE" >> ~/rule/progress5
