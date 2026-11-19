#!/bin/bash
# stopping-rule scheduler on a 192-core box: 13 primes of 768 bits (the ones run on the 16-core box),
# three of 1030 bits, then 13 more of 768 bits; one after the other, all cores each.
cd ~/rule/repo
OUT=~/rule/out
mkdir -p $OUT
run() {  # tag prime timeout
  rm -f work/shortcamp/$1.*
  /usr/bin/time -f "RUN $1 wall=%e user=%U sys=%S" timeout -k 20 $3 python3 ecpp/shortECPP.py p=$2 threads=192 v2=1 tag=$1 out=$OUT/$1.cert.txt sched=rule > $OUT/$1.log 2>&1
  echo "EXIT $?" >> $OUT/$1.log
  pkill -x ecm; pkill -x dscan; pkill -x classpoly; pkill -x cm_method; pkill -x smoothtest; pkill -x gp
  echo "$(date +%T) $1 done" >> $OUT/progress
  sleep 2
}
echo "$(date +%T) start on $(hostname)" >> $OUT/progress
k=0; while read P; do [ $k -le 12 ] && run a768_$k $P 1500; k=$((k+1)); done < ~/rule/p768.txt
k=0; while read P; do [ $k -le 2 ] && run b1030_$k $P 3600; k=$((k+1)); done < ~/rule/p1030.txt
k=0; while read P; do [ $k -ge 13 ] && run c768_$k $P 1500; k=$((k+1)); done < ~/rule/p768.txt
echo "$(date +%T) ALL_DONE" >> $OUT/progress
