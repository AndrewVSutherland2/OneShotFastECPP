#!/bin/bash
# the supply with two arms (wider bound / higher cap): gate chains, then the 768-bit set and three primes
# of about 1030 bits again, to compare with the runs before the change
cd ~/rule/repo
. ./setenv.sh > /dev/null
Y=~/yb/yafu/yafu
OUT=~/rule/out7
mkdir -p $OUT
run() {  # tag prime timeout
  rm -f work/shortcamp/$1.*
  /usr/bin/time -f "RUN $1 wall=%e user=%U sys=%S" timeout -k 20 $3 python3 ecpp/shortECPP.py p=$2 threads=192 v2=1 tag=$1 out=$OUT/$1.cert.txt sched=rule engine=yafu yafu=$Y > $OUT/$1.log 2>&1
  echo "EXIT $?" >> $OUT/$1.log
  pkill -x ecm; pkill -x dscan; pkill -x classpoly; pkill -x cm_method; pkill -x smoothtest; pkill -x gp
  echo "$(date +%T) $1 done: $(grep -o 'RUN.*' $OUT/$1.log) verified=$(grep -c 'verified twice' $OUT/$1.log)" >> $OUT/progress
  sleep 2
}
echo "$(date +%T) start on $(hostname)" >> $OUT/progress
while read C P; do case $C in 50|100|130|160|200) run c$C $P 1200;; esac; done < ~/rule/pall.txt
k=0; while read P; do [ $k -le 12 ] && run y768_$k $P 1500; k=$((k+1)); done < ~/rule/p768.txt
while read C P; do case $C in 300|310|320) run c$C $P 3600;; esac; done < ~/rule/pall.txt
echo "$(date +%T) ALL_DONE" >> $OUT/progress
