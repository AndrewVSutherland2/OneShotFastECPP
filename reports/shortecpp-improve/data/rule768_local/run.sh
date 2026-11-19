#!/bin/bash
# the rule scheduler on fresh 768-bit primes, one after the other, all local threads
cd /home/claude/OneShotFastECPP
k=0
while read P; do
  k=$((k+1))
  rm -f work/shortcamp/rb768_$k.*
  /usr/bin/time -f "RUN $k wall=%e user=%U sys=%S" timeout 5400 python3 ecpp/shortECPP.py p=$P threads=32 v2=1 tag=rb768_$k out=$OUT/cert_$k.txt sched=rule > $OUT/p$k.log 2>&1
  echo "EXIT $?" >> $OUT/p$k.log
done < $OUT/primes.txt
echo BATCH_DONE > $OUT/done
