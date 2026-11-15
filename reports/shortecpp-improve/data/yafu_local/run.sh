#!/bin/bash
# yafu as the ECM engine on the 16-core driver box (32 threads): the 13 primes of 768 bits of section 9.3,
# then nextprime(10^310); one after the other
cd /home/claude/OneShotFastECPP
OUT=/tmp/claude-1000/-home-claude-OneShotFastECPP/f1db9bbb-8f84-4630-ad02-28f9e40f5dc3/scratchpad/ylocal
Y=/tmp/claude-1000/-home-claude-OneShotFastECPP/f1db9bbb-8f84-4630-ad02-28f9e40f5dc3/scratchpad/yafu/yafu/yafu
run() {  # tag prime timeout
  rm -f work/shortcamp/$1.*
  /usr/bin/time -f "RUN $1 wall=%e user=%U sys=%S" timeout -k 20 $3 python3 ecpp/shortECPP.py p=$2 threads=32 v2=1 tag=$1 out=$OUT/$1.cert.txt sched=rule engine=yafu yafu=$Y > $OUT/$1.log 2>&1
  echo "EXIT $?" >> $OUT/$1.log
  echo "$(date +%T) $1 done" >> $OUT/progress
}
echo "$(date +%T) start" >> $OUT/progress
k=0; while read P; do [ $k -le 12 ] && run y768_$k $P 5400; k=$((k+1)); done < /tmp/claude-1000/-home-claude-OneShotFastECPP/f1db9bbb-8f84-4630-ad02-28f9e40f5dc3/scratchpad/spot3/p768.txt
k=0; while read P; do [ $k -le 0 ] && run y1030_$k $P 14400; k=$((k+1)); done < /tmp/claude-1000/-home-claude-OneShotFastECPP/f1db9bbb-8f84-4630-ad02-28f9e40f5dc3/scratchpad/spot3/p1030.txt
echo "$(date +%T) ALL_DONE" >> $OUT/progress
