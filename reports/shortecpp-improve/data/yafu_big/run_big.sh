#!/bin/bash
# run_big.sh C: nextprime(10^C) on this box.  resume=1 takes up whatever is in work/shortcamp
# for the tag: the levels done, the winners in hand and the search checkpoint.
C=$1
cd ~/rule/repo
. ./setenv.sh > /dev/null
Y=~/yb/yafu/yafu
OUT=~/rule/big
mkdir -p $OUT
P=$(awk -v c=$C '$1==c {print $2}' ~/rule/pbig.txt)
echo "$(date -u +%F_%T) start c$C on $(hostname) ($(lscpu | grep 'Model name' | sed 's/.*: *//'))" >> $OUT/progress
export SHORTECPP_CKPT=900
/usr/bin/time -f "RUN c$C wall=%e user=%U sys=%S" timeout -k 20 432000 python3 ecpp/shortECPP.py p=$P threads=192 v2=1 tag=c$C out=$OUT/c$C.cert.txt sched=rule engine=yafu yafu=$Y resume=1 >> $OUT/c$C.log 2>&1
echo "EXIT $?" >> $OUT/c$C.log
echo "$(date -u +%F_%T) c$C done: $(grep -o 'RUN.*' $OUT/c$C.log | tail -1) verified=$(grep -c 'verified twice' $OUT/c$C.log)" >> $OUT/progress
echo "$(date -u +%F_%T) ALL_DONE" >> $OUT/progress
