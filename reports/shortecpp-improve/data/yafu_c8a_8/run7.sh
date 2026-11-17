#!/bin/bash
# nextprime(10^410) again, with the supply that raises the factor-base cap and the longer 3e7 ladder
cd ~/rule/repo
. ./setenv.sh > /dev/null
Y=~/yb/yafu/yafu
OUT=~/rule/out8
mkdir -p $OUT
P=$(awk '$1==410 {print $2}' ~/rule/pbig.txt)
rm -f work/shortcamp/c410.* work/rootcache/*
echo "$(date +%T) start c410 on $(hostname)" >> $OUT/progress
/usr/bin/time -f "RUN c410 wall=%e user=%U sys=%S" timeout -k 20 259200 python3 ecpp/shortECPP.py p=$P threads=192 v2=1 tag=c410 out=$OUT/c410.cert.txt sched=rule engine=yafu yafu=$Y resume=1 > $OUT/c410.log 2>&1
echo "EXIT $?" >> $OUT/c410.log
echo "$(date +%T) c410 done: $(grep -o 'RUN.*' $OUT/c410.log) verified=$(grep -c 'verified twice' $OUT/c410.log)" >> $OUT/progress
echo "$(date +%T) ALL_DONE" >> $OUT/progress
