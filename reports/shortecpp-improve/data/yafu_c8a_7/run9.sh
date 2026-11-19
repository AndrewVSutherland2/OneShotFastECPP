#!/bin/bash
# on the cores the 10^410 build leaves idle: the 768-bit prime whose winner once needed a 90-second build,
# twelve times over, to see the reserve step act (the timings are not clean: the box is shared with the build)
cd ~/rule/repo
. ./setenv.sh > /dev/null
Y=~/yb/yafu/yafu
OUT=~/rule/out9
mkdir -p $OUT
P=$(sed -n 1p ~/rule/p768.txt)
for k in 1 2 3 4 5 6 7 8 9 10 11 12; do
  tag=r768_$k
  rm -f work/shortcamp/$tag.*
  /usr/bin/time -f "RUN $tag wall=%e user=%U sys=%S" timeout -k 20 1500 python3 ecpp/shortECPP.py p=$P threads=150 v2=1 tag=$tag out=$OUT/$tag.cert.txt sched=rule engine=yafu yafu=$Y > $OUT/$tag.log 2>&1
  echo "$(date +%T) $tag done: $(grep -o 'RUN.*' $OUT/$tag.log) verified=$(grep -c 'verified twice' $OUT/$tag.log) reserve=$(grep -c 'kept in reserve' $OUT/$tag.log) D=$(grep '769 bits: built' $OUT/$tag.log | sed 's/.*built D=-\([0-9]*\).* in \([0-9.]*\)s.*/\1 build \2s/')" >> $OUT/progress
done
echo "$(date +%T) ALL_DONE" >> $OUT/progress
