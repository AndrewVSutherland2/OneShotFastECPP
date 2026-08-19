#!/bin/bash
# drive.sh <tag> <prime-expr> <workers>  — first worker to find a cert wins; others killed by PID
tag=$1; expr=$2; W=$3
dir=run.$tag; mkdir -p $dir; : > $dir/pids
for i in $(seq 1 $W); do
  (echo "setrand($i); printrestricted($expr)" | gp -q restricted.gp > $dir/out.$i.txt 2>/dev/null) &
  echo $! >> $dir/pids
done
while :; do
  sleep 10
  if grep -h '^CERT' $dir/out.*.txt 2>/dev/null | head -1 | grep -q CERT; then
    for pid in $(cat $dir/pids); do kill $pid 2>/dev/null; done
    sleep 1
    grep -h '^CERT' $dir/out.*.txt | head -1 > $dir/RESULT
    grep -h '^CURVES' $dir/out.*.txt | head -1 >> $dir/RESULT
    awk '/^PROG/{if($2>m[FILENAME])m[FILENAME]=$2} END{t=0;for(f in m)t+=m[f];print "TOTALPROG",t}' $dir/out.*.txt >> $dir/RESULT
    date >> $dir/RESULT
    exit 0
  fi
  if ! kill -0 $(head -1 $dir/pids) 2>/dev/null && [ $(cat $dir/pids | while read p; do kill -0 $p 2>/dev/null && echo live; done | wc -l) -eq 0 ]; then
    echo "ALL WORKERS DIED" > $dir/RESULT; exit 1
  fi
done
