#!/bin/bash
# Candidate supply for the policy experiment: scan, strip, factor the smooth parts.
R=$HOME/OneShotFastECPP; E=$HOME/exp768; mkdir -p $E; cd $E
P=$(cat p768.txt)
{
date
/usr/bin/time -f "dscan: wall=%e user=%U" $R/ecpp/dscan p=$P B=200000000000 maxfb=100000000 threads=192 dump > scan.txt
wc -l scan.txt
$R/ecpp/smoothtest pbuild y=67108864 save=$E/P26.bin threads=64 | tail -1
python3 - <<'PY'
p = int(open("p768.txt").read())
with open("orders.txt", "w") as f, open("orders_d.txt", "w") as g:
    seen = set()
    for line in open("scan.txt"):
        d, t, v = map(int, line.split())
        if d <= 4:
            continue
        for N in (p + 1 - t, p + 1 + t):
            if N % 4 == 0 and N not in seen:
                seen.add(N)
                f.write("%d\n" % N)
                g.write("%d %d\n" % (d, v))
PY
wc -l orders.txt
/usr/bin/time -f "strip: wall=%e user=%U" $R/ecpp/smoothtest parts load=$E/P26.bin threads=192 n2=589824 prp=1 < orders.txt > strip.txt
awk '{print $4}' strip.txt > s.txt
cat > fac.gp <<'GP'
v=readvec("s.txt");
{for(i=1,#v, f=factor(v[i]); s=""; for(j=1,#f~, s=Str(s,f[j,1],"^",f[j,2]," ")); print(s))}
GP
/usr/bin/time -f "factor s: wall=%e" gp -q < fac.gp > strip.txt.fac
wc -l strip.txt strip.txt.fac
# compact file for the model: tail bits, prp flag
python3 - <<'PY'
import math
with open("compact.txt", "w") as f:
    for line in open("strip.txt"):
        w = line.split()
        T = int(w[2])
        f.write("%.4f %s %s\n" % (math.log2(T) if T > 1 else 0.0, w[4], w[3]))
PY
date
echo PREP_DONE
} > prep.log 2>&1
