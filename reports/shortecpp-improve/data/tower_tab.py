import re, sys
rows=[]; cur=None; retry=0
for ln in open(sys.argv[1]):
    if 'retrying' in ln: retry=1
    m=re.match(r"\[H_D\^inv\(deg (\d+) = (\d+) x (\d+), inv (\d+)\): ([\d.]+)s\s+roots: ([\d.]+)s \+ ([\d.]+)s",ln)
    if m: cur=tuple(m.groups())
    m2=re.match(r"\[H_D\^inv\(deg (\d+), inv (\d+)\): ([\d.]+)s\s+root: ([\d.]+)s",ln)
    if m2: cur=(m2.group(1),'-','-',m2.group(2),m2.group(3),m2.group(4),'0')
    m=re.match(r"D=(\d+) wall=([\d.]+) user=([\d.]+)",ln)
    if m and float(m.group(2))>0:
        rows.append((int(m.group(1)),)+cur+(float(m.group(2)),float(m.group(3)),retry)); retry=0
print("%14s %7s %11s %4s %7s %7s %7s %8s %9s"%("|D|","deg","m x n","inv","crt s","root m","root n","wall s","core-s"))
tw=tu=0
for r in rows:
    print("%14d %7s %11s %4s %7.1f %7.1f %7.1f %8.1f %9.0f %s"%(r[0],r[1],"%s x %s"%(r[2],r[3]),r[4],float(r[5]),float(r[6]),float(r[7]),r[8],r[9],"retry" if r[10] else ""))
    tw+=r[8]; tu+=r[9]
print("%d constructions: total wall %.0f s (median %.1f), total %.0f core-s"%(len(rows),tw,sorted(r[8] for r in rows)[len(rows)//2],tu))
