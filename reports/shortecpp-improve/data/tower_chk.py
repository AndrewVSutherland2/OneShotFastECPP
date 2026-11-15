#!/usr/bin/env python3
"""Check a decomposed class polynomial file (F/G format) against the plain class polynomial.

usage: chk.py tower.txt plain.txt [q ...]
For each test prime q that splits completely (found automatically if none given): V mod q must have m distinct roots,
each U_y must divide H mod q, and the product of the U_y must equal H mod q.  Uses PARI.
"""
import subprocess, sys
sys.set_int_max_str_digits(0)

def load(fn):
    hdr, co = {}, []
    for ln in open(fn):
        ln = ln.strip()
        if not ln: continue
        if ln[0].isalpha():
            k, v = ln.split('=', 1); hdr[k] = v; continue
        co.append(int(ln.split('*')[0]))
    return hdr, co

th, tc = load(sys.argv[1]); ph, pc = load(sys.argv[2])
m, n = int(th['F']), int(th['G'])
assert len(tc) == m*n + 1 and tc[m] == 1, (len(tc), m, n)
assert len(pc) == m*n + 1 and pc[-1] == 1, (len(pc), m*n)
assert th['I'] == ph['I'] and th['D'] == ph['D']
P = int(th.get('P', '0'))
assert P == int(ph.get('P', '0'))
V = tc[:m+1]
W = [tc[m+1+k*m: m+1+(k+1)*m] for k in range(n-1)]
qs = [int(a) for a in sys.argv[3:]]
gp = []
gp.append("V=Pol(%s,'y);" % str(V[::-1]))
gp.append("H=Pol(%s,'x);" % str(pc[::-1]))
gp.append("W=%s;" % str(W))
gp.append("n=%d; m=%d; D=%s;" % (n, m, th['D']))
gp.append("chk(q)={my(r,Vq,Hq,pr,ok,U,d,y);Vq=V*Mod(1,q);Hq=H*Mod(1,q);r=polrootsmod(V,q);if(#r!=m,return(0));pr=Mod(1,q);ok=1;for(i=1,m,y=r[i];d=subst(deriv(Vq),'y,y);if(d==0,return(0));U='x^n+y*'x^(n-1)+sum(k=0,n-2,subst(Pol(Vecrev(W[k+1]),'y)*Mod(1,q),'y,y)/d*'x^k);if(Hq%U!=0,ok=0);pr*=U);ok&&(pr==Hq)};")
if P:
    gp.append("print(chk(%d));" % P)
else:
    if not qs:
        # primes q with 4q = t^2 + |D| split completely in the ring class field, so H splits completely mod q
        gp.append("qs=List();for(t=20001,10^9,for(v=1,2,q=(t^2-D*v^2)/4;if(denominator(q)==1&&isprime(q),listput(qs,q)));if(#qs>=3,break));")
        gp.append("for(i=1,#qs,print(qs[i],\" \",chk(qs[i])));")
    else:
        gp.append("".join("print(%d,\" \",chk(%d));" % (q, q) for q in qs))
r = subprocess.run(["gp", "-q", "--default", "parisizemax=8G"], input="\n".join(gp) + "\n", capture_output=True, text=True)
print(r.stdout.strip(), r.stderr.strip()[:300])
