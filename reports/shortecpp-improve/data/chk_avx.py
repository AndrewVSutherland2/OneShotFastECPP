"""Check avx-ecm stage-1 residues (its save_b1.txt) against an independent
stage 1 modulo a known prime factor: Suyama curve from sigma, the same prime
powers avx-ecm uses, plain Montgomery ladder.
usage: chk_avx.py <save_b1.txt> <file holding the prime factor> <B1>"""
import sys
def xdbl(X,Z,a24,p):
    s=(X+Z)%p; d=(X-Z)%p; s2=s*s%p; d2=d*d%p; t=(s2-d2)%p
    return s2*d2%p, t*(d2+a24*t)%p
def xadd(X1,Z1,X2,Z2,Xd,Zd,p):
    a=(X1-Z1)*(X2+Z2)%p; b=(X1+Z1)*(X2-Z2)%p
    return Zd*(a+b)**2%p, Xd*(a-b)**2%p
def ladder(k,X,Z,a24,p):
    R0=(X,Z); R1=xdbl(X,Z,a24,p)
    for bit in bin(k)[3:]:
        if bit=='1': R0=xadd(R0[0],R0[1],R1[0],R1[1],X,Z,p); R1=xdbl(R1[0],R1[1],a24,p)
        else: R1=xadd(R0[0],R0[1],R1[0],R1[1],X,Z,p); R0=xdbl(R0[0],R0[1],a24,p)
    return R0
def primes(n):
    s=bytearray([1])*(n+1); s[0:2]=b'\0\0'
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
def check(savefile, p, B1):
    ok=n=0
    for line in open(savefile):
        d=dict(kv.strip().split('=',1) for kv in line.strip().strip(';').split(';') if '=' in kv)
        sig=int(d['SIGMA']); X=int(d['X'],16); Z=int(d['Z'],16) if 'Z' in d else 1
        u=(sig*sig-5)%p; v=4*sig%p
        x0=pow(u,3,p)*pow(pow(v,3,p),-1,p)%p
        A=((v-u)**3*(3*u+v)*pow(4*u**3*v,-1,p)-2)%p
        a24=(A+2)*pow(4,-1,p)%p
        P=(x0,1)
        q=2
        while q<B1: P=xdbl(P[0],P[1],a24,p); q*=2
        for q in primes(B1-1)[1:]:
            c=q
            P=ladder(q,P[0],P[1],a24,p)
            while c*q<B1: P=ladder(q,P[0],P[1],a24,p); c*=q
        want = P[0]*pow(P[1],-1,p)%p if P[1]%p else None
        got = X*pow(Z,-1,p)%p if Z%p else None
        n+=1; ok+= (want==got)
    return ok,n
if __name__=="__main__":
    p=int(open(sys.argv[2]).read()); print("residues matching an independent stage 1: %d of %d"%check(sys.argv[1],p,int(sys.argv[3])))
