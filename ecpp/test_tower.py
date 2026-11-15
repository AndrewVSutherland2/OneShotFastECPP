#!/usr/bin/env python3
"""Validate the class-group decomposition (Algorithm 1 of arXiv:1009.1082):
classpoly's CLASSPOLY_TOWER output and cm_method's use of it.

  A. Over Z.  For a list of (D, invariant, subgroup) the decomposition (V, W_k)
     is compared with the plain class polynomial H computed by the same binary:
     modulo three primes q that split completely in the ring class field, V must
     have m distinct roots y, every U_y = X^n + y X^(n-1) + sum W_k(y)/V'(y) X^k
     must divide H, and their product must be H (PARI).  Covers j and other
     invariants, square roots of class polynomials (ramified level), subgroups of
     the paper's prefix form and of the wider form <a_i^c_i>, sign-free
     decompositions (Weber f of even degree), and the height bound.
  B. Mod P (explicit CRT), one job and several: the two outputs must be
     identical, and pass the same check against H mod P.
  C. A decomposition whose V has a repeated root (s = e_1 fails): classpoly must
     notice and write the plain class polynomial instead.
  D. cm_method: with and without the decomposition, and with a damaged
     decomposition (fault=1,2: the retry with a larger height bound must succeed;
     fault=3,4: must fall back to the plain polynomial), the j it returns must be a root
     of the Hilbert class polynomial mod p (PARI polclass), or for large |D|
     the j-invariant of a curve with the right trace (PARI ellap).

Needs the classpoly environment (. ./setenv.sh) and gp.  One to two minutes.
"""
import math
import os
import random
import re
import subprocess
import sys
import tempfile

sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__))
CM = os.path.join(HERE, "cm_method")
FAILS = []


def check(ok, what):
    print("  %-4s %s" % ("ok" if ok else "FAIL", what), flush=True)
    if not ok:
        FAILS.append(what)


def classpoly(D, inv, P, out, env=None, verb=0):
    e = dict(os.environ)
    for k in ("CLASSPOLY_TOWER", "CLASSPOLY_TOWER_N", "CLASSPOLY_TOWER_C", "CLASSPOLY_JOBS", "CLASSPOLY_JOBID"):
        e.pop(k, None)
    e.update(env or {})
    r = subprocess.run(["classpoly", str(D), str(inv), str(P), out, str(verb)], capture_output=True, text=True, env=e,
                       cwd=os.path.dirname(out))
    return r.stdout + r.stderr


def load(fn):
    """classpoly file -> (header dict, coefficient list in file order)"""
    hdr, co = {}, []
    for ln in open(fn):
        ln = ln.strip()
        if not ln:
            continue
        if ln[0].isalpha():
            k, v = ln.split("=", 1)
            hdr[k] = v
        else:
            co.append(int(ln.split("*")[0]))
    return hdr, co


# chk(q): V must have m distinct roots mod q, each U_y must divide the class polynomial and their product must equal it.
# For a sign-free decomposition (sf=1) the coefficients marked odd are first divided by a square root of S, and the
# result may be the decomposition of H(X) or of H(-X).
GP_CHK = ("chk(q)={my(r,Vq,Hq,H2,pr,ok,ok2,U,d,y,ti,Wq);ti=Mod(1,q);"
          "if(sf,if(!issquare(Mod(S,q))||S%q==0,return(0));ti=1/sqrt(Mod(S,q)));"
          "Vq=sum(j=0,m,V[j+1]*if(PV[j+1],ti,Mod(1,q))*'y^j);"
          "Wq=vector(n-1,k,sum(j=0,m-1,W[k][j+1]*if(PW[k][j+1],ti,Mod(1,q))*'y^j));"
          "Hq=H*Mod(1,q);H2=subst(Hq,'x,-'x)*(-1)^(m*n);r=polrootsmod(lift(Vq),q);if(#r!=m,return(0));"
          "pr=Mod(1,q);ok=1;ok2=1;for(i=1,m,y=r[i];d=subst(deriv(Vq),'y,y);if(d==0,return(0));"
          "U='x^n+y*'x^(n-1)+sum(k=0,n-2,subst(Wq[k+1],'y,y)/d*'x^k);"
          "if(Hq%U!=0,ok=0);if(H2%U!=0,ok2=0);pr*=U);(ok&&pr==Hq)||(sf&&ok2&&pr==H2)};")


def gp(script):
    r = subprocess.run(["gp", "-q", "--default", "parisizemax=4G"], input=script, capture_output=True, text=True)
    return r.stdout


def decomposition_ok(tower_file, plain_file):
    """True iff the decomposition in tower_file is consistent with the class polynomial in plain_file."""
    th, tc = load(tower_file)
    ph, pc = load(plain_file)
    if "F" not in th or "G" not in th:
        return False
    m, n = int(th["F"]), int(th["G"])
    if len(tc) != m * n + 1 or tc[m] != 1 or len(pc) != m * n + 1 or pc[-1] != 1:
        return False
    if th["I"] != ph["I"] or th["D"] != ph["D"] or th.get("P", "0") != ph.get("P", "0"):
        return False
    P = int(th.get("P", "0"))
    V = tc[:m + 1]
    W = [tc[m + 1 + k * m: m + 1 + (k + 1) * m] for k in range(n - 1)]
    # parities under f -> -f: coefficient of Y^j in V is odd iff m-j is, that of Y^j in W_k iff n-k+m-1-j is
    PV = [(m - j) & 1 for j in range(m + 1)]
    PW = [[(n - k + m - 1 - j) & 1 for j in range(m)] for k in range(n - 1)]
    s = ["V=%s;" % V, "H=Pol(%s,'x);" % pc[::-1], "W=%s;" % W, "PV=%s;" % PV, "PW=%s;" % PW,
         "n=%d; m=%d; D=%s; sf=%d; S=%s;" % (n, m, th["D"], 1 if "S" in th else 0, th.get("S", "0")), GP_CHK]
    if P:
        s.append("print(chk(%d));" % P)
        want = ["1"]
    else:
        # primes q with 4q = t^2 + |D| v^2 split completely in the ring class field
        s.append("qs=List();for(t=20001,10^9,for(v=1,2,q=(t^2-D*v^2)/4;if(denominator(q)==1&&isprime(q),listput(qs,q)));if(#qs>=3,break));")
        s.append("for(i=1,3,print(chk(qs[i])));")
        want = ["1", "1", "1"]
    got = [t for t in gp("\n".join(s) + "\n").split() if t in ("0", "1")]
    return got == want


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def norm_prime(D, bits, seed):
    """a prime p of the given size with 4p = t^2 + |D| v^2"""
    rng = random.Random(seed)
    lo, hi = math.isqrt(1 << (bits + 1)) + 1, math.isqrt(1 << (bits + 2)) - (1 << 40)     # t^2/4 has exactly `bits` bits
    while True:
        t = rng.randrange(lo, hi)
        for v in (1, 2):
            if (t * t - D * v * v) % 4 == 0:
                p = (t * t - D * v * v) // 4
                if p.bit_length() == bits and is_prime(p):
                    return p


# (D, invariant, subgroup): the vector c defining <a_i^c_i>, an order "n=<order>", or None for the automatic choice
CASES_Z = [
    (-971, 0, "n=5"),               # the worked example of the paper: n = 5, m = 3 with j
    (-23447, 0, None),              # j, cyclic of order 123
    (-23447, 1, "n=41"),            # Weber f
    (-1000003, 171, "n=15"),        # Atkin A_71
    (-200004, 26, "2,2,2"),         # w_{2,13}, three generators, not of prefix form
    (-200004, 26, "1,8,2"),
    (-4000004, 10, None),           # w_{2,5}, degree 1032
    (-200008, 123, "1,6"),          # A_23 with 23 ramified: square root of the class polynomial
    (-300000003, 117, "1,4,3,2"),   # A_17 ramified, four generators, degree 1536
    (-3000008, 23, None),
    (-100319, 1, "sf"),             # Weber f of even degree: sign-free (CLASSPOLY_TOWER=2), m = 10 and n = 38
    (-200447, 1, "sf"),             # the same with m odd
    (-200255, 1, "sf"),
]


def tower_env(sub):
    if sub is None:
        return {"CLASSPOLY_TOWER": "1"}
    if sub == "sf":
        return {"CLASSPOLY_TOWER": "2"}
    return {"CLASSPOLY_TOWER_N": sub[2:]} if sub.startswith("n=") else {"CLASSPOLY_TOWER_C": sub}


def test_over_Z(tmp):
    print("A. decomposition over Z against the plain class polynomial")
    for D, inv, c in CASES_Z:
        plain = os.path.join(tmp, "p.txt")
        tower = os.path.join(tmp, "t.txt")
        for f in (plain, tower, tower + ".bad"):
            if os.path.exists(f):
                os.remove(f)
        classpoly(D, inv, 0, plain)
        o = classpoly(D, inv, 0, tower, tower_env(c), verb=1)
        mt = re.search(r"Decomposition: F of degree (\d+) \(cosets\), G of degree (\d+)", o)
        mh = re.search(r"height (\d+) \((\d+)\)", o)
        if not (os.path.exists(tower) and mt and mh):
            check(False, "D=%d inv=%d %s: no decomposition written" % (D, inv, c))
            continue
        ok = decomposition_ok(tower, plain) and (("S" in load(tower)[0]) == (c == "sf"))
        check(ok and int(mh.group(1)) <= int(mh.group(2)),
              "D=%d inv=%d subgroup %s: m=%s n=%s, height %s <= bound %s" % (D, inv, c or "auto", mt.group(1), mt.group(2), mh.group(1), mh.group(2)))


def test_mod_P(tmp):
    print("B. decomposition mod P: one job, four jobs, and the plain polynomial mod P")
    for D, inv, bits, tw in ((-1000003, 171, 256, "1"), (-300000003, 117, 384, "1"), (-4000004, 10, 521, "1"), (-200255, 1, 384, "2")):
        P = norm_prime(D, bits, 7)
        plain = os.path.join(tmp, "pp.txt")
        one = os.path.join(tmp, "t1.txt")
        four = os.path.join(tmp, "t4.txt")
        for f in (plain, one, four):
            if os.path.exists(f):
                os.remove(f)
        classpoly(D, inv, P, plain)
        classpoly(D, inv, P, one, {"CLASSPOLY_TOWER": tw})
        # the automatic choice depends on the number of jobs, so pin the subgroup order to the one-job choice
        sub = {"CLASSPOLY_TOWER": tw, "CLASSPOLY_TOWER_N": load(one)[0].get("G", "0")}
        edir = tempfile.mkdtemp(dir=tmp)
        procs = []
        for j in (1, 2, 3, 4):
            e = dict(os.environ)
            e.update(sub)
            e.update({"CLASSPOLY_JOBS": "4", "CLASSPOLY_JOBID": str(j), "CLASSPOLY_ECRT_DIR": edir})
            procs.append(subprocess.Popen(["classpoly", str(D), str(inv), str(P), four + ".w%d" % j, "0"], env=e, cwd=tmp,
                                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL))
        rcs = [p.wait() for p in procs]
        classpoly(D, inv, P, four, dict(sub, CLASSPOLY_JOBS="4", CLASSPOLY_JOBID="0", CLASSPOLY_ECRT_DIR=edir))
        same = os.path.exists(one) and os.path.exists(four) and open(one).read() == open(four).read()
        check(rcs == [0, 0, 0, 0] and same, "D=%d inv=%d, %d-bit P: four-job output identical to one-job output" % (D, inv, bits))
        check(os.path.exists(one) and decomposition_ok(one, plain) and (("S" in load(one)[0]) == (tw == "2")),
              "D=%d inv=%d, %d-bit P: decomposition%s consistent with H mod P" % (D, inv, bits, " (sign-free)" if tw == "2" else ""))


def test_repeated_root(tmp):
    print("C. a decomposition with a repeated coset value reverts to the plain class polynomial")
    D, inv, n = -1603, 171, 3
    plain = os.path.join(tmp, "rp.txt")
    tower = os.path.join(tmp, "rt.txt")
    classpoly(D, inv, 0, plain)
    o = classpoly(D, inv, 0, tower, {"CLASSPOLY_TOWER_N": str(n)})
    check("repeated coset value" in o, "D=%d inv=%d n=%d: repeated value detected" % (D, inv, n))
    check(os.path.exists(tower) and open(tower).read() == open(plain).read(), "D=%d inv=%d n=%d: output is the plain class polynomial" % (D, inv, n))


def test_cm_method():
    print("D. cm_method: the j-invariant is a root of the Hilbert class polynomial mod p")
    env = dict(os.environ)
    for D, bits in ((-23447, 256), (-1000003, 256), (-200004, 384), (-100319, 320), (-300000003, 256)):
        p = norm_prime(D, bits, 11)
        js = {}
        for name, args in (("decomposed", []), ("plain", ["tower=0"]), ("fault=1", ["fault=1"]), ("fault=2", ["fault=2"]),
                           ("fault=3", ["fault=3"]), ("fault=4", ["fault=4"])):
            r = subprocess.run([CM, "D=%d" % D, "p=%d" % p] + args, capture_output=True, text=True, env=env)
            mj = re.search(r"^j (\d+)$", r.stdout, re.M)
            js[name] = (int(mj.group(1)) if mj else None, r.stderr)
        used = " x " in js["decomposed"][1]
        # a damaged first attempt is retried with a larger height bound; if the retry is damaged too, the plain polynomial is used
        fb = (all("second attempt" in js[k][1] for k in ("fault=1", "fault=2"))
              and all("(after a failed decomposition)" in js[k][1] for k in ("fault=3", "fault=4")))
        vals = [js[k][0] for k in js]
        if None in vals:
            check(False, "D=%d: cm_method returned no j (%s)" % (D, [k for k in js if js[k][0] is None]))
            continue
        if -D <= 2000000:
            # root of the Hilbert class polynomial
            out = gp("default(parisizemax,\"4G\");\nH=polclass(%d);\np=%d;\n" % (D, p) + "".join("print(subst(H,x,Mod(%d,p))==0);\n" % v for v in vals))
            how = "roots of H_D"
        else:
            # H_D is too big for a quick test: check the trace of the curve instead (t^2 - 4p = v^2 D), with PARI's point count
            out = gp("p=%d;\nD=%d;\n" % (p, D) + "".join(
                "j=Mod(%d,p);\nk=j/(1728-j);\nt=ellap(ellinit([3*k,2*k]));\nprint(issquare((t^2-4*p)/D));\n" % v for v in vals))
            how = "curves whose trace t has (t^2-4p)/D a square"
        good = out.split() == ["1"] * len(vals)
        check(good and used and fb, "D=%d, %d-bit p: decomposed, plain, retried and fallback runs give %s (decomposition used: %s, fell back: %s)"
              % (D, bits, how, used, fb))


def main():
    if not os.environ.get("CLASSPOLY_PHI_DIR"):
        print("source setenv.sh first")
        return 2
    with tempfile.TemporaryDirectory(prefix="tower_test_") as tmp:
        test_over_Z(tmp)
        test_mod_P(tmp)
        test_repeated_root(tmp)
    test_cm_method()
    print("test_tower: %s" % ("ALL OK" if not FAILS else "%d FAILURES" % len(FAILS)))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
