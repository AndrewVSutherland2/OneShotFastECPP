#!/usr/bin/env python3
"""test_vshort2.py -- differential test of ecpp/vshort2 (C; faithful and --fast modes)
against vshort2.py.

Cases: vshort2.py's invalid vectors (incl. the CRT split attack); every chain in
certs/short2/certs.csv; random multi-level chains built with PARI (fixed n, B and radical
cap from p_0; non-terminal levels with B^2 < p_next < sqrt(p) prime, a terminal level
with p_next = 1 or < B^2); CRT forgeries over composite p = l1*l2 (one level whose
point has order o over Z/pZ, split between the two primes -- all must be False); and
per-level tampering of all of these (A, x, o; x -> x(rP) for small r | o, which only the
exactness tree catches; dropped / duplicated / swapped levels; p_0 + 2).  The three
verifiers must agree on every case.

    make -C ecpp vshort2 && python3 ecpp/test_vshort2.py [nchains] [seed]
"""
import csv, os, subprocess, sys
from math import ceil, gcd, log2

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "ecpp"))
import vshort2

CBIN = os.path.join(ROOT, "ecpp", "vshort2")

GP = r"""
default(parisizemax, "1G");
lev(p, n, B, radlim, want, P0) = {
  my(B2 = B^2, sp = sqrtint(p), L = sp + 1 + sqrtint(4*sp));
  for (tries = 1, 300,
    my(A = random(p));
    if (Mod(A, p)^2 == 4, next);
    my(E = ellinit([0, A, 0, 1, 0], p), P = random(E));
    if (P == [0], next);
    my(o = ellorder(E, P), cand = List());
    fordiv(o, d,
      if (d <= L, next);
      my(f = factor(d), m = 1, rad = 1);
      for (j = 1, #f~, if (f[j,1] <= B, m *= f[j,1]^f[j,2]; rad *= f[j,1]));
      if (m == 1 || #binary(rad) > radlim, next);
      if (d >= vecmin(f[,1]) * L, next);
      my(pn = d / m);
      if (want, if (pn > B2 && pn^2 < p && isprime(pn), listput(cand, [d, pn])),
                if (pn == 1 || pn < B2, listput(cand, [d, pn]))));
    if (#cand == 0, next);
    my(c = cand[random(#cand) + 1], Q = ellmul(E, P, o / c[1]));
    return([A, lift(Q[1]), c[1], c[2]]));
  0
}
chain(bits) = {
  for (t = 1, 40,
    my(p0 = randomprime([2^(bits-1), 2^bits - 1]), n = #binary(p0), lg = log(n)/log(2));
    my(B = ceil(n^2/lg), radlim = ceil(n/lg), p = p0, out = [p0], ok = 1);
    while (1,
      my(nonterm = (p > B^4), v = lev(p, n, B, radlim, nonterm, p0));
      if (v == 0 && nonterm, nonterm = 0; v = lev(p, n, B, radlim, 0, p0));
      if (v == 0, ok = 0; break);
      out = concat(out, [v[1], v[2], v[3]]);
      if (!nonterm, break);
      p = v[4]);
    if (ok, return(out)));
  0
}
forge(bits) = {
  for (tries = 1, 3000,
    my(l1 = randomprime([2^(bits\2 - 1), 2^(bits\2)]), l2 = randomprime([2^(bits - bits\2 - 1), 2^(bits - bits\2)]));
    if (l1 == l2, next);
    my(p = l1*l2, n = #binary(p), lg = log(n)/log(2), B = ceil(n^2/lg), radlim = ceil(n/lg));
    my(A = random(p), x0 = random(p), y1, y2);
    if (gcd(A^2 - 4, p) != 1, next);
    my(rhs = x0^3 + A*x0^2 + x0);
    if (!issquare(Mod(rhs, l1), &y1) || !issquare(Mod(rhs, l2), &y2), next);
    my(E1 = ellinit([0, A, 0, 1, 0], l1), E2 = ellinit([0, A, 0, 1, 0], l2));
    my(P1 = [Mod(x0, l1), y1], P2 = [Mod(x0, l2), y2]);
    my(o = lcm(ellorder(E1, P1), ellorder(E2, P2)), sp = sqrtint(p), L = sp + 1 + sqrtint(4*sp), cand = List());
    fordiv(o, d,
      if (d <= L, next);
      my(f = factor(d), m = 1, rad = 1);
      for (j = 1, #f~, if (f[j,1] <= B, m *= f[j,1]^f[j,2]; rad *= f[j,1]));
      if (m == 1 || #binary(rad) > radlim || d >= vecmin(f[,1]) * L, next);
      my(pn = d / m);
      if (pn == 1 || pn < B^2, listput(cand, d)));
    if (#cand == 0, next);
    my(d = cand[random(#cand) + 1], Q1 = ellmul(E1, P1, o / d), Q2 = ellmul(E2, P2, o / d));
    if (Q1 == [0] || Q2 == [0], next);
    return([p, A, lift(chinese(Q1[1], Q2[1])), d]));
  0
}
"""


def gp_generate(nchains, seed):
    lines = [f"setrand({seed});"]
    for i in range(nchains):
        lines.append(f'v = chain({24 + (i * 7) % 77}); if (v, print("V ", v));')
    for i in range(nchains // 2):
        lines.append(f'v = forge({24 + (i * 5) % 37}); if (v, print("C ", v));')
    out = subprocess.run(["gp", "-q", "-D", "colors=no"], input=GP + "\n".join(lines) + "\nquit\n",
                         capture_output=True, text=True, check=True).stdout
    valid, forged = [], []
    for ln in out.splitlines():
        tag, _, vec = ln.partition(" ")
        if tag in ("V", "C"):               # gp also echoes the function definitions
            vals = [int(t) for t in vec.strip(" []").replace(",", " ").split()]
            (valid if tag == "V" else forged).append(vals)
    return valid, forged


def multiple_x(k, x, A, p):
    X, Z = vshort2.ladder(k, x, 1, A, p)
    return X * pow(Z, -1, p) % p if gcd(Z, p) == 1 else None


def moduli(seq):
    """p_0, p_1, ... as the verifier recovers them from the o_i."""
    n = seq[0].bit_length()
    primes = vshort2.sieve_primes(ceil(n * n / log2(n)))
    P = vshort2.balanced_product(primes)
    ps, p = [], seq[0]
    for i in range(1, len(seq), 3):
        ps.append(p)
        o = seq[i + 2]
        g, m = gcd(P % o, o), 1
        for q in primes:
            if q > g:
                break
            if g % q == 0:
                while (o // m) % q == 0:
                    m *= q
        p = o // m
    return ps


def tampers(seq):
    out = [("p0+2", [seq[0] + 2] + seq[1:])]
    k = (len(seq) - 1) // 3
    if k >= 2:
        out += [("drop last level", seq[:-3]), ("swap levels 0,1", seq[:1] + seq[4:7] + seq[1:4] + seq[7:])]
    out.append(("dup last level", seq + seq[-3:]))
    ps = moduli(seq)
    for lev in range(k):
        i = 1 + 3 * lev
        p = ps[lev] if lev < len(ps) else seq[0]
        A, x, o = seq[i:i + 3]
        for name, j, v in (("A+1", i, (A + 1) % p), ("x+1", i + 1, (x + 1) % p), ("o+1", i + 2, o + 1),
                           ("o-1", i + 2, o - 1), ("2o", i + 2, 2 * o), ("x=p", i + 1, p)):
            s = list(seq)
            s[j] = v
            out.append((f"L{lev} {name}", s))
        for r in [r for r in (2, 3, 5, 7, 11, 13) if o % r == 0][:2]:
            xr = multiple_x(r, x, A, p)
            if xr is not None:
                s = list(seq)
                s[i + 1] = xr
                out.append((f"L{lev} x(rP)", s))
    return out


def run_c(seq, fast):
    r = subprocess.run([CBIN] + (["--fast"] if fast else []) + [str(t) for t in seq],
                       capture_output=True, text=True)
    if r.returncode not in (0, 1) or r.stdout.strip() not in ("True", "False") \
            or (r.returncode == 0) != (r.stdout.strip() == "True"):
        raise RuntimeError(f"C verifier error {r.returncode}: {r.stderr.strip()} on {seq}")
    return r.stdout.strip() == "True"


def main():
    nchains = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    cases = [("vshort2 invalid", [int(t) for t in s.split()], False) for s in vshort2._INVALID]
    rows = [[int(v) for v in r if v.strip()] for r in csv.reader(open(os.path.join(ROOT, "certs/short2/certs.csv")))]
    cases += [("certs/short2", r, True) for r in rows]
    valid, forged = gp_generate(nchains, seed)
    cases += [("pari chain", v, None) for v in valid]
    cases += [("CRT forgery", v, False) for v in forged]
    for label, v, _ in [c for c in cases if c[0] in ("certs/short2", "pari chain", "CRT forgery")]:
        cases += [(f"{label} / {t}", w, None) for t, w in tampers(v)]

    bad, tally, levels = 0, {}, {}
    for label, v, exp in cases:
        py = vshort2.verify(v)
        c0, c1 = run_c(v, False), run_c(v, True)
        key = label.split(" / ")[-1]
        key = key.split(" ", 1)[1] if key.startswith("L") and key[1].isdigit() else key
        t = tally.setdefault(key, [0, 0])
        t[0 if py else 1] += 1
        if label == "pari chain" and py:
            levels[(len(v) - 1) // 3] = levels.get((len(v) - 1) // 3, 0) + 1
        if not (c0 == c1 == py) or (exp is not None and py != exp):
            bad += 1
            print(f"MISMATCH [{label}] faithful={c0} fast={c1} py={py} expected={exp}: {v}")
    for k in sorted(tally):
        print(f"  {k:20s} True {tally[k][0]:4d}  False {tally[k][1]:4d}")
    print(f"  PARI chains accepted, by number of levels: {dict(sorted(levels.items()))}")
    print(f"{len(cases)} cases ({len(valid)} PARI chains, {len(forged)} CRT forgeries), {bad} mismatches")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
