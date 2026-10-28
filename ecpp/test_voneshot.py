#!/usr/bin/env python3
"""test_voneshot.py -- differential test of ecpp/voneshot (C) against voneshot.py.

Cases: voneshot.py's own --test vectors; every one-shot certificate in certs/ (n^4,
restricted n^2-smooth, SEA); random small valid certificates built with PARI (prime p,
point of exact order m in the size window, repeated factors allowed); CRT forgeries
over composite p = l1*l2 (a point whose order over Z/pZ is m, with the order split
between the two primes, the shape of the counterexample in ecpp-varieties Rem. 2.3);
and tampered variants of all of these.  Every case must get the same True/False (and
exit status) from both verifiers; the composite-p cases must all be False.

    make -C ecpp voneshot && python3 ecpp/test_voneshot.py [ngen] [seed]
"""
import glob, os, subprocess, sys, tempfile
from math import gcd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, ROOT)
import voneshot

CBIN = os.path.join(ROOT, "ecpp", "voneshot")

GP = r"""
default(parisizemax, "1G");
gencert(bits) = {
  for (tries = 1, 400,
    my(p = randomprime([2^(bits-1), 2^bits - 1]), n = #binary(p), n2 = n^2, n4 = n2^2);
    my(A = random(p));
    if (Mod(A, p)^2 == 4, next);
    my(E = ellinit([0, A, 0, 1, 0], p), P = random(E));
    if (P == [0], next);
    my(o = ellorder(E, P), sp = sqrtint(p), L = sp + 1 + sqrtint(4*sp), cand = List());
    foreach(divisors(o), m, if (m > L, my(f = factor(m)[,1]); if (vecmax(f) < n4 && m < L*f[1], listput(cand, m))));
    if (#cand == 0, next);
    my(m = cand[random(#cand) + 1], Q = ellmul(E, P, o / m), f = factor(m)[,1]);
    return(concat([p, A, lift(Q[1]), m], [q | q <- f~, q > n2])));
  0
}
gencomp(bits) = {
  for (tries = 1, 2000,
    my(l1 = randomprime([2^(bits\2 - 1), 2^(bits\2)]), l2 = randomprime([2^(bits - bits\2 - 1), 2^(bits - bits\2)]));
    if (l1 == l2, next);
    my(p = l1*l2, n = #binary(p), n2 = n^2, n4 = n2^2, A = random(p), x0 = random(p), y1, y2);
    if (gcd(A^2 - 4, p) != 1, next);
    my(rhs = x0^3 + A*x0^2 + x0);
    if (!issquare(Mod(rhs, l1), &y1) || !issquare(Mod(rhs, l2), &y2), next);
    my(E1 = ellinit([0, A, 0, 1, 0], l1), E2 = ellinit([0, A, 0, 1, 0], l2));
    my(P1 = [Mod(x0, l1), y1], P2 = [Mod(x0, l2), y2]);
    my(o = lcm(ellorder(E1, P1), ellorder(E2, P2)), sp = sqrtint(p), L = sp + 1 + sqrtint(4*sp), cand = List());
    foreach(divisors(o), m, if (m > L && m <= p + 1 + sqrtint(4*p), my(f = factor(m)[,1]); if (vecmax(f) < n4 && m < L*f[1], listput(cand, m))));
    if (#cand == 0, next);
    my(m = cand[random(#cand) + 1], Q1 = ellmul(E1, P1, o / m), Q2 = ellmul(E2, P2, o / m));
    if (Q1 == [0] || Q2 == [0], next);
    my(f = factor(m)[,1]);
    return(concat([p, A, lift(chinese(Q1[1], Q2[1])), m], [q | q <- f~, q > n2])));
  0
}
"""


def gp_generate(ngen, seed):
    lines = [f"setrand({seed});"]
    for i in range(ngen):
        lines.append(f'v = gencert({20 + i % 45}); if (v, print("V ", v));')
    for i in range(ngen // 2):
        lines.append(f'v = gencomp({24 + i % 41}); if (v, print("C ", v));')
    with tempfile.NamedTemporaryFile("w", suffix=".gp", delete=False) as f:
        f.write(GP + "\n".join(lines) + "\nquit\n")
        path = f.name
    try:
        out = subprocess.run(["gp", "-q", "-D", "colors=no"], stdin=open(path),
                             capture_output=True, text=True, check=True).stdout
    finally:
        os.unlink(path)
    valid, comp = [], []
    for ln in out.splitlines():
        tag, _, vec = ln.partition(" ")
        if tag in ("V", "C"):               # gp also echoes the function definitions
            vals = [int(t) for t in vec.strip(" []").replace(",", " ").split()]
            (valid if tag == "V" else comp).append(vals)
    return valid, comp


def cert_files():
    out = []
    for f in sorted(glob.glob(os.path.join(ROOT, "certs", "*.txt")) +
                    glob.glob(os.path.join(ROOT, "certs", "restricted", "*.txt")) +
                    glob.glob(os.path.join(ROOT, "certs", "sea", "*.txt"))):
        for ln in open(f):
            if ln.strip() and not ln.startswith("#"):
                out.append((os.path.relpath(f, ROOT), [int(t) for t in ln.split()]))
    return out


def multiple_x(k, x0, A, p):
    X, Z = voneshot.ladder(k, x0, 1, A, p)
    return X * pow(Z, -1, p) % p if gcd(Z, p) == 1 else None


def tampers(v):
    p, A, x0, m, qs = v[0], v[1], v[2], v[3], v[4:]
    n2 = p.bit_length() ** 2
    out = [("x0+1", [p, A, (x0 + 1) % p, m] + qs), ("A+1", [p, (A + 1) % p, x0, m] + qs),
           ("m+1", [p, A, x0, m + 1] + qs), ("m-1", [p, A, x0, m - 1] + qs),
           ("2m", [p, A, x0, 2 * m] + qs), ("p+2", [p + 2, A, x0, m] + qs),
           ("x0=p", [p, A, p, m] + qs), ("A=p", [p, p, x0, m] + qs)]
    if qs:
        out += [("drop q", [p, A, x0, m] + qs[:-1]), ("dup q", [p, A, x0, m] + qs + qs[-1:]),
                ("rev qs", [p, A, x0, m] + qs[::-1]), ("q+2", [p, A, x0, m] + qs[:-1] + [qs[-1] + 2]),
                ("fold 2q", [p, A, x0, m] + qs[:-1] + [2 * qs[-1]]),
                ("extra q", [p, A, x0, m] + qs + [qs[-1] + 2])]
    else:
        out.append(("bogus q", [p, A, x0, m, n2 + 1]))
    # x0 -> x(rP) for each prime r | m: m*P' = O still, but (m/r)P' = O, so the tree
    # must catch it.  Primes <= n^2 are found by trial division of m, as the verifier does.
    rs = [r for r in range(2, min(n2, 5000) + 1) if m % r == 0 and all(r % s for s in range(2, int(r**0.5) + 1))]
    for r in rs[:3] + qs[:2]:
        x = multiple_x(r, x0, A, p)
        if x is not None:
            out.append((f"x({r}P)", [p, A, x, m] + qs))
    return out


def run_c(v):
    r = subprocess.run([CBIN] + [str(t) for t in v], capture_output=True, text=True)
    if r.returncode not in (0, 1) or r.stdout.strip() not in ("True", "False"):
        raise RuntimeError(f"C verifier error {r.returncode}: {r.stderr.strip()} on {v}")
    if (r.returncode == 0) != (r.stdout.strip() == "True"):
        raise RuntimeError(f"exit status disagrees with output on {v}")
    return r.stdout.strip() == "True"


def main():
    ngen = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    cases = []          # (label, vals, expected or None)
    for v, e in [([101, 3, 24, 24], True), ([1003, 3, 24, 24], False), ([2, 0, 0, 1], False),
                 ([3, 0, 0, 4], False), ([101, 104, 24, 24], False), ([101, 3, 125, 24], False),
                 ([101, 3, 24, 123], False), ([101, 3, 24, 24, 5], False),
                 ([10000019, 6873344, 1, 3308, 1654], False), ([10000019, 6873344, 1, 3308, 827], False),
                 ([4410667997551, 1365834658413, 107710304518, 4200232, 525029], False)]:
        cases.append(("voneshot --test", v, e))
    for ln in voneshot._TEST_VECTORS.strip().splitlines():
        cases.append(("voneshot vector", [int(t) for t in ln.split()], True))
    for f, v in cert_files():
        cases.append((f, v, True))
    valid, comp = gp_generate(ngen, seed)
    cases += [("pari valid", v, True) for v in valid]
    cases += [("CRT forgery", v, False) for v in comp]
    base = [c for c in cases if c[2] is True]
    for label, v, _ in base + [c for c in cases if c[0] == "CRT forgery"]:
        cases += [(f"{label} / {t}", w, None) for t, w in tampers(v)]

    bad = 0
    tally = {}
    for label, v, exp in cases:
        py = voneshot.verify(v[0], v[1], v[2], v[3], tuple(v[4:]))
        c = run_c(v)
        key = label.split(" / ")[-1] if " / " in label else label.split("/")[0]
        key = "x(rP), r | m" if key.startswith("x(") else key
        t = tally.setdefault(key, [0, 0])
        t[0 if c else 1] += 1
        if c != py or (exp is not None and py != exp):
            bad += 1
            print(f"MISMATCH [{label}] C={c} py={py} expected={exp}: {v}")
    for k in sorted(tally):
        print(f"  {k:24s} True {tally[k][0]:4d}  False {tally[k][1]:4d}")
    print(f"{len(cases)} cases ({len(valid)} PARI-valid, {len(comp)} CRT forgeries), {bad} mismatches")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
