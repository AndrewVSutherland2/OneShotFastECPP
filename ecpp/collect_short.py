#!/usr/bin/env python3
"""collect_short.py -- gather short.gp output files into certs/short/nextprime10c.txt,
verifying each certificate with vsmallECPP and cross-checking every level in PARI.

usage: collect_short.py <dir-with-out<c>.txt files> [out=certs/short/nextprime10c.txt]
"""
import os
import re
import subprocess
import sys
from math import isqrt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vsmallECPP import verify, sieve_primes


def pari_check(seq, n2):
    """Independent PARI check of every level: curve nonsingular, point on curve,
    order exactly o_i, o_i in the window, descent and smoothness conditions."""
    p = seq[0]
    lines = []
    for i in range(1, len(seq), 3):
        A, x, o = seq[i], seq[i + 1], seq[i + 2]
        lines.append(
            "E=ellinit([0,%d,0,1,0],%d); if(#E==0,print(\"SINGULAR\"));"
            "y2=Mod(%d^3+%d*%d^2+%d,%d); B=y2; P=[Mod(%d,%d)/B, Mod(1,%d)/B];"
            "E2=ellinit([0,%d/B,0,1/B^2,0],%d);"
            "print(ellisoncurve(E2,P), \" \", ellorder(E2,P)==%d, \" \", %d%%ellorder(E2,P)==0);"
            % (A, p, x, A, x, x, p, x, p, p, A, p, o, o))
        # next modulus = n^2-rough part of o
        rest = o
        for q in sieve_primes(n2):
            while rest % q == 0:
                rest //= q
        p = rest
    r = subprocess.run(["gp", "-q"], input="default(parisizemax,\"2G\");\n" + "\n".join(lines) + "\n",
                       capture_output=True, text=True)
    out = [l for l in r.stdout.splitlines() if l and "Warning" not in l]
    return all(l.strip() == "1 1 1" for l in out) and len(out) == (len(seq) - 1) // 3


def main():
    d = sys.argv[1]
    out = "certs/short/nextprime10c.txt"
    for a in sys.argv[2:]:
        if a.startswith("out="):
            out = a[4:]
    rows = []
    for c in range(10, 101, 10):
        path = os.path.join(d, "out%d.txt" % c)
        if not os.path.exists(path):
            print("c=%3d: MISSING" % c)
            continue
        txt = open(path).read()
        m = re.search(r"RESULT c=(\d+) curves=(\d+) time=([\d.]+) levels=(\d+)", txt)
        s = re.search(r"^\[([0-9, ]+)\]", txt, re.M)
        if not m or not s:
            print("c=%3d: no result yet" % c)
            continue
        seq = [int(t) for t in s.group(1).split(",")]
        n = seq[0].bit_length()
        ok = verify(seq)
        okp = pari_check(seq, n * n)
        bits = sum(v.bit_length() for v in seq)
        print("c=%3d: %4d bits, %d level(s), %s curves, %6.1fs   vsmallECPP=%s PARI=%s  cert %d bits"
              % (c, n, int(m.group(4)), m.group(2), float(m.group(3)), ok, okp, bits))
        if not (ok and okp):
            sys.exit("VERIFICATION FAILED at c=%d" % c)
        rows.append((c, m.group(2), m.group(3), seq))
    outdir = os.path.dirname(out)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    with open(out, "w") as f:
        f.write("# short ECPP certificates for nextprime(10^c), c = 10, 20, ..., 100\n")
        f.write("# produced by ecpp/short.gp (toy SEA prover); verified by ecpp/vsmallECPP.py\n")
        f.write("# format: p A_0 x_0 o_0 A_1 x_1 o_1 ... A_k x_k o_k\n")
        for (c, curves, t, seq) in rows:
            f.write("# c=%d  (%d bits, %d levels, %s curves, %ss)\n"
                    % (c, seq[0].bit_length(), (len(seq) - 1) // 3, curves, t))
            f.write(" ".join(map(str, seq)) + "\n")
    print("\nwrote %d certificates to %s" % (len(rows), out))


if __name__ == "__main__":
    main()
