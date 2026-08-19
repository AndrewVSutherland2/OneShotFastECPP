#!/usr/bin/env python3
"""Dickman rho values used in the ecpp-varieties write-up (tables in sections
4-6): solve u rho'(u) = -rho(u-1) by RK in LOG space with linear interpolation
of ln(rho) on a fine grid (interpolating rho itself in absolute terms loses all
relative accuracy once rho < ~1e-9).  Accurate to <0.1% for u <= 6.5 and ~3% at
u = 8 with h = 1/2048; beyond u ~ 9 use published table values (the write-up
used rho(9) = 1.0162e-9, rho(10) = 2.7702e-11, rho(11) = 6.4426e-13,
rho(12) = 1.3116e-14, log-interpolated).

The write-up's u-values: smooth part of #E must exceed L ~ sqrt(p), so
u/2 = (ln p / 2) / ln y with y = n^2 or n^4, n = ceil(log2 p)."""
import math

H = 1.0 / 2048
N = int(15 / H)
lr = [0.0] * (N + 1)                      # ln rho on the grid
for i in range(N + 1):
    u = i * H
    if u <= 1:
        lr[i] = 0.0
    elif u <= 2 + 1e-12:
        lr[i] = math.log(1 - math.log(u))

def rho_interp(x):
    j = x / H; j0 = int(j); fr = j - j0
    return math.exp(lr[j0] * (1 - fr) + lr[min(j0 + 1, N)] * fr)

i2 = int(2 / H)
cur = math.exp(lr[i2])
for i in range(i2 + 1, N + 1):
    u0 = (i - 1) * H
    k1 = -rho_interp(u0 - 1.0) / u0
    k2 = -rho_interp(u0 + H / 2 - 1.0) / (u0 + H / 2)
    k4 = -rho_interp(u0 + H - 1.0) / (u0 + H)
    cur = cur + H * (k1 + 4 * k2 + k4) / 6.0
    lr[i] = math.log(cur)

def rho(u):
    return rho_interp(u)

if __name__ == "__main__":
    print("anchors (literature):")
    for u, v in [(4, 4.9109e-3), (6, 1.9650e-5), (8, 3.2320e-8)]:
        print(f"  rho({u}) = {rho(u):.5g}  (lit {v:.5g}, ratio {rho(u)/v:.4f})")
    print("write-up values (u/2 for smooth-part > sqrt(p)):")
    for tag, u in [("10^20 n2", 2.738), ("10^30 n2", 3.750), ("10^40 n2", 4.708),
                   ("10^50 n2", 5.624), ("10^60 n2", 6.519), ("256b n2", 8.000),
                   ("256b n4", 4.000), ("384b n4", 5.591), ("416b n4", 5.966),
                   ("448b n4", 6.360)]:
        print(f"  {tag:10s} u = {u:6.3f}  rho = {rho(u):.3g}  1/rho = {1/rho(u):.4g}")
