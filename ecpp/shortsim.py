#!/usr/bin/env python3
"""shortsim.py -- a heuristic model of one level of the short-ECPP search.

What the model is
-----------------
A candidate order N of a level with modulus p is stripped of its y0-smooth
part, leaving a tail T.  The level is won when T = (primes that ECM removes) * q
with q a prime the certificate can descend to, i.e. with L/m < q < r(m) L/m for
some divisor m of the n^2-smooth part s of N (L ~ sqrt p, r(m) the least prime
of m).  Write lam = log2 L and, for a prime q, delta = lam - log2 q; the set of
delta that some divisor covers is the candidate's *coverage* Cov(s).

Prime factors of T between y0 and 2^D (D = 160 bits, beyond anything ECM
finds) are modelled as a Poisson process with intensity db/b in the bit size b,
independent of the rest; the cofactor R of T by those primes is a 2^D-rough
number, prime with the density Buchstab's function gives.  A prime of b bits
survives a tier of c curves with probability (1-pi(b))^c, pi measured for
gmp-ecm (pi_calib.json) or taken from the two-stage Dickman function.  All sums
of prime sizes are then compound-Poisson laws, computed exactly on a 1-bit
grid, which gives for a tail of t bits known to be composite after tier j:

  Z_j(t)       weight of that observation,
  transitions  to a smaller tail (composite or prime) at tier j+1,
  the chance   that the tail becomes a prime inside the coverage (a win).

Two things are built on it:
  * emulate(): Monte Carlo of the reference scheduler's tiers and admission
    caps on a pool, for comparison with campaign logs (split rates, deaths,
    winners per tier, core-seconds);
  * optimal(): the policy that minimises expected cost per winner when fresh
    candidates can always be bought at a fixed price -- a stopping problem per
    candidate, solved by dynamic programming over (tier, tail size) with the
    value of a win found by bisection.

Nothing here touches the prover; it reads candidate data and prints tables.
"""

import json
import math
import os
import sys

import numpy as np

LN2 = math.log(2.0)
D = 160                       # largest "medium" prime, bits
NB = 1536                     # grid: sizes 0..NB-1 bits
NFFT = 4096


def set_grid(nbits):
    """Size the bit grid for moduli of nbits bits (the default holds 1535 bits).
    Everything built afterwards (Ladder, IndexTable, ...) uses the new grid; a
    transform of twice the grid length keeps the convolutions free of wrap-around."""
    global NB, NFFT
    NFFT = 1 << max(11, (2 * (nbits + 8) - 1).bit_length())
    NB = NFFT // 2


# --------------------------------------------------------------- special functions
def _dickman_table(umax=40.0, h=1.0 / 2048):
    n = int(umax / h) + 1
    u = np.arange(n) * h
    rho = np.ones(n)
    k1 = int(round(1 / h))
    for i in range(k1 + 1, n):
        rho[i] = rho[i - 1] - 0.5 * h * (rho[i - 1 - k1] / u[i - 1] + rho[i - k1] / u[i])
    return u, rho, h


_U, _RHO, _H = _dickman_table()
_RHO_CUM = None


def rho(u):
    """Dickman's function."""
    u = np.asarray(u, dtype=float)
    return np.where(u <= 1, 1.0, np.maximum(np.interp(u, _U, _RHO, right=0.0), 0.0))


def mu(alpha, beta):
    """Probability that a random integer x has all prime factors below x^(1/alpha)
    except possibly one below x^(beta/alpha) (stage 1 bound, stage 2 bound)."""
    if alpha <= 1:
        return 1.0
    lo, hi = max(alpha - beta, 0.0), alpha - 1.0
    val = float(rho(alpha))
    if hi > lo:
        t = np.linspace(lo, hi, max(64, int((hi - lo) / _H) + 1))
        f = rho(t) / (alpha - t)
        val += float(np.sum((f[1:] + f[:-1]) * 0.5 * (t[1] - t[0])))
    return min(1.0, val)


def _buchstab_table(umax=40.0, h=1.0 / 1024):
    n = int(umax / h) + 1
    u = np.arange(n) * h
    w = np.zeros(n)
    k1 = int(round(1 / h))
    w[k1:2 * k1 + 1] = 1.0 / u[k1:2 * k1 + 1]
    uw = u * w
    for i in range(2 * k1 + 1, n):            # (u w(u))' = w(u-1)
        uw[i] = uw[i - 1] + 0.5 * h * (w[i - 1 - k1] + w[i - k1])
        w[i] = uw[i] / u[i]
    return u, w


_BU, _BW = _buchstab_table()


def buchstab(u):
    """Buchstab's function: the density of z-rough integers near z^u is w(u)/ln z."""
    u = np.asarray(u, dtype=float)
    return np.where(u < 1, 0.0, np.interp(u, _BU, _BW))


# --------------------------------------------------------------- ECM success model
# gmp-ecm 7 default stage-2 bounds for the B1 values the prover uses
B2_DEFAULT = {2000: 147396, 11000: 1873422, 50000: 12746592, 250000: 128992510,
              1000000: 1045563762, 3000000: 5706890290, 10000000: 35133391030,
              30000000: 144289975846}
# gmp-ecm P-1 default stage-2 bounds
B2_PM1 = {500000: 1.3e9, 10000000: 1.0e11}


def ecm_pi_theory(bits, b1, b2=None, lndelta=math.log(24.0)):
    """Two-stage smoothness probability of a curve order near 2^bits/delta."""
    b2 = b2 or B2_DEFAULT.get(b1, 100.0 * b1)
    out = np.empty(len(bits))
    for i, b in enumerate(bits):
        a = max(b * LN2 - lndelta, 0.0) / math.log(b1)
        out[i] = mu(a, math.log(b2) / math.log(b1))
    return out


class EcmModel:
    """pi(b) per tier: measured cells where the calibration has enough hits,
    the theoretical curve (rescaled to join the last good cell) beyond them."""

    def __init__(self, calib=None, lndelta=math.log(24.0)):
        self.cells = {}
        self.lndelta = lndelta
        if calib and os.path.exists(calib):
            for c in json.load(open(calib)):
                self.cells.setdefault((c["method"], c["b1"]), []).append(c)

    def pi(self, method, b1, bits):
        bits = np.asarray(bits, dtype=float)
        if method == "pm1":
            th = ecm_pi_theory(bits, b1, B2_PM1.get(b1, 100.0 * b1), math.log(2.0))
        else:
            th = ecm_pi_theory(bits, b1, None, self.lndelta)
        cells = sorted(self.cells.get((method, b1), []), key=lambda c: c["bits"])
        good = [c for c in cells if c["hits"] >= 8 and c["curves"] > 0]
        if len(good) < 3:
            return th
        xb = np.array([c["bits"] for c in good], dtype=float)
        pm = np.array([min(c["hits"] / c["curves"], 0.999) for c in good])
        th_at = ecm_pi_theory(xb, b1, B2_PM1.get(b1) if method == "pm1" else None,
                              math.log(2.0) if method == "pm1" else self.lndelta)
        # work with the log of the ratio to theory: smooth, and extrapolates sensibly
        lr = np.log(pm / np.maximum(th_at, 1e-300))
        ratio = np.exp(np.interp(bits, xb, lr))        # constant beyond the measured range
        return np.clip(th * ratio, 0.0, 1.0)


# --------------------------------------------------------------- compound Poisson tools
def expconv(a):
    """sum_k a^{*k}/k! on the size grid (index 0 carries the k = 0 term)."""
    f = np.fft.rfft(a, NFFT)
    g = np.fft.irfft(np.exp(f), NFFT)[:NB]
    return np.maximum(g, 0.0)


def conv(a, b):
    return np.fft.irfft(np.fft.rfft(a, NFFT) * np.fft.rfft(b, NFFT), NFFT)[:NB]


class Ladder:
    """Tiers = list of (method, B1, curves); precomputes everything that does not
    depend on the candidate."""

    def __init__(self, tiers, ecm, y0bits=32, cost=None, pm1_ratio=8.0):
        self.tiers = list(tiers)
        self.y0 = y0bits
        bits = np.arange(NB, dtype=float)
        w = np.zeros(NB)                       # Poisson mass of each 1-bit cell: int db/b
        lo = np.maximum(bits - 0.5, y0bits)
        hi = bits + 0.5
        cell = (bits > y0bits) & (bits <= D)
        w[cell] = np.log(hi[cell] / lo[cell])
        self.w = w
        mb = np.arange(NB)
        surv = np.ones(NB)
        self.S = [surv.copy()]                 # S[j] = survival through the first j tiers
        pm1_seen = np.zeros(NB)                # best P-1 success so far (deterministic method)
        self.pi = []
        for (method, b1, curves) in self.tiers:
            pi = np.zeros(NB)
            pi[cell] = ecm.pi(method, b1, mb[cell])
            self.pi.append(pi)
            if method == "pm1":
                new = np.maximum(pi, pm1_seen)
                surv = surv * (1 - new) / np.maximum(1 - pm1_seen, 1e-12)
                pm1_seen = new
            else:
                surv = surv * (1 - pi) ** curves
            self.S.append(surv.copy())
        J = len(self.tiers)
        self.nu = [self.S[j] * w for j in range(J + 1)]            # remaining-prime intensity
        self.Lam = [float(v.sum()) for v in self.nu]
        self.G = [expconv(v) for v in self.nu]
        self.phi = [self.nu[j] - self.nu[j + 1] for j in range(J)]  # found at tier j+1
        self.Gphi = [expconv(v) for v in self.phi]
        self.Lphi = [float(v.sum()) for v in self.phi]
        t = np.arange(NB, dtype=float)
        self.omega = buchstab(t / D)                                # rough weight of a cofactor of t bits
        self.prime = np.where(t >= D, D / np.maximum(t, 1.0), 0.0)  # prime weight (= omega for t < 2D)
        self.Z = []
        for j in range(J + 1):
            z = math.exp(-self.Lam[j]) * (conv(self.G[j], self.omega) - self.prime)
            self.Z.append(np.maximum(z, 1e-300))
        self.cost_fn = cost
        self.pm1_ratio = pm1_ratio

    def tier_cost(self, j, tbits):
        method, b1, curves = self.tiers[j]
        c = self.cost_fn(b1, tbits)
        return c * curves / (self.pm1_ratio if method == "pm1" else 1.0)

    def split_rate(self, j):
        """P(at least one factor found at tier j+1 | reached it), ignoring the
        (weak) dependence on the tail size."""
        return 1.0 - math.exp(-self.Lphi[j])

    def p_prime_fresh(self, tbits):
        """P(a fresh stripped tail of tbits bits is prime)."""
        t = int(round(tbits))
        tot = math.exp(-self.Lam[0]) * conv(self.G[0], self.omega)[t]
        return math.exp(-self.Lam[0]) * self.prime[t] / max(tot, 1e-300)


# --------------------------------------------------------------- cost model
def make_cost(sec1030, exponent=1.5):
    """sec1030: {B1: core-seconds per gmp-ecm curve on a 1030-bit modulus}."""
    def cost(b1, tbits):
        return sec1030[b1] * (np.asarray(tbits, dtype=float) / 1030.0) ** exponent
    return cost


# Measured on a c8a.48xlarge (Zen 5, 192 cores, every core busy) with the
# gmp-ecm 7.0.5 that Ubuntu ships, default B2: core-seconds per curve on a
# 1030-bit modulus.  B1 <= 11000 and 3e6 are the idle-core times scaled by the
# loaded/idle ratio 1.33 seen at the three B1 that were run loaded; 1e7 and 3e7
# are extrapolated (x3 per step) and have not been measured.
SEC1030_LOADED = {2000: 0.0134, 11000: 0.064, 50000: 0.263, 250000: 1.215,
                  1000000: 4.95, 3000000: 13.5, 10000000: 40.0, 30000000: 120.0}
# Effective cost in the August 2026 campaign's 1030-bit run (wall x 192 cores /
# curves): includes launch overhead and idle cores.
SEC1030_CAMPAIGN = {2000: 0.048, 11000: 0.101, 50000: 0.336, 250000: 1.36,
                    1000000: 4.68, 3000000: 12.0}


# --------------------------------------------------------------- candidates
def coverage_cells(sfac, lam_floor_bits=0.0, bcap=None, radlim=None, maxdiv=200000, width=192):
    """Coverage of a candidate with n^2-smooth part s = prod p^e (dict sfac):
    cov[k] = measure of {delta in [k-1/2, k+1/2) : some divisor m > 1 of s has
    log2 m - log2 r(m) < delta < log2 m}, r(m) the least prime of m.  With the
    revised format's limits, m uses only primes <= bcap and log2 rad(m) <= radlim."""
    primes = sorted(q for q in sfac if (bcap is None or q <= bcap))
    res = 8                                    # sub-cells per bit
    diff = np.zeros(width * res + 2)
    nd = 0
    # divisors grouped by least prime: m = lp^a * (divisor of the part above lp)
    for i, lp in enumerate(primes):
        llp = math.log2(lp)
        ups = [(0.0, 0.0)]                     # (log2 of divisor of the larger primes, log2 of its radical)
        for q in primes[i + 1:]:
            lq = math.log2(q)
            new = []
            for (lm, lr) in ups:
                for e in range(1, sfac[q] + 1):
                    if radlim is not None and lr + lq > radlim:
                        break
                    new.append((lm + e * lq, lr + lq))
            ups += new
            if len(ups) > maxdiv:
                break
        nd += len(ups) * sfac[lp]
        if nd > maxdiv:                        # very smooth s: call everything up to log2 s covered
            tot = sum(sfac[q] * math.log2(q) for q in primes)
            diff[:] = 0
            diff[0] += 1
            diff[min(int(tot * res), width * res)] -= 1
            break
        arr = np.array(ups)
        for a in range(1, sfac[lp] + 1):
            if radlim is not None and np.any(arr[:, 1] + llp > radlim):
                top = arr[arr[:, 1] + llp <= radlim, 0] + a * llp
            else:
                top = arr[:, 0] + a * llp
            lo = np.clip(((top - llp) * res).astype(int), 0, width * res)
            hi = np.clip((top * res).astype(int), 0, width * res)
            np.add.at(diff, lo, 1)
            np.add.at(diff, hi, -1)
    covered = (np.cumsum(diff)[:width * res] > 0).astype(float)
    # integer cells [k-1/2, k+1/2)
    pad = np.concatenate([np.zeros(res // 2), covered, np.zeros(res)])
    cells = pad[:width * res].reshape(width, res).mean(axis=1)
    return cells


class Pool:
    """Candidate ensemble: integer tail sizes, coverage classes, per-candidate data."""

    def __init__(self, tbits, cov, lam, prp=None, lns=None):
        self.t = np.asarray(tbits, dtype=float)
        self.lns = None if lns is None else np.asarray(lns, dtype=float)   # ln s
        self.cov = np.asarray(cov)             # (ncand, width) coverage cells
        self.lam = lam                         # log2 L
        self.prp = np.zeros(len(self.t), dtype=bool) if prp is None else np.asarray(prp, dtype=bool)
        self.w = self.cov.sum(axis=1)          # covered measure, bits
        self.cls = np.minimum(np.round(self.w).astype(int), self.cov.shape[1] - 1)

    def class_profiles(self):
        out = {}
        for c in np.unique(self.cls):
            sel = self.cls == c
            out[int(c)] = (self.cov[sel].mean(axis=0), int(sel.sum()))
        return out

    def window(self, covcells):
        """Coverage cells -> weight per tail size r (bits): r = lam - delta."""
        win = np.zeros(NB)
        lam = int(round(self.lam))
        k = np.arange(len(covcells))
        r = lam - k
        ok = (r >= 0) & (r < NB)
        win[r[ok]] = covcells[k[ok]]
        return win


# --------------------------------------------------------------- optimal stopping
def solve_class(lad, win, price, reward):
    """Value of a composite tail by tier and size for one coverage class.

    win[r]   coverage weight of a prime tail of r bits
    price    not used here (charged once per candidate by the caller)
    reward   value of a win, in core-seconds
    Returns (V, go): V[j][t] = value of holding a composite tail of t bits after
    j tiers under the best continuation, go[j][t] = whether to run tier j+1."""
    J = len(lad.tiers)
    t = np.arange(NB, dtype=float)
    V = [None] * (J + 1)
    go = [None] * J
    V[J] = np.zeros(NB)
    for j in range(J - 1, -1, -1):
        e1 = math.exp(-lad.Lam[j + 1])
        primeterm = e1 * lad.prime * reward * win            # tail became a prime of that size
        U = lad.Z[j + 1] * V[j + 1] + primeterm
        cont = math.exp(-lad.Lphi[j]) * (conv(lad.Gphi[j], U) - lad.Gphi[j][0] * primeterm) / lad.Z[j]
        cont = cont - lad.tier_cost(j, np.maximum(t, 64.0))
        cont[:D] = -1.0                                       # tails below D bits cannot hold a q
        go[j] = cont > 0
        V[j] = np.where(go[j], cont, 0.0)
    return V, go


def fresh_value(lad, pool, cls_profile, tvals, V0, reward):
    """Expected value of fresh candidates of one class with the given tail sizes
    (prime tails are settled at once, for free)."""
    win = pool.window(cls_profile)
    ti = np.clip(np.round(tvals).astype(int), 0, NB - 1)
    tot = math.exp(-lad.Lam[0]) * conv(lad.G[0], lad.omega)
    pprime = math.exp(-lad.Lam[0]) * lad.prime[ti] / np.maximum(tot[ti], 1e-300)
    return pprime * reward * win[ti] + (1 - pprime) * V0[ti]


def optimal(lad, pool, price, lo=1.0, hi=1e9, iters=40):
    """Dinkelbach/bisection: the cost per win c* is the reward at which a fresh
    candidate is worth exactly its price.  Returns (c*, policy) where policy[cls]
    = go arrays, plus per-class fresh values at c*."""
    profiles = pool.class_profiles()

    def net(reward):
        tot = 0.0
        for c, (prof, cnt) in profiles.items():
            win = pool.window(prof)
            V, _ = solve_class(lad, win, price, reward)
            sel = pool.cls == c
            tot += float(fresh_value(lad, pool, prof, pool.t[sel], V[0], reward).sum())
        return tot / len(pool.t) - price
    for _ in range(iters):
        mid = math.sqrt(lo * hi)
        if net(mid) > 0:
            hi = mid
        else:
            lo = mid
    reward = math.sqrt(lo * hi)
    policy = {}
    for c, (prof, cnt) in profiles.items():
        V, go = solve_class(lad, pool.window(prof), price, reward)
        policy[c] = go
    return reward, policy


def corr(a, b):
    """c[t] = sum_y a[y] b[t+y]."""
    return np.fft.irfft(np.conj(np.fft.rfft(a, NFFT)) * np.fft.rfft(b, NFFT), NFFT)[:NB]


def evaluate(lad, pool, policy, classes=None):
    """Exact expectations, per fresh candidate, of a policy go[cls][j][t] (run
    tier j+1 on a composite tail of t bits in coverage class cls): per-tier
    probability of running, cost and win probability, plus wins at the strip.
    Returns dict(run=[..], cost=[..], win=[..], win0=..), each a pool average."""
    J = len(lad.tiers)
    run = np.zeros(J)
    cost = np.zeros(J)
    win = np.zeros(J)
    win0 = 0.0
    t = np.arange(NB, dtype=float)
    tot0 = math.exp(-lad.Lam[0]) * conv(lad.G[0], lad.omega)
    profiles = pool.class_profiles()
    for c, (prof, cnt) in profiles.items():
        if classes is not None and c not in classes:
            continue
        w = pool.window(prof)
        sel = pool.cls == c
        ti = np.clip(np.round(pool.t[sel]).astype(int), 0, NB - 1)
        hist = np.bincount(ti, minlength=NB).astype(float)        # fresh tails by size
        pprime = math.exp(-lad.Lam[0]) * lad.prime / np.maximum(tot0, 1e-300)
        win0 += float((hist * pprime * w).sum())
        fj = hist * (1 - pprime)                                  # live composite tails after 0 tiers
        go = policy[c]
        for j in range(J):
            m = fj * go[j]
            run[j] += m.sum()
            cost[j] += float((m * lad.tier_cost(j, np.maximum(t, 64.0))).sum())
            cc = corr(lad.Gphi[j], m / lad.Z[j])
            e = math.exp(-lad.Lphi[j])
            win[j] += e * math.exp(-lad.Lam[j + 1]) * float((lad.prime * w * (cc - lad.Gphi[j][0] * m / lad.Z[j])).sum())
            fj = np.maximum(e * lad.Z[j + 1] * cc, 0.0)
    ncand = float(len(pool.t))
    return dict(run=run / ncand, cost=cost / ncand, win=win / ncand, win0=win0 / ncand)


def reference_policy(lad, pool, fractions, slope=0.03):
    """The reference scheduler's admission as a per-class threshold policy for
    one pass over a fresh pool: tier j+1 takes the fraction fractions[j] of the
    pool with the largest index ln ln s - slope * (tail bits above lam) among
    the tails still alive; 1 = everything alive."""
    J = len(lad.tiers)
    t = np.arange(NB, dtype=float)
    profiles = pool.class_profiles()
    classes = sorted(profiles)
    ncand = float(len(pool.t))
    tot0 = math.exp(-lad.Lam[0]) * conv(lad.G[0], lad.omega)
    pprime = math.exp(-lad.Lam[0]) * lad.prime / np.maximum(tot0, 1e-300)
    index, f = {}, {}
    for c in classes:
        sel = pool.cls == c
        lnln = math.log(max(float(np.mean(pool.lns[sel])), 1.1))
        index[c] = lnln - slope * np.maximum(0.0, t - pool.lam)
        ti = np.clip(np.round(pool.t[sel]).astype(int), 0, NB - 1)
        f[c] = np.bincount(ti, minlength=NB).astype(float) * (1 - pprime)   # live composite tails
    policy = {c: [] for c in classes}
    for j in range(J):
        th = -np.inf
        if fractions[j] < 1:
            idx = np.concatenate([index[c] for c in classes])
            m = np.concatenate([f[c] for c in classes])
            order = np.argsort(-idx, kind="stable")
            cum = np.cumsum(m[order])
            k = int(np.searchsorted(cum, fractions[j] * ncand))
            th = idx[order[min(k, len(order) - 1)]]
        e = math.exp(-lad.Lphi[j])
        for c in classes:
            go = index[c] >= th
            policy[c].append(go)
            f[c] = np.maximum(e * lad.Z[j + 1] * corr(lad.Gphi[j], f[c] * go / lad.Z[j]), 0.0)
    return policy


# --------------------------------------------------------------- Monte Carlo emulation
def sample_tails(lad, pool, rng, idx=None):
    """For each candidate with a composite tail, draw its medium primes (bit
    sizes, as a list per candidate) and whether the remaining cofactor is prime."""
    idx = np.arange(len(pool.t)) if idx is None else idx
    cellb = np.nonzero(lad.w > 0)[0]
    wts = lad.w[cellb]
    Lam = wts.sum()
    pw = wts / Lam
    primes, rprime = [], []
    for i in idx:
        t = pool.t[i]
        while True:
            k = rng.poisson(Lam)
            bs = rng.choice(cellb, size=k, p=pw) + rng.uniform(-0.5, 0.5, size=k) if k else np.zeros(0)
            r = t - bs.sum()
            if r < D:
                continue                              # not a D-rough cofactor (or nothing left)
            om = float(buchstab(r / D))
            if rng.random() > om:                     # omega <= 1: accept with the rough weight
                continue
            isp = rng.random() < (D / r) / om
            if k == 0 and isp:
                continue                              # the tail is known to be composite
            break
        primes.append(np.sort(bs))
        rprime.append(isp)
    return primes, np.array(rprime)


def emulate(lad, pool, rng, admit, hb=None, verbose=True, idx=None):
    """One pass of a tiered ladder over a pool of fresh candidates.

    admit(j, live_idx, tail_bits) -> indices admitted to tier j+1 (the caller
    encodes caps, ranking or a stopping policy).  Returns per-tier statistics:
    tails run, split, wins, deaths, core-seconds."""
    idx = np.arange(len(pool.t)) if idx is None else np.asarray(idx)
    comp = idx[~pool.prp[idx]]
    primes, rprime = sample_tails(lad, pool, rng, comp)
    tail = {int(i): float(pool.t[i]) for i in comp}
    rem = {int(i): list(pr) for i, pr in zip(comp, primes)}
    isp = {int(i): bool(v) for i, v in zip(comp, rprime)}
    live = set(int(i) for i in comp)
    lam = pool.lam
    stats = []
    wins_fresh = 0
    for i in idx[pool.prp[idx]]:                      # prime tails: settled at classification
        k = int(round(lam - pool.t[i]))
        if 0 <= k < pool.cov.shape[1] and rng.random() < pool.cov[i, k]:
            wins_fresh += 1
    for j, (method, b1, curves) in enumerate(lad.tiers):
        liv = np.array(sorted(live), dtype=int)
        if len(liv) == 0:
            stats.append(dict(tier=j, run=0, split=0, wins=0, dead=0, cost=0.0))
            continue
        adm = admit(j, liv, np.array([tail[int(i)] for i in liv]))
        nsplit = wins = dead = 0
        cost = 0.0
        pi = lad.pi[j]
        for i in adm:
            i = int(i)
            cost += float(lad.tier_cost(j, tail[i]))
            if not rem[i]:
                continue
            bs = np.array(rem[i])
            pj = pi[np.clip(np.round(bs).astype(int), 0, NB - 1)]
            if method == "pm1":
                found = rng.random(len(bs)) < pj      # (a repeat at a larger B1 is handled by S only)
            else:
                found = rng.random(len(bs)) < 1 - (1 - pj) ** curves
            if not found.any():
                continue
            nsplit += 1
            tail[i] -= float(bs[found].sum())
            rem[i] = list(bs[~found])
            if not rem[i]:
                if isp[i]:                            # the tail is now a prime: win or dead end
                    live.discard(i)
                    k = int(round(lam - tail[i]))
                    if 0 <= k < pool.cov.shape[1] and rng.random() < pool.cov[i, k]:
                        wins += 1
                    else:
                        dead += 1
            if tail[i] < lam - pool.w[i] - 2:          # s * tail <= L: nothing left to find
                if i in live:
                    live.discard(i)
                    dead += 1
        stats.append(dict(tier=j, run=len(adm), split=nsplit, wins=wins, dead=dead, cost=cost))
    return stats, wins_fresh


# --------------------------------------------------------------- break-even rewards
def class_window(c, lam):
    """Window weights of coverage class c: the n^2-smooth part s of the order has
    log2 s in [c - 1/2, c + 1/2).  The divisors m of s cover every delta in
    (0, log2 s) -- take the primes of s in descending order: the powers of the
    largest tile (0, e1 log p1), and each further prime extends every tile
    boundary to the right by its own power -- so the profile of a class is an
    initial segment, with half weight on the two end cells.  (The revised
    format's radical cap could cut it short; it binds for s above 2^80 only.)"""
    win = np.zeros(NB)
    lami = int(round(lam))
    for k in range(c + 1):
        r = lami - k
        if 0 <= r < NB:
            win[r] = 0.5 if k in (0, c) else 1.0
    return win


class IndexTable:
    """Break-even rewards of every state of a candidate.

    r(c, j, t) is the least value of a win, in core-seconds, at which running
    sweep j+1 on a composite tail of t bits in coverage class c pays under the
    best continuation; its reciprocal is the candidate's index (win probability
    per core-second under the best stopping rule), so sorting by r is the order
    in which work should be done, and "r <= R" is the stopping rule when a win is
    worth R.  The table does not depend on what fresh orders cost.  pwin(c, j, t)
    is the chance that this one sweep wins, which lets a scheduler add up the wins
    the model expects from the work it has actually done.  Alongside,
    for a fresh composite tail and each reward of the grid: value V0 and expected
    ECM cost C0 under the rule "r <= R", from which the price at which buying
    fresh orders breaks even (threshold) and the expected wins and cost per
    order follow for any mix of fresh orders.

    Solved on a geometric grid of rewards by the dynamic programme of
    solve_class, all classes at once; the break-even point between two grid
    values is interpolated linearly in the continuation value."""

    def __init__(self, lad, lam, cmax=100, rlo=30.0, rhi=3.0e7, ratio=1.3, threads=1):
        self.lad, self.lam, self.cmax = lad, lam, cmax
        G = int(math.ceil(math.log(rhi / rlo) / math.log(ratio))) + 1
        self.rgrid = rlo * ratio ** np.arange(G)
        J = len(lad.tiers)
        t = np.arange(NB, dtype=float)
        self._tc = np.array([lad.tier_cost(j, np.maximum(t, 64.0)) for j in range(J)])
        self._Gf = [np.fft.rfft(lad.Gphi[j], NFFT) for j in range(J)]
        self._bad = [(t < D) | (lad.Z[j] < 1e-200) for j in range(J)]
        wins = np.array([class_window(c, lam) for c in range(cmax + 1)])
        # a tail shorter than the window's far end holds no admissible prime: such
        # a candidate is dead, and the model's weights there are numerical noise
        dead = t[None, :] < (int(round(lam)) - np.arange(cmax + 1) - 1)[:, None]
        nchunk = max(1, min(threads, (cmax + 1) // 4))
        parts = np.array_split(np.arange(cmax + 1), nchunk)
        if nchunk > 1:
            import concurrent.futures
            with concurrent.futures.ThreadPoolExecutor(max_workers=nchunk) as ex:
                res = list(ex.map(lambda cs: self._solve(wins[cs], dead[cs]), parts))
        else:
            res = [self._solve(wins, dead)]
        self.rmin = np.concatenate([x[0] for x in res], axis=0)       # (class, sweep, bits)
        self.V0 = np.concatenate([x[1] for x in res], axis=1)         # (reward, class, bits)
        self.C0 = np.concatenate([x[2] for x in res], axis=1)
        # chance that one sweep turns the tail into a prime inside the window
        self.pw = np.zeros((cmax + 1, J, NB), dtype=np.float32)
        for j in range(J):
            pt = math.exp(-lad.Lam[j + 1]) * lad.prime[None, :] * wins
            Y = np.fft.irfft(np.fft.rfft(pt, NFFT, axis=1) * self._Gf[j][None, :], NFFT, axis=1)[:, :NB]
            w = math.exp(-lad.Lphi[j]) * (Y - lad.Gphi[j][0] * pt) / lad.Z[j][None, :]
            w[:, self._bad[j]] = 0.0
            w[dead] = 0.0
            self.pw[:, j, :] = np.clip(w, 0.0, 1.0)

    def _solve(self, wins, dead):
        lad = self.lad
        J = len(lad.tiers)
        n = wins.shape[0]
        G = len(self.rgrid)
        rmin = np.full((n, J, NB), np.inf, dtype=np.float32)
        prev = np.zeros((n, J, NB))
        V0 = np.zeros((G, n, NB), dtype=np.float32)
        C0 = np.zeros((G, n, NB), dtype=np.float32)
        rp = 0.0
        for g, r in enumerate(self.rgrid):
            V = np.zeros((n, NB))
            C = np.zeros((n, NB))
            for j in range(J - 1, -1, -1):
                pt = (math.exp(-lad.Lam[j + 1]) * r) * lad.prime[None, :] * wins
                zn = lad.Z[j + 1][None, :]
                X = np.concatenate([zn * V + pt, zn * C])
                Y = np.fft.irfft(np.fft.rfft(X, NFFT, axis=1) * self._Gf[j][None, :], NFFT, axis=1)[:, :NB]
                e = math.exp(-lad.Lphi[j])
                iz = 1.0 / lad.Z[j][None, :]
                cont = e * (Y[:n] - lad.Gphi[j][0] * pt) * iz - self._tc[j][None, :]
                cont[:, self._bad[j]] = -1.0
                cont[dead] = -1.0
                np.clip(cont, -1e30, 1e30, out=cont)
                go = cont > 0
                new = go & ~np.isfinite(rmin[:, j, :])
                if new.any():
                    if g:
                        pc = prev[:, j, :][new]
                        rmin[:, j, :][new] = rp + (r - rp) * (-pc) / np.maximum(cont[new] - pc, 1e-300)
                    else:
                        rmin[:, j, :][new] = r
                prev[:, j, :] = cont
                V = np.where(go, cont, 0.0)
                C = np.where(go, np.clip(e * Y[n:] * iz, 0.0, 1e30) + self._tc[j][None, :], 0.0)
            V0[g] = V
            C0[g] = C
            rp = r
        return rmin, V0, C0

    def r(self, cls, j, tbits):
        return float(self.rmin[min(cls, self.cmax), j, min(tbits, NB - 1)])

    def pwin(self, cls, j, tbits):
        """Probability that sweep j+1 wins from this state."""
        return float(self.pw[min(cls, self.cmax), j, min(tbits, NB - 1)])

    def fresh(self, H):
        """H[c][t] = number of fresh composite tails of t bits in class c.  Returns
        (F, E): summed value and summed expected ECM cost at each grid reward."""
        return (np.tensordot(self.V0, H, axes=([1, 2], [0, 1])),
                np.tensordot(self.C0, H, axes=([1, 2], [0, 1])))

    def threshold(self, price, F):
        """The reward at which a fresh order (mean value F by grid reward) is worth
        its price: the cost per win that the supply can sustain."""
        ok = F > 0
        if not ok.any():
            return float(self.rgrid[-1])
        return float(np.exp(np.interp(math.log(max(price, 1e-12)), np.log(F[ok]), np.log(self.rgrid[ok]))))

    def at(self, R, A):
        """Interpolate a per-reward array (F or E) at the reward R."""
        return float(np.interp(math.log(R), np.log(self.rgrid), A))


# --------------------------------------------------------------- data loading
def load_candidates(path, p, n2, v2=False, limit=None):
    """Lines 'd S T s prp' (shortECPP's strip output) plus a file path + '.fac'
    holding the factorisation of each s as 'q^e q^e ...' -> Pool."""
    n = p.bit_length()
    rt = math.isqrt(p)
    L = rt + 1 + math.isqrt(4 * rt)
    lam = math.log2(L)
    bcap = radlim = None
    if v2:
        bcap = math.ceil(n * n / math.log2(n))
        radlim = math.ceil(n / math.log2(n))
    tb, cov, prp, lns = [], [], [], []
    with open(path) as f, open(path + ".fac") as g:
        for k, (line, fl) in enumerate(zip(f, g)):
            if limit and k >= limit:
                break
            w = line.split()
            T = int(w[2])
            sfac = {}
            for tok in fl.split():
                q, e = tok.split("^")
                sfac[int(q)] = int(e)
            tb.append(math.log2(T) if T > 1 else 0.0)
            cov.append(coverage_cells(sfac, bcap=bcap, radlim=radlim))
            prp.append(w[4] == "1")
            lns.append(math.log(max(int(w[3]), 1)))
    return Pool(tb, np.array(cov), lam, prp, lns)


if __name__ == "__main__":
    print(__doc__)
    sys.exit(0)
