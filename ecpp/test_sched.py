#!/usr/bin/env python3
"""test_sched.py -- the stopping-rule scheduler (shortsim.IndexTable, shortsched).

  A  the break-even table against the reference dynamic programme: for several
     rewards R, "r(c, j, t) <= R" must be solve_class's go/no-go (up to the
     states within one grid step of the boundary), and the fresh value V0 must
     match its value function;
  B  monotonicity: more coverage never makes a state worse, and the value of a
     fresh tail grows with the reward;
  C  the coverage claim behind the classes: for random smooth parts s the
     divisors m of s cover every delta in (0, log2 s), i.e. coverage_cells(s)
     is an initial segment of total length log2 s;
  D  queue mechanics on a toy ladder with a hand-made table: order of work,
     the stopping rule, chunking, sweeps of one B1 in flight together but no
     larger B1 before the smaller has reported, cancel; with an engine that
     runs 8 curves at a time (yafu) every unit is a multiple of 8;
  E  simulated levels end to end (virtual time, tails from the model): they
     terminate with a win, the work done is of the size the model predicts and
     the cores are kept busy; a winner whose build fails resumes the search; a
     scan that cannot go further ends in a win or a clean failure;
  F  the yafu engine: its report is parsed into the same splits as gmp-ecm's
     (on recorded output, and against gmp-ecm on planted factors when a yafu
     binary is given in SHORTECPP_YAFU), and a simulated level with its tiers,
     costs and success table ends with a win.

usage: test_sched.py            (about 30 s)
"""
import math
import os
import random
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import shortsim as ss            # noqa: E402
import shortsched as sd          # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print("  %-62s %s%s" % (name, "ok" if ok else "FAIL", ("  " + detail) if detail else ""))
    if not ok:
        FAILS.append(name)


def make_ladder(bits):
    ss.set_grid(bits)
    cost = ss.make_cost(ss.SEC1030_LOADED)
    tiers = sd.default_tiers()
    lad = ss.Ladder(tiers, ss.EcmModel(os.path.join(HERE, "ecm_pi.json")), y0bits=26, cost=cost)
    return tiers, lad, cost


def test_table():
    print("A/B: break-even table (768-bit level)")
    tiers, lad, cost = make_ladder(769)
    lam = 385.0
    tab = ss.IndexTable(lad, lam, cmax=60, threads=4)
    J = len(tiers)
    ratio = tab.rgrid[1] / tab.rgrid[0]
    t = np.arange(ss.NB)
    for R in (900.0, 4000.0, 20000.0):
        bad = tot = worst = 0
        for c in (3, 12, 25, 40, 60):
            V, go = ss.solve_class(lad, ss.class_window(c, lam), 0, R)
            live = t >= lam - c                  # tails that can still hold a q of the window
            for j in range(J):
                ok = (lad.Z[j] >= 1e-200) & live
                a = tab.rmin[c, j] <= R
                diff = (a != go[j]) & ok
                tot += int(go[j][ok].sum())
                bad += int(diff.sum())
                # a disagreement is allowed only where r is within one grid step of R
                far = diff & ((tab.rmin[c, j] < R / ratio) | (tab.rmin[c, j] > R * ratio))
                worst += int(far.sum())
        check("R = %-6.0f go/no-go agrees away from the boundary" % R, worst == 0 and bad <= 0.03 * tot,
              "%d go states, %d boundary disagreements, %d elsewhere" % (tot, bad, worst))
    # V0 at a grid reward equals the reference value function
    g = len(tab.rgrid) // 2
    R = float(tab.rgrid[g])
    err = 0.0
    for c in (5, 20, 50):
        V, _ = ss.solve_class(lad, ss.class_window(c, lam), 0, R)
        ok = (lad.Z[0] >= 1e-200) & (t >= lam - c)
        err = max(err, float(np.max(np.abs(tab.V0[g, c][ok] - V[0][ok]) / np.maximum(1.0, np.abs(V[0][ok])))))
    check("fresh value V0 equals the reference at a grid reward", err < 1e-3, "max rel. error %.1e" % err)
    # monotone in the class at fixed sweep and size, where both are finite
    viol = 0
    for j in (0, 5, 12):
        a = tab.rmin[10:59, j, 380:]             # sizes alive in every class compared
        b = tab.rmin[11:60, j, 380:]
        viol += int((b > a * 1.0001).sum())
    check("more coverage never raises the break-even reward", viol == 0, "%d violations" % viol)
    F, E = tab.fresh(np.ones((tab.cmax + 1, ss.NB)))
    check("fresh value and ECM cost grow with the reward", bool(np.all(np.diff(F) >= -1e-6) and np.all(np.diff(E) >= -1e-6)))
    R1 = tab.threshold(0.01, F / F.size * 0 + F / 1e6)
    R2 = tab.threshold(0.10, F / 1e6)
    check("a dearer supply raises the rule's reward", R2 > R1, "%.0f -> %.0f" % (R1, R2))
    return tiers, lad, cost, tab, lam


def test_coverage():
    print("C: divisors of s cover (0, log2 s)")
    rng = random.Random(5)
    primes = [q for q in range(2, 2000) if all(q % d for d in range(2, int(q ** 0.5) + 1))]
    worst = 0.0
    gaps = 0
    for _ in range(400):
        sfac = {}
        for _ in range(rng.randint(1, 6)):
            q = rng.choice(primes[:rng.choice((5, 30, len(primes)))])
            sfac[q] = sfac.get(q, 0) + rng.randint(1, 3)
        ls = sum(e * math.log2(q) for q, e in sfac.items())
        if ls > 150:
            continue
        cells = ss.coverage_cells(sfac)
        worst = max(worst, abs(float(cells.sum()) - ls))
        k = int(ls)
        if k >= 3 and not np.all(cells[1:k - 1] > 0.999):
            gaps += 1
    check("coverage is an initial segment of length log2 s", gaps == 0 and worst < 0.26,
          "max |measure - log2 s| = %.3f bits (cell resolution 1/8), %d with gaps" % (worst, gaps))


class FakeTable:
    """r = base[class] * (1 + sweep) -- later sweeps and lower classes are worse."""
    cmax = 9

    def __init__(self, inf_class=None):
        self.inf_class = inf_class

    def r(self, cls, j, tb):
        if cls == self.inf_class:
            return float("inf")
        return (10 - cls) * 100.0 * (1 + j) + tb * 1e-3

    def pwin(self, cls, j, tb):
        return 0.0


class C:
    def __init__(self, cls, tb=400):
        self.dead = False
        sd.sched_init(self, cls, tb)


def test_queue():
    print("D: queue mechanics")
    tiers = [("ecm", 2000, 4), ("ecm", 50000, 6), ("ecm", 50000, 6), ("ecm", 1000000, 8)]
    secs = {2000: 0.01, 50000: 1.0, 1000000: 10.0}
    sch = sd.Sched(FakeTable(), tiers, 8, lambda b1, tb: secs[b1], unit_min=1.5, chunk_secs=1.5)
    check("stages group the sweeps of one B1", sch.stages == [(0, [4]), (1, [6, 12]), (3, [8])], str(sch.stages))
    cs = [C(k) for k in (9, 9, 5, 5, 1)]
    for c in cs:
        sch.push(c)
    u = sch.take(1, 250.0)
    check("best first, cheap sweeps share a process", len(u) == 1 and u[0].cands == cs[:2] and u[0].n == 4 and u[0].j == 0,
          "%d cands, %d curves" % (len(u[0].cands), u[0].n))
    check("stopping rule holds the rest back", sch.take(4, 250.0) == [] and sch.best() == 500.4)
    sch.finished(u[0])
    check("stage advances when its curves are in", cs[0].st == 1 and cs[0].sent == 0 and cs[0].done == 0)
    u2 = sch.take(8, 250.0)
    # stage 1: 6 + 6 curves at 1 s each, unit of 1.5 s -> 1 curve per unit; sweep 1 has r = 200, sweep 2 r = 300
    n1 = [x for x in u2 if x.cands == [cs[0]]]
    check("an expensive sweep is split over the cores", len(u2) == 8 and all(x.n == 1 for x in u2)
          and len(n1) == 4 and all(x.j == 1 for x in u2), "%d units, %d on the first candidate" % (len(u2), len(n1)))
    u3 = sch.take(8, 250.0)
    check("... up to the sweep's curves, the next sweep only if the rule allows",
          len(u3) == 4 and cs[0].sent == 6 and cs[1].sent == 6 and sch.take(8, 250.0) == [])
    u4 = sch.take(3, 350.0)
    check("sweeps of one B1 run together when the rule allows", len(u4) == 3 and all(x.j == 2 for x in u4))
    for x in u2 + u3:
        sch.finished(x)
    check("a larger B1 waits for the smaller", cs[0].st == 1 and cs[0].done == 6 and sch.take(8, 1e9)[0].j in (0, 2))
    before = (cs[0].sent, cs[0].infl)
    for x in u4:
        sch.cancel(x)
    check("cancel returns the curves", cs[0].infl + cs[1].infl <= before[1] + cs[1].infl and cs[0].sent + cs[1].sent <= 24)
    sch2 = sd.Sched(FakeTable(inf_class=3), tiers, 4, lambda b1, tb: secs[b1])
    d = C(3)
    sch2.push(d)
    check("a state no reward makes worthwhile is abandoned", d.dead and sch2.nabandoned == 1)
    # an engine that runs 8 curves at a time
    t8 = [("ecm", 2000, 16), ("ecm", 50000, 16), ("ecm", 50000, 16), ("ecm", 1000000, 32)]
    sch3 = sd.Sched(FakeTable(), t8, 8, lambda b1, tb: secs[b1], unit_min=1.5, chunk_secs=1.5, quantum=8)
    e = [C(9), C(9)]
    for c in e:
        sch3.push(c)
    seen = []
    for _ in range(40):
        us = sch3.take(8, 1e9)
        if not us:
            break
        seen += us
        for x in us:
            sch3.finished(x)
    check("with quantum 8 every unit is a multiple of 8 curves", bool(seen) and all(x.n % 8 == 0 and x.n >= 8 for x in seen)
          and all(c.st == 3 for c in e), "%d units, sizes %s" % (len(seen), sorted({x.n for x in seen})))
    try:
        sd.Sched(FakeTable(), tiers, 8, lambda b1, tb: 1.0, quantum=8)
        bad = False
    except ValueError:
        bad = True
    check("... and sweeps that are not are refused", bad)


def test_sim(tiers, lad, cost, tab, lam):
    print("E: a simulated level")

    class Pool:
        pass
    rng = np.random.default_rng(11)
    n = 60000
    pool = Pool()
    # tail sizes and smooth parts with roughly the law of real 768-bit orders
    pool.w = np.minimum(rng.exponential(13.0, size=n) + 1.0, 90.0)
    pool.t = 769.0 - pool.w - rng.exponential(22.0, size=n)
    pool.prp = rng.random(n) < 0.06
    pool.lam = lam
    walls, orders, work = [], [], []
    for k in range(12):
        drv = sd.SimDriver(pool, lad, tab, 64, np.random.default_rng(100 + k), cost)
        st = sd.level_loop(drv, tab, tiers, cost)
        walls.append(st["wall"])
        orders.append(st["orders"])
        work.append(st["scan_s"] + st["ecm_s"])
        if not (st["orders"] > 0 and st["wall"] > 0 and drv.wins >= 1):
            check("level %d ends with a win" % k, False)
    mw = float(np.mean(work))
    check("12 simulated levels end with a win", len(walls) == 12, "mean %.0f core-s, %.0f orders, %.0f s wall on 64 threads"
          % (mw, float(np.mean(orders)), float(np.mean(walls))))
    # the rule's own estimate of the cost per win at the mean final price is of the same size
    check("work per level is of the modelled size", 500 < mw < 60000, "%.0f core-s" % mw)
    util = sum(work) / (sum(walls) * 64)
    check("cores are kept busy", util > 0.7, "utilisation %.0f%%" % (100 * util))

    # a winner whose build fails: the search must resume and end at a later win
    class FailOnce(sd.SimDriver):
        def __init__(self, *a, **k):
            sd.SimDriver.__init__(self, *a, **k)
            self.builds = 0
            self.first = None

        def build(self, stats):
            self.builds += 1
            if self.builds == 1:
                self.first = dict(stats)
                self.wins = 0
                return None
            return stats
    drv = FailOnce(pool, lad, tab, 64, np.random.default_rng(7), cost)
    st = sd.level_loop(drv, tab, tiers, cost)
    check("a failed build resumes the search", drv.builds == 2 and st["orders"] >= drv.first["orders"]
          and st["ecm_s"] > drv.first["ecm_s"] and st["wall"] > drv.first["wall"],
          "first win after %.0f core-s, second after %.0f" % (drv.first["ecm_s"], st["ecm_s"]))

    # a scan that cannot go further: everything queued is then in play, and the level
    # either wins or reports that it is exhausted
    ends = []
    for k in range(6):
        drv = sd.SimDriver(pool, lad, tab, 64, np.random.default_rng(300 + k), cost, B0=1e6, Bmax=1e6)
        try:
            st = sd.level_loop(drv, tab, tiers, cost)
            ends.append("win" if drv.wins else "?")
        except RuntimeError as e:
            ends.append("exhausted" if "exhausted" in str(e) else "?")
    check("an exhausted scan ends in a win or a clean failure", all(x in ("win", "exhausted") for x in ends)
          and "exhausted" in ends, " ".join(ends))


YAFU_SAMPLE = """
>> === Starting work on batchfile expression #1 ===
ecm(%(n1)d,24)
=============================================
ecm: 8/24 curves on C312 @ B1=11000, B2=1446869\recm: 8/24 curves on C312 @ B1=11000, B2=1446869
ecm: process took 0.1561 seconds.
ecm: found c22 factor = %(f23)d
ecm: found prp11 factor = %(f1)d

***factors found***
C22 = %(f23)d
PRP11 = %(f1)d
C280 = %(g)d

ans = %(g)d

=== Starting work on batchfile expression #2 ===
ecm(%(n2)d,16)
***factors found***
C300 = %(n2)d

ans = %(n2)d

=== Starting work on batchfile expression #3 ===
ecm(%(n3)d,16)
ecm: found prp10 factor = %(f1)d
ecm: found prp297 (co)factor = %(h)d

***factors found***
PRP10 = %(f1)d
PRP297 = %(h)d

ans = 1

=== Starting work on batchfile expression #4 ===
ecm(%(n4)d,16)
***factors found***
PRP10 = %(f1)d
C200 = %(g)d
eof; done processing batchfile
"""


def test_yafu(tiers_gmp, lam):
    print("F: the yafu engine")
    import shortECPP as sp
    f1, f2, f3 = 20699021657, 1546280705083, 2920621559
    g = (2 ** 521 - 1) * (2 ** 607 - 1)
    h = 2 ** 607 - 1
    v = {"f1": f1, "f23": f2 * f3, "g": g, "h": h, "n1": f1 * f2 * f3 * g, "n2": g * 3 + 2, "n3": f1 * h,
         "n4": f1 * g * 7}
    nums = [v["n1"], v["n2"], v["n3"], v["n4"]]
    res = sp.parse_yafu(YAFU_SAMPLE % v, nums)
    check("yafu report: factors, composite cofactor", res.get(v["n1"]) == ([f2 * f3, f1], g, False))
    check("yafu report: nothing found is no split", v["n2"] not in res)
    check("yafu report: probable-prime cofactor", res.get(v["n3"]) == ([f1], h, True))
    check("yafu report: parts that do not multiply back are dropped", v["n4"] not in res and len(res) == 2)
    ybin = os.environ.get("SHORTECPP_YAFU", "")
    if ybin and os.access(ybin, os.X_OK) and sp.ECM_BIN:
        sp.YAFU_BIN = ybin
        rng = random.Random(5)

        def prime(bits):
            while True:
                n = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
                if sp.is_prp(n, 8):
                    return n
        big = prime(400) * prime(420)
        tails = [prime(30) * prime(33) * big, prime(31) * prime(700), big, prime(28) * prime(36) * prime(34) * big]
        ry = sp.EcmRunner("yafu").run(tails, 11000, 24, 300)
        rg = sp.EcmRunner("gmp").run(tails, 11000, 24, 300)
        same = all((t in ry) == (t in rg) and (t not in ry or (ry[t][1], ry[t][2]) == (rg[t][1], rg[t][2])) for t in tails)
        check("yafu and gmp-ecm split planted factors alike", same and len(ry) == 3,
              "%d and %d of 4 split" % (len(ry), len(rg)))
    else:
        print("  (no yafu binary in SHORTECPP_YAFU: live comparison skipped)")
    if not os.path.exists(os.path.join(HERE, ss.ENGINES["yafu"]["pi"])):
        print("  (no yafu success table: simulated level skipped)")
        return
    em, sec, cost, quantum = ss.engine_model("yafu", 1.0, HERE)
    cost1 = ss.engine_model("yafu", 1.0, HERE, scalar=True)[2]
    grid = [(b1, tb) for b1 in sec for tb in (30, 64, 300, 1030, 1040, 1041, 1300)]
    same = all(abs(cost1(b1, tb) - float(cost(b1, max(tb, 64)))) <= 1e-12 * cost1(b1, tb) for b1, tb in grid)
    st = ss.ENGINES["yafu"]["step"]
    check("cost of a yafu curve: scalar and array forms agree, step above %d bits" % st[0],
          same and abs(cost1(50000, 1030) - sec[50000]) < 1e-12 and cost1(50000, st[0] + 1) > st[1] * cost1(50000, st[0]) * 0.999
          and quantum == 8)
    tiers = sd.default_tiers("yafu")
    lad = ss.Ladder(tiers, em, y0bits=26, cost=cost)
    tab = ss.IndexTable(lad, lam, cmax=60, threads=4)
    py = em.pi("ecm", 50000, np.array([50.0, 70.0, 90.0]))
    check("yafu success table loads", bool(np.all(py > 0) and np.all(np.diff(py) < 0)), "pi(50, 70, 90 bits at B1=5e4) = %s" % np.round(py, 5))

    class Pool:
        pass
    rng = np.random.default_rng(12)
    n = 60000
    pool = Pool()
    pool.w = np.minimum(rng.exponential(13.0, size=n) + 1.0, 90.0)
    pool.t = 769.0 - pool.w - rng.exponential(22.0, size=n)
    pool.prp = rng.random(n) < 0.06
    pool.lam = lam
    work, sizes = [], set()
    for k in range(8):
        drv = sd.SimDriver(pool, lad, tab, 64, np.random.default_rng(300 + k), cost)
        start = drv.start_unit

        def spy(u, start=start):
            sizes.add(u.n % 8)
            start(u)
        drv.start_unit = spy
        st = sd.level_loop(drv, tab, tiers, cost, {"quantum": quantum})
        work.append(st["scan_s"] + st["ecm_s"])
    check("8 simulated levels with yafu end with a win, units in batches of 8", len(work) == 8 and sizes == {0},
          "mean %.0f core-s" % (sum(work) / len(work)))


def main():
    tiers, lad, cost, tab, lam = test_table()
    test_coverage()
    test_queue()
    test_sim(tiers, lad, cost, tab, lam)
    test_yafu(tiers, lam)
    print("ALL OK" if not FAILS else "FAILED: %s" % ", ".join(FAILS))
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
