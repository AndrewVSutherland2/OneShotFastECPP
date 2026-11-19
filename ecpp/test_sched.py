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
import time

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
    # candidates that are each a fair chance do not share a process
    class Hot(FakeTable):
        def pwin(self, cls, j, tb):
            return 0.03
    sch4 = sd.Sched(Hot(), tiers, 2, lambda b1, tb: secs[b1], unit_min=1.5, chunk_secs=1.5, unit_pw=0.05)
    for c in [C(9) for _ in range(40)]:
        sch4.push(c)
    u = sch4.take(1, 1e9)
    sch5 = sd.Sched(FakeTable(), tiers, 2, lambda b1, tb: secs[b1], unit_min=1.5, chunk_secs=1.5, unit_pw=0.05)
    for c in [C(9) for _ in range(40)]:
        sch5.push(c)
    v = sch5.take(1, 1e9)
    check("a shared process stops at unit_pw expected wins", len(u[0].cands) == 2 and len(v[0].cands) > 2,
          "%d candidates at 0.03 each, %d when wins are remote" % (len(u[0].cands), len(v[0].cands)))
    # a search taken up with the other ECM engine: the two ladders at 1362 bits
    ty, tg = sd.default_tiers("yafu", 1362), sd.default_tiers("gmp", 1362)
    sy, sg = sd.ladder_stages(ty), sd.ladder_stages(tg)
    qg = sd.Sched(FakeTable(), tg, 8, lambda b1, tb: 1.0)
    x = C(9)
    x.st, x.sent, x.done = 7, 64, 64        # 64 of yafu's 192 curves at B1 = 3e7, where gmp-ecm's ladder has 48
    qg.push(x)
    limbo = not x.dead and qg.seq == 0      # read as it stands it is neither queued nor through
    ex = (sd.restage(7, 64, ty, tg), sd.restage(7, 32, ty, tg), sd.restage(4, 40, ty, tg), sd.restage(4, 72, ty, tg),
          sd.restage(3, 47, tg, ty, 8), sd.restage(8, 0, ty, tg),
          sd.restage(0, 8, [("ecm", 3000000, 16)], [("ecm", 2000, 8), ("ecm", 10000000, 8)]))
    lost = 0
    for (old, new, q) in ((ty, tg, 1), (tg, ty, 8)):
        so, sn = sd.ladder_stages(old), sd.ladder_stages(new)
        sch6 = sd.Sched(FakeTable(), new, 8, lambda b1, tb: 1.0, quantum=q)
        for st in range(len(so) + 1):
            for done in range(so[st][2] if st < len(so) else 1):
                y = C(9)
                y.st, y.done = sd.restage(st, done, old, new, q)
                y.sent = y.done
                n0 = sch6.seq
                sch6.push(y)
                lost += not (y.dead or sch6.seq == n0 + 1) or y.done % q != 0 or y.st > len(sn)
    check("a candidate's place on one engine's ladder is read onto the other's: queued or through, never neither",
          limbo and lost == 0 and ex == ((8, 0), (7, 32), (4, 40), (5, 0), (3, 40), (8, 0), (1, 0))
          and [t[2] for t in sy] == [16, 24, 32, 48, 80, 96, 96, 192] and [t[2] for t in sg] == [12, 20, 32, 48, 72, 96, 96, 48],
          "%d lost; %s" % (lost, ex))


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

    # a search that is cut off and taken up again from its checkpoint: the counters
    # and every candidate's place on the ladder carry over, what was in flight is run again
    class Stop(Exception):
        pass

    class Ck(sd.SimDriver):
        ckpt_every = 4.0                    # virtual seconds

        def __init__(self, *a, **k):
            sd.SimDriver.__init__(self, *a, **k)
            self.all, self.saved, self.given, self.nck, self.stop_at = [], None, None, 0, 0

        def take_orders(self, payload):
            out, n, cst = sd.SimDriver.take_orders(self, payload)
            self.all += out
            return out, n, cst

        def checkpoint(self, loop):
            self.nck += 1
            if self.nck == self.stop_at:
                self.saved = (loop, [c for c in self.all if not c.dead], self.B, self.Bnext, self.t)
                raise Stop()

        def restore(self):
            g, self.given = self.given, None
            return g
    drv = None
    for seed in range(41, 60):              # a level that lasts until its third checkpoint
        drv = Ck(pool, lad, tab, 64, np.random.default_rng(seed), cost)
        drv.stop_at = 3
        try:
            sd.level_loop(drv, tab, tiers, cost)
        except Stop:
            break
    ok = drv.saved is not None
    detail = "no level lasted three checkpoints"
    if ok:
        loop, cands, B, Bnext, t = drv.saved
        infl = sum(1 for c in cands if c.infl)
        for c in cands:
            c.sent, c.infl = c.done, 0
        d2 = Ck(pool, lad, tab, 64, np.random.default_rng(77), cost)
        d2.B, d2.Bnext, d2.t, d2.all, d2.given = B, Bnext, t, list(cands), (cands, loop)
        st = sd.level_loop(d2, tab, tiers, cost)
        ok = (d2.wins >= 1 and st["orders"] >= loop["K"] > 0 and st["ecm_s"] > loop["spent_ecm"] > 0
              and st["scan_s"] >= loop["spent_scan"] and st["wall"] > loop["elapsed"] and st["expwin"] > loop["expwin"] > 0
              and st["units"] > loop["nunits"])
        detail = ("cut at %.0f s with %d orders, %d live candidates (%d in flight), %.2f expected wins; won at %.0f s with %d orders, %.2f"
                  % (loop["elapsed"], loop["K"], len(cands), infl, loop["expwin"], st["wall"], st["orders"], st["expwin"]))
    check("a search cut off resumes from its checkpoint", ok, detail)


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
    check("yafu report: an input dropped after its first batch says how many curves ran",
          res.ran == {v["n1"]: 8}, str({k.bit_length(): x for k, x in res.ran.items()}))
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
        # yafu alone stops at the first batch with a factor; the runner must spend the rest of the
        # sweep on the cofactor: six numbers with a 22- and a 44-bit prime, 64 curves at B1 = 2000
        # (complete sweeps find all twelve primes 99.9% of the time, first batches only 15%)
        two = [prime(22) * prime(44) * big for _ in range(6)]
        r2 = sp.EcmRunner("yafu").run(two, 2000, 64, 300)
        check("a yafu sweep goes on with the cofactor after a factor", all(t in r2 and r2[t][1] == big for t in two),
              "%d of 6 fully stripped" % sum(1 for t in two if t in r2 and r2[t][1] == big))
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
    # thin cells are pooled: no single cell of a dozen hits sets the rate beyond the measured sizes
    deep = np.arange(100.0, 172.0, 4.0)
    rat = [em.pi("ecm", b1, deep) / ss.ecm_pi_theory(deep, b1, ss.B2_YAFU.get(b1), em.lndelta) for b1 in (3000000, 10000000, 30000000)]
    check("deep success rates stay near the two-stage Dickman curve", all(0.8 < r.min() and r.max() < 1.25 for r in rat),
          "measured/theory at B1 = 3e6, 1e7, 3e7: %s" % ", ".join("%.2f..%.2f" % (r.min(), r.max()) for r in rat))
    cells = [{"method": "ecm", "b1": 50000, "bits": 40 + 4 * i, "curves": 4000,
              "hits": h} for i, h in enumerate((2000, 900, 400, 150, 60, 20, 9, 1, 0, 0))]
    cells[6]["hits"] = 27                   # one thin cell three times too high
    m2 = ss.EcmModel(b2=ss.B2_YAFU)
    m2.cells[("ecm", 50000)] = cells
    x = np.array([64.0, 76.0])
    t2 = ss.ecm_pi_theory(np.array([c["bits"] for c in cells], dtype=float), 50000, ss.B2_YAFU.get(50000), m2.lndelta)
    pooled = float(sum(c["hits"] for c in cells[4:]) / (4000 * t2[4:].sum()))
    got = float(m2.pi("ecm", 50000, x)[1] / ss.ecm_pi_theory(x, 50000, ss.B2_YAFU.get(50000), m2.lndelta)[1])
    alone = float(27 / (4000 * t2[6]))
    check("a thin cell does not set the rate beyond the measured sizes", abs(got / pooled - 1.0) < 0.35 and got < 0.7 * alone,
          "ratio to theory beyond the cells %.2f; pooled thin cells %.2f; the high cell alone %.2f" % (got, pooled, alone))

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


def test_small():
    print("G: small levels")
    check("medium-prime bound by size", [ss.medium_bits(b) for b in (140, 300, 432, 485, 769, 1330)] == [48, 99, 143, 160, 160, 160])
    tops = [max(b1 for (_, b1, _) in sd.default_tiers("yafu", b)) for b in (200, 300, 400, 499, 500)]
    deep = [sum(c for (_, b1, c) in sd.default_tiers("yafu", b) if b1 == 30000000) for b in (769, 1299, 1362)]
    check("the ladder stops earlier at small sizes, and is longer at 3e7 from 1300 bits",
          tops == [11000, 50000, 250000, 1000000, 30000000] and deep == [64, 64, 192]
          and sd.default_tiers("gmp", 769) == sd.default_tiers("gmp") and sd.default_tiers("yafu", 769) == sd.default_tiers("yafu"),
          "%s, %s curves at 3e7" % (tops, deep))
    import shortECPP as sp
    check("level parameters by size", sp.level_params(167, 192, 3778) == (16, 10 ** 5, 1)
          and sp.level_params(300, 192, 10938) == (20, 10 ** 6, 8) and sp.level_params(769, 192, 61752)[1:] == (None, 192)
          and sp.level_params(210, 32, 1 << 21)[0] == 22)
    check("smallest CM level follows the recursion floor", sp.cm_min_bits(3778 ** 2) == 72 and sp.cm_min_bits(170455 ** 2) == 94)
    bits, lam = 300, 150.0
    ss.set_grid(bits)
    em, sec, cost, quantum = ss.engine_model("yafu", 1.0, HERE)
    tiers = sd.default_tiers("yafu", bits)
    lad = ss.Ladder(tiers, em, y0bits=20, cost=cost)
    scale = (bits / 500.0) ** 8
    tab = ss.IndexTable(lad, lam, rlo=3.0 * scale, rhi=3.0e5 * scale, threads=2)
    fin = np.isfinite(tab.rmin[5:40, 0, 200:290])
    check("300 bits: grid and D shrink, first sweeps have a price", ss.NB == 512 and ss.D == 99 and fin.mean() > 0.9
          and float(np.nanmin(np.where(fin, tab.rmin[5:40, 0, 200:290], np.nan))) < 5.0,
          "cheapest first sweep %.2f core-s per win" % float(np.nanmin(np.where(fin, tab.rmin[5:40, 0, 200:290], np.nan))))

    class Pool:
        pass
    rng = np.random.default_rng(3)
    n = 4000
    pool = Pool()
    pool.w = np.minimum(rng.exponential(13.0, size=n) + 1.0, 45.0)
    pool.t = 300.0 - 20.0 - rng.exponential(12.0, size=n)
    pool.prp = rng.random(n) < 0.06
    pool.lam = lam
    work = []
    for k in range(6):
        drv = sd.SimDriver(pool, lad, tab, 8, np.random.default_rng(700 + k), ss.engine_model("yafu", 1.0, HERE, scalar=True)[2])
        st = sd.level_loop(drv, tab, tiers, cost, {"quantum": quantum})
        work.append(st["scan_s"] + st["ecm_s"])
    check("6 simulated 300-bit levels end with a win", len(work) == 6 and max(work) < 2000, "mean %.1f core-s" % (sum(work) / len(work)))
    # a winner in reserve: the loop goes on until the driver says build
    drv = sd.SimDriver(pool, lad, tab, 8, np.random.default_rng(1), ss.engine_model("yafu", 1.0, HERE, scalar=True)[2])
    calls = []

    def hold(R, drv=drv):
        calls.append(drv.wins)
        return len(calls) > 25
    drv.should_build = hold
    st = sd.level_loop(drv, tab, tiers, cost, {"quantum": quantum})
    drv0 = sd.SimDriver(pool, lad, tab, 8, np.random.default_rng(1), ss.engine_model("yafu", 1.0, HERE, scalar=True)[2])
    st0 = sd.level_loop(drv0, tab, tiers, cost, {"quantum": quantum})
    check("a winner kept in reserve: the search goes on, the first win's figures are kept",
          len(calls) == 26 and st["wall"] > st0["wall"] and abs(st["first_win_wall"] - st0["wall"]) < 1e-9
          and st["expwin"] <= st0["expwin"] + 1e-9 and calls[-1] >= calls[0] >= 1,
          "%d more wins while waiting, %.1f s against %.1f s" % (calls[-1] - calls[0], st["wall"], st0["wall"]))
    # the last scan chunk brings a winner and nothing to work on, and the driver
    # would keep it in reserve: with no search left it is built all the same
    class Last:
        threads = 4

        def __init__(self, good):
            self.t, self.scanned, self.won, self.built, self.asked, self.good, self.ev = 0.0, False, False, 0, 0, good, []

        def now(self):
            return self.t

        def log(self, m):
            pass

        def next_chunk(self):
            return None if self.scanned else (1.0, 1.0)

        def start_scan(self, ns):
            self.scanned = True
            self.ev.append(("orders", None))

        def wait(self, timeout):
            self.t += 1.0
            ev, self.ev = self.ev, []
            return ev

        def take_orders(self, payload):
            self.won = True
            return [], 1, 1.0

        def has_winner(self):
            return self.won

        def abort(self):
            pass

        def should_build(self, R):
            self.asked += 1
            return False

        def build(self, st):
            self.built += 1
            self.won = False
            return st if self.good else None
    last = Last(True)
    st = sd.level_loop(last, tab, tiers, cost, {"quantum": quantum})
    bad = Last(False)
    try:
        sd.level_loop(bad, tab, tiers, cost, {"quantum": quantum})
        msg = "returned"
    except RuntimeError as e:
        msg = str(e)
    check("a winner in reserve is built when the scan is exhausted and no work is left",
          last.built == 1 and last.asked == 1 and st["orders"] == 1 and bad.built == 1 and msg.startswith("scan exhausted"),
          "built %d, asked %d; after a failed build: %s" % (last.built, last.asked, msg))
    ss.set_grid(1527)
    # build estimates and the accept-or-pass step of the prover's driver
    est = sp.build_estimate([450490639, 140861442439, 4315], 769)
    check("build estimate from the class number", len(est) == 3 and 40 < est[450490639] < 60 and est[4315] < 1.0
          and est[140861442439] > 100, "%.1f s for h = 2 x 6379 at 769 bits" % est.get(450490639, 0))
    sp.BUILD_EST.clear()

    class W:
        def __init__(self, d):
            self.d = d

    class D:
        pass
    d = D()
    d.winners, d.nb, d.threads, d.build_ref, d.pass_t0, d.pass_n = [(0, W(450490639), 0, [])], 769, 192, 2.5, None, 0
    quiet, sp.log = sp.log, (lambda m: None)
    hold1 = sp.RuleDriver.should_build(d, 1600.0)          # 51 s against 8 + 2.5 s
    d.winners.append((0, W(4315), 0, []))
    hold2 = sp.RuleDriver.should_build(d, 1600.0)          # a cheap winner has arrived
    d2 = D()
    d2.winners, d2.nb, d2.threads, d2.build_ref, d2.pass_t0, d2.pass_n = [(0, W(450490639), 0, [])], 769, 32, 2.5, None, 0
    hold3 = sp.RuleDriver.should_build(d2, 1600.0)         # on 32 threads the next winner is 50 s away
    d.winners = [(0, W(450490639), 0, [])]
    d.pass_t0 = time.time() - 60.0
    hold4 = sp.RuleDriver.should_build(d, 1600.0)          # waited as long as the build would have taken
    sp.log = quiet
    sp.BUILD_EST.clear()
    check("a costly winner waits for a cheap one, but not on few cores or for long", (hold1, hold2, hold3, hold4) == (False, True, True, True),
          str((hold1, hold2, hold3, hold4)))
    # checkpoint record of a candidate: one large integer, everything else rebuilt
    tt, vv, dd = 2 * 3 ** 80 + 1, 2 * 5 ** 50 + 1, 7
    pp = (tt * tt + dd * vv * vv) // 4
    c = sp.Cand(dd, tt, vv, pp + 1 + tt)
    c.s = c.N & -c.N
    c.mid, c.splits, c.prp_checked = 3 ** 0, 2, True
    c.tail = c.N // c.s
    sd.sched_init(c, 11, int(round(math.log2(c.tail))))
    c.st, c.sent, c.done, c.infl = 2, 48, 16, 2
    rec = sp.ckpt_pack(c, pp)
    e = sp.ckpt_unpack(rec, pp)
    same = all(getattr(c, k) == getattr(e, k) for k in ("d", "t", "v", "N", "s", "mid", "tail", "splits", "prp_checked", "cls", "tb", "st", "done"))
    bad = False
    try:
        sp.ckpt_unpack(rec, pp + 2)
    except ValueError:
        bad = True
    check("checkpoint record of a candidate: round trip, curves in flight are not counted, another prime is refused",
          same and e.sent == 16 and e.infl == 0 and not e.dead and bad and len(rec) == 11)
    # a checkpoint written with yafu and read by a run with gmp-ecm, and by one with yafu
    import shutil
    import tempfile
    tmpd = tempfile.mkdtemp(prefix="test_sched_ckpt_")
    keep = (sp.CKPT, sp.CKPT_RESUME, sp.log)
    sp.CKPT, sp.CKPT_RESUME, sp.log = os.path.join(tmpd, "x.search.pkl"), True, (lambda m: None)

    def stub(engine, nsm=10 ** 6):
        d = D()
        d.nsm = nsm
        d._CKPT_FIELDS = sp.RuleDriver._CKPT_FIELDS
        for k in d._CKPT_FIELDS:
            setattr(d, k, 0)
        d.act, d.est = {30000000: 5.0}, {30000000: 4.0}
        d.p, d.nb, d.stats, d.pool, d.engine = pp, 1362, {}, [], engine
        d.tiers, d.quantum = [tuple(t) for t in sd.default_tiers(engine, 1362)], ss.ENGINES[engine]["quantum"]
        return d
    w = stub("yafu")
    c.st, c.sent, c.done, c.dead = 7, 80, 64, False
    c2 = sp.ckpt_unpack(sp.ckpt_pack(c, pp), pp)
    c2.st, c2.sent, c2.done, c2.dead = 4, 40, 40, False
    w.pool = [c, c2]
    sp.RuleDriver.checkpoint(w, {"K": 2, "spent_scan": 1.0, "spent_ecm": 2.0, "elapsed": 3.0})
    rg, ry = stub("gmp"), stub("yafu")
    got_g = sp.RuleDriver._load_ckpt(rg)
    got_y = sp.RuleDriver._load_ckpt(ry)
    got_f = sp.RuleDriver._load_ckpt(stub("yafu", nsm=10 ** 5))      # the other certificate format
    sp.CKPT, sp.CKPT_RESUME, sp.log = keep
    shutil.rmtree(tmpd, True)
    check("a checkpoint taken up with the other engine: places converted; same engine: as they were; other format: refused",
          got_g is not None and got_y is not None and got_f is None
          and [(k.st, k.sent, k.done) for k in rg.pool] == [(8, 0, 0), (4, 40, 40)] and rg.act == {}
          and [(k.st, k.sent, k.done) for k in ry.pool] == [(7, 64, 64), (4, 40, 40)] and ry.act == {30000000: 5.0},
          "gmp-ecm run: %s" % [(k.st, k.done) for k in rg.pool])
    # supply: the wider bound or the higher cap, whichever is cheaper per order
    def supply(B, cap, n, cs, lo, hi, tot_n, tot_s, fixed=0, capmax=2 * 10 ** 10, Bmax=10 ** 15):
        d = D()
        d.nb, d.Bmax, d.grow, d.fixed_cap, d.capmax = 1362, Bmax, 2.0, fixed, capmax
        d.c_root, d.raise_fit, d._nxt, d.arm = 1.0e-9 * 1362 * 1362, (0.7, 1.0), None, None
        d._density = lambda b, c: sp.RuleDriver._density(d, b, c)
        d.B, d.cap, d.tot_orders, d.tot_test_s = B, cap, tot_n, tot_s
        d.last = (n, cs, 0.02 * cs, lo, hi, d._density(hi, cap))
        nx = sp.RuleDriver.next_chunk(d)
        return nx, d.arm
    # the scan of nextprime(10^410) as it stood at 10^9 and at 4 * 10^12 (cap 10^8 since 10^11)
    early = supply(10 ** 9, 10 ** 6, 4510, 311, 0, 10 ** 9, 4510, 311)
    late = supply(4096 * 10 ** 9, 10 ** 8, 68188, 286821, 2048 * 10 ** 9, 4096 * 10 ** 9, 263290, 629000)
    capped = supply(4096 * 10 ** 9, 10 ** 8, 68188, 286821, 2048 * 10 ** 9, 4096 * 10 ** 9, 263290, 629000, fixed=10 ** 8)
    atend = supply(10 ** 13, 10 ** 8, 68188, 286821, 5 * 10 ** 12, 10 ** 13, 263290, 629000, Bmax=10 ** 13)
    done = supply(10 ** 13, 10 ** 8, 68188, 286821, 5 * 10 ** 12, 10 ** 13, 263290, 629000, Bmax=10 ** 13, capmax=10 ** 8)
    check("more orders: a wider bound early, a higher cap once the bound is dear",
          early[1][0] == "widen" and late[1] == ("raise", 4096 * 10 ** 9, 4 * 10 ** 8) and capped[1][0] == "widen"
          and atend[1][0] == "raise" and done[0] is None,
          "at 4e12: %.0f orders for %.0f core-s by raising the cap" % late[0])


def main():
    tiers, lad, cost, tab, lam = test_table()
    test_coverage()
    test_queue()
    test_sim(tiers, lad, cost, tab, lam)
    test_yafu(tiers, lam)
    test_small()
    print("ALL OK" if not FAILS else "FAILED: %s" % ", ".join(FAILS))
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
