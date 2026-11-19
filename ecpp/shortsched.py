#!/usr/bin/env python3
"""shortsched.py -- the stopping-rule scheduler of one short-ECPP level.

The search for one level is a bandit problem: every candidate order is an arm
whose state is (coverage class, sweeps done, tail size), a pull is a few ECM
curves, the level ends at the first win, and fresh arms can be bought from the
discriminant scan at a price that rises as the scan widens.  shortsim.IndexTable
gives each state its break-even reward r (core-seconds per win at which the next
sweep pays under the best continuation); 1/r is the state's index.  The rule
implemented here:

  * work is done in increasing order of r, as many units at once as there are
    cores;
  * a candidate is worked only while r <= R, where R is the reward at which a
    fresh order is worth the scan's current marginal price -- the stopping rule;
    nothing is abandoned for good, R rises with the price and a parked candidate
    comes back when R reaches it;
  * the scan runs when the work below R is about to run out, on the share of the
    cores that its price has in the cost of a fresh order.  Cores left without
    work below R while it runs take what lies just above (up to relax x R, see
    level_loop): with a fast ECM program on many cores that is a large part of
    the ECM work (40-100% of it on 192 cores with yafu), but it is work brought
    forward, not wasted, since R rises with every chunk.  Giving the scan more
    cores instead buys the next, dearer chunk sooner and costs more per win
    (simulated at 768 bits on 192 threads: 0.78 against 0.65 core-hours).

Two things are needed to keep a many-core machine honest.  A sweep of c curves
is split into units of a few seconds (curves are independent trials, so a
candidate's sweep runs on c cores at once), and the sweeps of one B1 may be in
flight together when their r allows it: otherwise the eventual winner's deep
sweeps would run one after the other while every other core works through
candidates that the sequential rule would never have touched.  A larger B1 is
not started before the smaller ones have reported.

Sched holds the queues and knows nothing about processes; level_loop runs a
level against a driver, which is either the prover (shortECPP.py: gmp-ecm or
yafu, dscan) or SimDriver below (virtual time, tails drawn from the model), so
the same code is exercised by simulation and by the real thing.

The ECM program enters in three places only: the cost of a curve, the success
table behind the break-even table (shortsim.ENGINES), and the number of curves
it runs at once (quantum: 1 for gmp-ecm, 8 for yafu's AVX-512 ECM), of which
every unit is a multiple.  A unit must run all its curves, going on with the
cofactor when a factor turns up: the state "j sweeps done, tail of t bits" means
a tail that has survived those sweeps in full, and a tail that lost a small
prime to the first few curves of a sweep is not a fresh tail of its new size.

The log adds up the wins the model expects from the units handed out, in total
and by B1, and records the B1 of the unit that completed the first winner: at a
level's first win the total is exponentially distributed with mean 1 if the
model's win rates are right.

usage (simulation): shortsched.py sim <pool.pkl> [threads=192] [runs=200]
                        [seed=1] [bits=769] [engine=gmp|yafu] [cost=1.0]
                        [tiers=B1xCURVESxSWEEPS,...] [key=value ...]
"""

import bisect
import heapq
import math
import os
import sys
import time

INF = float("inf")


class Unit:
    """One process worth of work: n curves of sweep j on each of cands."""
    __slots__ = ("cands", "j", "method", "b1", "n", "est", "t0", "uid", "pw", "ps")

    def __init__(self, cands, j, method, b1, n, est, pw=0.0, ps=0.0):
        self.cands, self.j, self.method, self.b1, self.n, self.est = cands, j, method, b1, n, est
        self.t0 = 0.0
        self.uid = 0
        self.pw = pw          # wins the model expects from it
        self.ps = ps          # candidates it is expected to split


def sched_init(c, cls, tb):
    """Scheduler fields of a candidate object (any class with these attributes)."""
    c.cls, c.tb = cls, tb
    c.st = 0          # stage: index of the B1 it is at
    c.sent = 0        # curves of this stage handed to units
    c.done = 0        # curves of this stage reported
    c.infl = 0        # units in flight
    c.ver = 0         # bumps whenever its queue key changes (stale entries are skipped)


class Sched:
    def __init__(self, tab, tiers, threads, cost, unit_frac=0.05, unit_min=1.5, chunk_secs=1.5,
                 chunk_max=32, quantum=1, unit_pw=0.05):
        """tab: shortsim.IndexTable; tiers: [(method, B1, curves)] as given to the
        Ladder behind tab; cost(B1, tail bits): core-seconds per curve; quantum:
        curves the engine runs at once (8 for yafu's AVX-512 ECM), every unit is a
        multiple of it; unit_pw: candidates stop being added to a shared process
        once the model expects that many wins from it (at a small level the best
        few candidates are each a fair chance, and what runs alongside the unit
        that wins is wasted)."""
        self.unit_pw = unit_pw
        self.tab, self.tiers, self.threads, self.cost = tab, list(tiers), threads, cost
        self.quantum = quantum
        if any(c % quantum for (_, _, c) in self.tiers):
            raise ValueError("every sweep must be a multiple of %d curves" % quantum)
        self.unit_frac, self.unit_min = unit_frac, unit_min
        self.chunk_secs, self.chunk_max = chunk_secs, chunk_max
        self.stages = []                    # (first sweep, [cumulative curve counts of its sweeps])
        for j, (m, b1, c) in enumerate(self.tiers):
            if j and self.tiers[j - 1][:2] == (m, b1):
                self.stages[-1][1].append(self.stages[-1][1][-1] + c)
            else:
                self.stages.append((j, [c]))
        self.heaps = [[] for _ in self.tiers]
        self.seq = 0
        self.unit_secs = unit_min
        self.nabandoned = 0
        self.expwin = 0.0                   # wins the model expects from the units handed out
        self.expwin_b1 = {}                 # ... by B1
        self.expsplit_b1 = {}               # candidates the model expects to lose a factor to them, by B1

    def set_reward(self, R):
        """Units are kept short against the time the machine needs for one win, so
        that what is in flight when the level ends is a small part of its cost."""
        self.unit_secs = max(self.unit_min, self.unit_frac * R / self.threads)

    def position(self, c):
        """(sweep of the next curve to hand out, curves left in that sweep)."""
        first, ends = self.stages[c.st]
        k = bisect.bisect_right(ends, c.sent)
        return first + k, ends[k] - c.sent

    def push(self, c):
        """(Re)queue c at its current position and tail size."""
        c.ver += 1
        if c.dead:
            return
        if c.st >= len(self.stages):
            if c.infl == 0:
                c.dead = True               # through the whole ladder
                self.nabandoned += 1
            return
        if c.sent >= self.stages[c.st][1][-1]:
            return                          # this stage is all in flight: wait for it
        j, _ = self.position(c)
        r = self.tab.r(c.cls, j, c.tb)
        if r == INF:
            if c.infl == 0:
                c.dead = True               # no reward on the grid makes the next sweep pay
                self.nabandoned += 1
            return
        heapq.heappush(self.heaps[j], (r, self.seq, c.ver, c))
        self.seq += 1

    def _top(self, h):
        while h and (h[0][3].dead or h[0][2] != h[0][3].ver):
            heapq.heappop(h)
        return h[0] if h else None

    def best(self):
        """Smallest r in the queues (INF if nothing is queued)."""
        b = INF
        for h in self.heaps:
            e = self._top(h)
            if e is not None and e[0] < b:
                b = e[0]
        return b

    def take(self, nfree, rcap):
        """Up to nfree units, best first, none with r above rcap."""
        units = []
        while nfree > 0:
            bj, br = -1, INF
            for j, h in enumerate(self.heaps):
                e = self._top(h)
                if e is not None and e[0] < br:
                    bj, br = j, e[0]
            if bj < 0 or br > rcap:
                break
            h = self.heaps[bj]
            c = heapq.heappop(h)[3]
            method, b1, _ = self.tiers[bj]
            _, left = self.position(c)
            e1 = self.cost(b1, c.tb)
            if e1 * left <= self.unit_secs:
                # a whole sweep is cheap: share the process with the next best
                # candidates of the same sweep, but leave work for the other cores
                cands, est = [c], e1 * left
                kmax = max(1, min(self.chunk_max, len(h) // max(1, nfree)))
                pws = self.tab.pwin(c.cls, bj, c.tb)
                while est < self.chunk_secs and len(cands) < kmax and pws < self.unit_pw:
                    e = self._top(h)
                    if e is None or e[0] > rcap or self.position(e[3])[1] != left:
                        break
                    heapq.heappop(h)
                    cands.append(e[3])
                    est += self.cost(b1, e[3].tb) * left
                    pws += self.tab.pwin(e[3].cls, bj, e[3].tb)
                n = left
            else:
                q = self.quantum
                n = max(q, min(left, q * int(self.unit_secs / (q * e1))))
                cands, est = [c], e1 * n
            pw = n / float(self.tiers[bj][2]) * sum(self.tab.pwin(x.cls, bj, x.tb) for x in cands)
            self.expwin += pw
            self.expwin_b1[b1] = self.expwin_b1.get(b1, 0.0) + pw
            lad = getattr(self.tab, "lad", None)
            ps = 0.0
            if lad is not None:             # a second check of the model, on far more events than the wins
                ps = len(cands) * -math.expm1(-lad.Lphi[bj] * n / float(self.tiers[bj][2]))
                self.expsplit_b1[b1] = self.expsplit_b1.get(b1, 0.0) + ps
            for x in cands:
                x.sent += n
                x.infl += 1
                self.push(x)
            units.append(Unit(cands, bj, method, b1, n, est, pw, ps))
            nfree -= 1
        return units

    def finished(self, u):
        """u has reported and the driver has updated its candidates' tails."""
        for c in u.cands:
            c.infl -= 1
            c.done += u.n
            if c.dead:
                continue
            if c.st < len(self.stages) and c.done >= self.stages[c.st][1][-1]:
                c.st += 1
                c.sent = c.done = 0
            self.push(c)

    def cancel(self, u):
        """u was killed before it reported."""
        self.expwin -= u.pw
        self.expwin_b1[u.b1] = self.expwin_b1.get(u.b1, 0.0) - u.pw
        for c in u.cands:
            c.infl -= 1
            c.sent -= u.n
            if not c.dead:
                self.push(c)

    def backlog(self, rcap, limit=INF):
        """Core-seconds of queued sweeps with r <= rcap (stops counting at limit)."""
        tot = 0.0
        for j, h in enumerate(self.heaps):
            b1 = self.tiers[j][1]
            for (r, _, ver, c) in h:
                if r <= rcap and ver == c.ver and not c.dead:
                    tot += self.cost(b1, c.tb) * self.position(c)[1]
                    if tot >= limit:
                        return tot
        return tot

    def census(self):
        """Live queue entries per stage (for the log)."""
        out = [0] * len(self.stages)
        sweep_stage = {}
        for si, (first, ends) in enumerate(self.stages):
            for k in range(len(ends)):
                sweep_stage[first + k] = si
        for j, h in enumerate(self.heaps):
            out[sweep_stage[j]] += sum(1 for (_, _, ver, c) in h if ver == c.ver and not c.dead)
        return out


# ----------------------------------------------------------------- the level loop
def level_loop(drv, tab, tiers, cost, par=None):
    """Run one level to its first buildable win.  Returns whatever drv.build()
    returns.  The driver supplies:
      threads, now(), log(msg)
      next_chunk() -> (orders expected, core-seconds expected) of the next scan
                      chunk, or None when the scan cannot go further
      start_scan(ns), start_unit(unit)
      wait(timeout) -> list of events ("orders", payload) / ("ecm", unit, payload)
                       / ("scanthreads", k): the running scan now uses only k threads
      take_orders(payload) -> (fresh composite candidates with sched_init done,
                               number of orders bought, core-seconds they cost)
      apply(unit, payload): update tails, mark resolved candidates dead, note wins
      has_winner(), abort() (kill everything in flight), build() -> result or None
      should_build(R) (optional): False keeps the search going with a winner in
                      hand, when its build would cost more than finding another
                      (not asked once the scan is exhausted and no work is left:
                      the winner in hand is then built)
      ckpt_every, checkpoint(counters), restore() -> (candidates, counters) or
                      None (optional): the search state is handed to the driver
                      every ckpt_every seconds, and a driver that was started
                      from such a state hands it back here
    """
    import numpy as np
    import shortsim as ss
    par = dict(par or {})
    relax = par.get("relax", 2.0)
    share_lo, share_hi = par.get("share_lo", 0.15), par.get("share_hi", 0.6)
    tick = par.get("tick", 10.0)
    sch = Sched(tab, tiers, drv.threads, cost, unit_frac=par.get("unit_frac", 0.05),
                unit_min=par.get("unit_min", 1.5), quantum=int(par.get("quantum", 1)))
    H = np.zeros((tab.cmax + 1, ss.NB))
    K = 0                  # orders bought
    spent_scan = 0.0       # core-seconds
    spent_ecm = 0.0
    R = INF                # the reward at which the next chunk of orders breaks even
    Fk = Ek = None
    share = 1.0
    scan = None            # threads of the running scan
    running = {}
    uid = 0
    t0 = drv.now()
    tlog = t0
    win_b1 = None          # B1 of the unit whose report completed the first winner (0: won at intake)
    blog = (0.0, -1.0)     # (backlog, time it was computed)
    bpass = 0              # loop passes since
    nunits = 0
    over = 0.0             # core-seconds handed out above the rule (cores that would have idled)
    first = None           # (expected wins, wall) at the first win in hand

    def reprice():
        nonlocal R, share
        nxt = drv.next_chunk()
        if K == 0 or Fk is None:
            return
        if nxt is None:
            R = INF
            share = 0.0
            return
        price = nxt[1] / max(nxt[0], 1.0)
        R = tab.threshold(price, Fk)
        e = max(tab.at(R, Ek), 1e-9)
        share = min(share_hi, max(share_lo, price / (price + e)))
        sch.set_reward(R)

    ckpt_every = getattr(drv, "ckpt_every", 0.0)
    tckpt = t0
    rs = getattr(drv, "restore", None)
    rs = rs() if rs is not None else None
    if rs is not None and rs[1]["H"].shape == H.shape and rs[1]["K"] > 0:
        # a search taken up again: the orders bought, what they cost and where
        # every candidate stands; whatever was in flight is run again
        cands, ls = rs
        H, K = ls["H"], ls["K"]
        spent_scan, spent_ecm, nunits, over = ls["spent_scan"], ls["spent_ecm"], ls["nunits"], ls["over"]
        t0 = tlog = drv.now() - ls["elapsed"]
        sch.expwin, sch.nabandoned = ls["expwin"], ls["nabandoned"]
        sch.expwin_b1, sch.expsplit_b1 = dict(ls["expwin_b1"]), dict(ls["expsplit_b1"])
        for c in cands:
            sch.push(c)
        F, E = tab.fresh(H)
        Fk, Ek = F / K, E / K
        reprice()

    while True:
        now = drv.now()
        # ---- supply: start the next chunk when the work below R is running out
        nxt = drv.next_chunk()
        if scan is None and nxt is not None:
            if K == 0:
                need, ns = True, drv.threads
            else:
                ns = max(1, min(drv.threads, int(math.ceil(share * drv.threads))))
                dur = nxt[1] / ns
                low = dur * max(1, drv.threads - ns)
                bpass += 1
                if now - blog[1] > 0.25 and (bpass >= par.get("bpass", 1) or blog[1] < 0):
                    blog, bpass = (sch.backlog(R, low), now), 0
                need = blog[0] + sum(u.est for u in running.values()) * 0.5 < low
                if not running and sch.best() > R:
                    need, ns = True, drv.threads
            if need:
                drv.start_scan(ns)
                scan = ns
                blog = (0.0, -1.0)
        # ---- work: fill the cores, best first
        if K:
            nfree = drv.threads - (scan or 0) - len(running)
            if nfree > 0:
                units = sch.take(nfree, R)
                if len(units) < nfree and (scan is not None or nxt is None):
                    # cores that would idle while the scan runs take what lies just
                    # above the rule; with the scan exhausted, everything is in play
                    more = sch.take(nfree - len(units), INF if nxt is None else relax * R)
                    over += sum(u.est for u in more)
                    units += more
                for u in units:
                    uid += 1
                    u.uid, u.t0 = uid, now
                    running[uid] = u
                    drv.start_unit(u)
                nunits += len(units)
        force = False
        if not running and scan is None:
            if nxt is not None:
                drv.start_scan(drv.threads)     # nothing below the relaxed rule either: buy
                scan = drv.threads
            elif drv.has_winner():
                force = True                    # no search left: what is in reserve gets built
            else:
                raise RuntimeError("scan exhausted and nothing left to work on (%d orders)" % K)
        # ---- events
        for ev in (() if force else drv.wait(tick)):
            had = drv.has_winner()
            if ev[0] == "orders":
                cands, n, cst = drv.take_orders(ev[1])
                if not had and drv.has_winner():
                    win_b1 = 0
                scan = None
                K += n
                spent_scan += cst
                for c in cands:
                    H[min(c.cls, tab.cmax), min(c.tb, ss.NB - 1)] += 1
                    sch.push(c)
                F, E = tab.fresh(H)
                Fk, Ek = F / K, E / K
                reprice()
            elif ev[0] == "scanthreads":
                if scan is not None:        # the rest of the scan (the strip) needs fewer threads
                    scan = min(scan, ev[1])
            else:
                u = ev[1]
                drv.apply(u, ev[2])
                if not had and drv.has_winner():
                    win_b1 = u.b1
                spent_ecm += drv.now() - u.t0
                sch.finished(u)
                running.pop(u.uid, None)
        if drv.has_winner():
            if first is None:
                # the model's count at the first win is the calibration figure;
                # what is in flight has not reported and does not belong to it
                infl = sum(u.pw for u in running.values())
                first = (sch.expwin - infl, drv.now() - t0)
            hold = getattr(drv, "should_build", None)
            # hold(R) false: a winner in reserve, its build costs more than finding
            # another; the search, its log and its checkpoints go on meanwhile
            if hold is None or force or hold(R):
                drv.abort()
                for u in running.values():
                    sch.cancel(u)
                    spent_ecm += drv.now() - u.t0
                running.clear()
                scan = None
                res = drv.build({"orders": K, "scan_s": spent_scan, "ecm_s": spent_ecm,
                                 "wall": drv.now() - t0, "units": nunits, "R": R,
                                 "expwin": first[0] if hold is not None else sch.expwin, "first_win_wall": first[1],
                                 "over_s": over, "expwin_b1": dict(sch.expwin_b1), "win_b1": win_b1,
                                 "expsplit_b1": dict(sch.expsplit_b1)})
                if res is not None:
                    return res
                win_b1 = None
                first = None
        now = drv.now()
        if ckpt_every and K and now - tckpt >= ckpt_every:
            tckpt = now
            infl_b1, infl_sp = {}, {}
            for u in running.values():      # what is in flight is run again after a restore
                infl_b1[u.b1] = infl_b1.get(u.b1, 0.0) + u.pw
                infl_sp[u.b1] = infl_sp.get(u.b1, 0.0) + u.ps
            drv.checkpoint({"H": H.copy(), "K": K, "spent_scan": spent_scan,
                            "spent_ecm": spent_ecm, "nunits": nunits - len(running),
                            "over": over, "elapsed": now - t0,
                            "expwin": sch.expwin - sum(infl_b1.values()), "nabandoned": sch.nabandoned,
                            "expwin_b1": {b: w - infl_b1.get(b, 0.0) for b, w in sch.expwin_b1.items()},
                            "expsplit_b1": {b: w - infl_sp.get(b, 0.0) for b, w in sch.expsplit_b1.items()}})
        if now - tlog >= tick and K:
            tlog = now
            exp = K * (tab.at(R, Fk) + tab.at(R, Ek)) / R if R < INF else float("nan")
            drv.log("    %5.0fs: %d orders, rule R = %.0f core-s (%.2f core-h); the model expects %.2f wins from the "
                    "work done (%.2f once all orders are worked to R); %d units running%s, queued by B1 %s, "
                    "%d abandoned; scan %.0f + ecm %.0f core-s (%.0f above the rule)"
                    % (now - t0, K, R, R / 3600, sch.expwin, exp, len(running), " + scan" if scan else "",
                       sch.census(), sch.nabandoned, spent_scan, spent_ecm, over))


# ----------------------------------------------------------------- simulation
class SimCand:
    __slots__ = ("i", "W", "tail", "rem", "isp", "dead", "cls", "tb", "st", "sent", "done", "infl", "ver")


class SimDriver:
    """Virtual-time driver: orders are drawn from a real pool (tail size and
    coverage of each), their hidden medium primes from the model, ECM successes
    from the same success table the rule was built with.  Costs: a scan to B
    yields a*sqrt(B) orders and costs rate*B core-seconds (plus 8% of the range
    below, which dscan enumerates again), the strip a fixed amount per order."""

    def __init__(self, pool, lad, tab, threads, rng, cost, a=0.558, rate=2.45e-7, strip=0.01,
                 B0=6e7, g=2.0, Bmax=2e11, noise=0.15, startup=0.02, verbose=False):
        self.pool, self.lad, self.tab, self.threads, self.rng, self.cost = pool, lad, tab, threads, rng, cost
        self.a, self.rate, self.strip, self.g, self.Bmax = a, rate, strip, g, Bmax
        self.B, self.Bnext = 0.0, B0
        self.noise, self.startup, self.verbose = noise, startup, verbose
        self.t = 0.0
        self.ev = []
        self.seq = 0
        self.epoch = 0
        self.wins = 0
        self.result = None

    def now(self):
        return self.t

    def log(self, msg):
        if self.verbose:
            print(msg, flush=True)

    def _chunk(self):
        n = self.a * (math.sqrt(self.Bnext) - math.sqrt(self.B))
        cst = self.rate * ((self.Bnext - self.B) + 0.08 * self.B) + self.strip * n
        return n, cst

    def next_chunk(self):
        if self.B >= self.Bmax:
            return None
        return self._chunk()

    def _post(self, dt, ev):
        self.seq += 1
        heapq.heappush(self.ev, (self.t + dt, self.seq, self.epoch, ev))

    def start_scan(self, ns):
        n, cst = self._chunk()
        self._post(cst / ns, ("orders", (int(round(n)), cst)))

    def start_unit(self, u):
        dur = self.startup + sum(u.n * float(self.cost(u.b1, c.tb)) for c in u.cands)
        if self.noise:
            dur *= math.exp(self.noise * self.rng.standard_normal())
        self._post(dur, ("ecm", u, None))

    def wait(self, timeout):
        while self.ev:
            t, _, ep, ev = heapq.heappop(self.ev)
            if ep != self.epoch:
                continue
            self.t = t
            return [ev]
        return []

    def take_orders(self, payload):
        n, cst = payload
        pool, lam = self.pool, self.pool.lam
        self.B, self.Bnext = self.Bnext, min(self.Bnext * self.g, self.Bmax)
        idx = self.rng.integers(0, len(pool.t), size=n)
        out = []
        for i in idx:
            t = float(pool.t[i])
            W = float(pool.w[i])
            if t <= 1 or t <= lam - W:
                continue
            if pool.prp[i]:
                if 0 < lam - t < W:
                    self.wins += 1
                continue
            c = SimCand()
            c.i, c.W, c.tail, c.rem, c.isp, c.dead = int(i), W, t, None, False, False
            sched_init(c, int(round(W)), int(round(t)))
            out.append(c)
        return out, n, cst

    def apply(self, u, payload):
        import shortsim as ss
        pi = self.lad.pi[u.j]
        lam = self.pool.lam
        for c in u.cands:
            if c.dead:
                continue
            if c.rem is None:
                pr, isp = ss.sample_tails(self.lad, self.pool, self.rng, [c.i])
                c.rem, c.isp = list(pr[0]), bool(isp[0])
            if not c.rem:
                continue
            keep = []
            for b in c.rem:
                p = float(pi[min(int(round(b)), ss.NB - 1)])
                if self.rng.random() < 1.0 - (1.0 - p) ** u.n:
                    c.tail -= b
                else:
                    keep.append(b)
            if len(keep) == len(c.rem):
                continue
            c.rem = keep
            if not keep and c.isp:
                c.dead = True
                if 0 < lam - c.tail < c.W:
                    self.wins += 1
            elif c.tail <= lam - c.W:
                c.dead = True
            c.tb = int(round(c.tail))

    def has_winner(self):
        return self.wins > 0

    def abort(self):
        self.epoch += 1
        self.ev = []

    def build(self, stats):
        return stats


def sim_main(argv):
    import pickle
    import numpy as np
    import shortsim as ss
    kv = dict(a.split("=", 1) for a in argv[1:] if "=" in a)
    pool = pickle.load(open(argv[0], "rb"))
    threads = int(kv.get("threads", 192))
    runs = int(kv.get("runs", 200))
    seed = int(kv.get("seed", 1))
    bits = int(kv.get("bits", 769))
    scale = float(kv.get("cost", 1.0))
    ss.set_grid(bits)
    here = os.path.dirname(os.path.abspath(__file__))
    engine = kv.get("engine", "gmp")
    em, sec, cost, quantum = ss.engine_model(engine, scale, here)
    tiers = default_tiers(engine)
    if "tiers" in kv:                         # tiers=B1xCURVESxSWEEPS,...
        tiers = []
        for t in kv["tiers"].split(","):
            b1, c, k = t.split("x")
            tiers += [("ecm", int(float(b1)), int(c))] * int(k)
    lad = ss.Ladder(tiers, em, y0bits=26, cost=cost)
    t0 = time.time()
    tab = ss.IndexTable(lad, pool.lam, threads=8)
    print("table built in %.1fs (%d rewards, %d sweeps)" % (time.time() - t0, len(tab.rgrid), len(tiers)))
    par = {k: float(kv[k]) for k in ("relax", "unit_frac", "unit_min", "share_lo", "share_hi", "bpass") if k in kv}
    par["quantum"] = quantum
    dk = {k: float(kv[k]) for k in ("a", "rate", "strip", "B0", "g", "Bmax", "noise") if k in kv}
    rows = []
    t0 = time.time()
    won, exp_b1 = {}, {}
    for k in range(runs):
        rng = np.random.default_rng(seed * 100003 + k)
        drv = SimDriver(pool, lad, tab, threads, rng, cost, verbose=(kv.get("v") == "1" and k == 0), **dk)
        st = level_loop(drv, tab, tiers, cost, par)
        rows.append((st["wall"], st["orders"], st["scan_s"], st["ecm_s"], st["units"], st["R"], st["expwin"],
                     st["over_s"]))
        won[st["win_b1"]] = won.get(st["win_b1"], 0) + 1
        for b1, v in st["expwin_b1"].items():
            exp_b1[b1] = exp_b1.get(b1, 0.0) + v
    a = np.array(rows)
    busy = a[:, 2] + a[:, 3]
    print("%d levels simulated in %.0fs, %d threads" % (runs, time.time() - t0, threads))
    print("  wall to the first win     mean %8.1f s   median %8.1f s" % (a[:, 0].mean(), np.median(a[:, 0])))
    print("  machine cost (wall x thr) mean %8.0f core-s = %.2f core-h" % (a[:, 0].mean() * threads, a[:, 0].mean() * threads / 3600))
    print("  work done                 mean %8.0f core-s = %.2f core-h  (scan %.0f, ecm %.0f); utilisation %.0f%%"
          % (busy.mean(), busy.mean() / 3600, a[:, 2].mean(), a[:, 3].mean(), 100 * busy.sum() / (a[:, 0].sum() * threads)))
    print("  orders bought             mean %8.0f   median %8.0f;  units %.0f;  rule R at the end: median %.0f core-s"
          % (a[:, 1].mean(), np.median(a[:, 1]), a[:, 4].mean(), np.median(a[:, 5])))
    print("  wins the model expected from the work done, at the first win: mean %.2f (1 if the rule's bookkeeping is right)"
          % a[:, 6].mean())
    print("  ECM handed out above the rule to cores that would have idled: %.0f%% of the ECM work"
          % (100 * a[:, 7].sum() / a[:, 3].sum()))
    print("  first wins by the B1 of the unit that completed them (0 = at intake), and the wins the model expected there:")
    print("    " + "  ".join("%g: %d / %.1f" % (b1, won.get(b1, 0), exp_b1.get(b1, 0.0))
                             for b1 in sorted(set(won) | set(exp_b1), key=lambda x: x or 0)))
    return a


def default_tiers(engine="gmp", nbits=None):
    """The ladder of sweeps the rule chooses from: (method, B1, curves per sweep).
    yafu runs 8 curves at a time, so its sweeps are multiples of 8.  For moduli
    of nbits bits the ladder stops at the B1 whose curves still find next to
    nothing at shortsim.medium_bits(nbits) bits (the model takes larger primes
    as out of reach): 11000 below 230 bits, 50000 below 320, 250000 below 410,
    10^6 below 500, 3*10^7 from there.  From 1300 bits a win is worth a thousand
    core-hours and more, and the best candidates run through 64 curves at 3*10^7
    with the rule still willing: yafu's ladder then offers 192."""
    if engine == "yafu":
        spec = [(2000, 16, 1), (11000, 24, 1), (50000, 16, 2), (250000, 16, 3),
                (1000000, 16, 5), (3000000, 16, 6), (10000000, 16, 6),
                (30000000, 16, 12 if (nbits or 0) >= 1300 else 4)]
    else:
        spec = [(2000, 12, 1), (11000, 20, 1), (50000, 16, 2), (250000, 12, 4),
                (1000000, 12, 6), (3000000, 12, 8), (10000000, 12, 8), (30000000, 12, 4)]
    top = 30000000
    if nbits is not None:
        top = 11000 if nbits < 230 else 50000 if nbits < 320 else 250000 if nbits < 410 else \
            1000000 if nbits < 500 else 30000000
    spec = [x for x in spec if x[0] <= top]
    tiers = []
    for (b1, c, k) in spec:
        tiers += [("ecm", b1, c)] * k
    return tiers


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "sim":
        sim_main(sys.argv[2:])
    else:
        print(__doc__)
