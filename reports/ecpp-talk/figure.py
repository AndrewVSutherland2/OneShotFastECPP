#!/usr/bin/env python3
"""figure.py -- "Are short proofs hard to find?": find-vs-verify cost of primality
certificates, after AVS's original PNG, with staged builds for the talk.

    python3 figure.py            # writes figs/gap-1.pdf .. figs/gap-4.pdf and figs/gap.pdf

Stages: 1 = AKS / randomized-AKS / Goldwasser-Kilian ECPP; 2 = + Pomerance;
3 = + one-shot; 4 = + short ECPP (the full figure, also written as gap.pdf, with
title; the staged builds omit the title so a frame title can carry it).

Data changed vs the original PNG (Oct 2026): largest Pomerance triple 87 -> 90 bits
(10^27+103, DANGER3, 2026-08-28); largest short ECPP 1,030 -> 1,329 bits (10^400+69);
deterministic AKS plotted at exponent 6 (Lenstra-Pomerance) rather than 6.5.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BLUE, RED, GREEN, GREY = "#2068b0", "#d0352b", "#2e8b2e", "#777777"
LABELS = ["AKS\n(Lenstra–\nPomerance)", "randomized\nAKS\n(Bernstein)",
          "Goldwasser–\nKilian ECPP\n(Atkin–Morain)", "Pomerance\nECPP",
          "one-shot\nECPP", "short\nECPP"]
BOLD = {4, 5}
BITS = [None, "1,025 bits", "363,074 bits", "90 bits", "399 bits", "1,329 bits"]

# y-coordinates: polynomial region is the true exponent c of (log p)^c; the band
# above Y_BREAK is schematic ("not to scale").
Y_POLY_TOP, Y_BREAK = 7.45, 7.95
Y_SUB_TOP, Y_TOP = 10.6, 14.0
FIND = [6.0, 4.0, 4.0, 13.45, 12.7, 9.5]       # p^{1/2}, p^{1/8}, L(1/3) are schematic
VERIFY = [None, 4.0, 3.0, 2.0, 2.0, 2.0]


def gap_label(ax, text, at, top, bottom, weight):
    """Centered label at `at` with arrows leaving its top edge to `top` and its bottom
    edge to `bottom` (the second copy of the text is invisible, just an arrow anchor)."""
    for target, rel, color in ((top, (0.5, 1.0), "k"), (bottom, (0.5, 0.0), "none")):
        ax.annotate(text, xy=target, xytext=at, ha="center", va="center", fontsize=12.5,
                    fontweight=weight, color=color,
                    arrowprops=dict(arrowstyle="-|>", color="k", lw=1.5, mutation_scale=14,
                                    relpos=rel, shrinkA=2, shrinkB=2))


def draw(stage, title):
    npts = {1: 3, 2: 4, 3: 5, 4: 6}[stage]
    fig, ax = plt.subplots(figsize=(8.4, 4.6 if title else 4.0))
    ax.set_xlim(-0.55, 6.15)
    ax.set_ylim(0.15, Y_TOP)

    # bands
    ax.axhspan(Y_SUB_TOP, Y_TOP - 0.05, xmin=0.0, xmax=1.0, color="#ececec", zorder=0)
    ax.axhspan(Y_BREAK, Y_SUB_TOP, color="#f6f6f6", zorder=0)
    ax.axhline(Y_SUB_TOP, color="#888888", ls="--", lw=1.3, zorder=1)
    for y in range(2, 8):
        ax.axhline(y, color="#e6e6e6", lw=0.8, zorder=0)
    ax.text(-0.45, Y_TOP - 0.95, "exponential", style="italic", fontsize=11, color="#444444")
    ax.text(-0.45, Y_TOP - 1.75, "(not to scale)", fontsize=8.5, color=GREY)
    ax.text(-0.45, Y_SUB_TOP - 0.95, "subexponential", style="italic", fontsize=11, color="#444444")
    ax.text(-0.45, Y_SUB_TOP - 1.75, "(not to scale)", fontsize=8.5, color=GREY)
    ax.text(-0.45, 6.95, "polynomial", style="italic", fontsize=11, color="#444444")

    # axes: left spine only over the polynomial range, with break marks
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_bounds(0.15, Y_POLY_TOP)
    ax.set_yticks(range(2, 8))
    ax.tick_params(axis="y", labelsize=10.5)
    ax.set_ylabel(r"cost exponent  $c$  in  $(\log p)^c$", fontsize=11.5)
    ax.yaxis.set_label_coords(-0.05, 0.33)
    kw = dict(transform=ax.transData, color="k", clip_on=False, lw=1.2)
    for dy in (0.0, 0.25):
        ax.plot([-0.62, -0.48], [Y_POLY_TOP + 0.08 + dy, Y_POLY_TOP + 0.30 + dy], **kw)
    ax.set_xticks(range(6))
    ax.set_xticklabels(LABELS, fontsize=10)
    for i, t in enumerate(ax.get_xticklabels()):
        if i in BOLD:
            t.set_fontweight("bold")
    ax.tick_params(axis="x", length=0, pad=8)

    # largest p found (by search)
    ax.text(-0.45, 0.62, "largest $p$ found\n(by search):", fontsize=8.5, style="italic",
            color="#444444", linespacing=0.95)
    for i in range(1, npts):
        ax.text(i, 0.75, BITS[i], ha="center", fontsize=9.5, color="#444444")

    # FIND curve (red): straight segments, with an exponential rise GK -> Pomerance
    xs, ys = [0.0, 1.0, 2.0], [FIND[0], FIND[1], FIND[2]]
    if npts >= 4:
        t = np.linspace(0, 1, 80)
        a = 4.2
        xs += list(2 + t)
        ys += list(FIND[2] + (FIND[3] - FIND[2]) * (np.exp(a * t) - 1) / (np.exp(a) - 1))
    for i in range(4, npts):
        xs.append(float(i)); ys.append(FIND[i])
    ax.plot(xs, ys, color=RED, lw=3.4, zorder=3, solid_capstyle="round")
    for i in range(npts):
        ax.plot(i, FIND[i], "o", color=RED, ms=17 if i == 5 else 12, zorder=4)

    # VERIFY curve (blue)
    vx = [i for i in range(1, npts)]
    ax.plot(vx, [VERIFY[i] for i in vx], color=BLUE, lw=3.4, zorder=3)
    for i in vx:
        ax.plot(i, VERIFY[i], "o", color=BLUE, ms=12, zorder=4)
    # the small find/verify gap of the polynomial methods
    ax.fill_between([1, 2], [4, 4], [4, 3], color=RED, alpha=0.08, lw=0, zorder=2)

    arrow = dict(arrowstyle="-|>", color="k", lw=1.5, mutation_scale=14)
    ax.annotate("for AKS the proof IS\nthe computation\n(verify = find)", xy=(0.5, 5.0),
                xytext=(1.0, 1.55), ha="center", fontsize=9, color=GREY,
                arrowprops=dict(arrowstyle="-|>", color=GREY, lw=1.3, mutation_scale=12))

    if npts >= 4:
        ax.text(2.72, FIND[3] - 0.5, r"$\sim\!p^{1/2}$", color=RED, fontsize=13.5, ha="right")
        gap_label(ax, "exponential\ngap", (3.62, 5.6), (3.06, FIND[3] - 0.6), (3.04, 2.25), "normal")
    if npts >= 5:
        ax.text(3.72, FIND[4] - 1.75, r"$\sim\!p^{1/8}$", color=RED, fontsize=13.5, ha="center")
    if npts >= 6:
        ax.text(4.78, FIND[5] - 0.4, r"$L(1/3)$", color=RED, fontsize=14.5, ha="right")
        ax.text(4.15, 2.75, "quasi-quadratic:\n" r"$(\log p)^{2+o(1)}$", color=BLUE, fontsize=10,
                ha="center")
        # where the red segment one-shot -> short crosses the barrier
        xc = 4 + (FIND[4] - Y_SUB_TOP) / (FIND[4] - FIND[5])
        ax.annotate("breaks the\nexponential barrier!", xy=(xc + 0.03, Y_SUB_TOP + 0.12),
                    xytext=(5.3, Y_TOP - 0.2), ha="center", va="top", fontsize=11,
                    fontweight="bold", color=GREEN,
                    arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=2, mutation_scale=16))
        gap_label(ax, "subexponential\ngap", (5.5, 5.3), (5.1, FIND[5] - 0.5), (5.07, 2.2), "bold")

    handles = [Line2D([], [], color=BLUE, lw=3.4, marker="o", ms=10, label="cost to VERIFY proof"),
               Line2D([], [], color=RED, lw=3.4, marker="o", ms=10,
                      label="cost to FIND proof (heuristic)")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.31), ncol=2,
              frameon=False, fontsize=11, handlelength=2.2, columnspacing=3)
    if title:
        fig.suptitle("Are short proofs hard to find?", fontsize=17, fontweight="bold", y=0.985)
        fig.text(0.5, 0.905, "certifying that an integer $p$ is prime", ha="center",
                 fontsize=11, style="italic", color="#444444")
        fig.subplots_adjust(left=0.085, right=0.985, top=0.87, bottom=0.255)
    else:
        fig.subplots_adjust(left=0.085, right=0.985, top=0.99, bottom=0.29)
    return fig


if __name__ == "__main__":
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs")
    os.makedirs(out, exist_ok=True)
    for st in (1, 2, 3, 4):
        draw(st, title=False).savefig(os.path.join(out, f"gap-{st}.pdf"))
    f = draw(4, title=True)
    f.savefig(os.path.join(out, "gap.pdf"))
    f.savefig(os.path.join(out, "gap.png"), dpi=200)
