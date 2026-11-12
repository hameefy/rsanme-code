"""Redraw the Dolan-More performance profiles from results/results_monotone.csv.

Same definition as perfprof.m (failures = infinite cost, ratio to the best solver, log axis), but one
full-width figure per metric, distinct markers and line styles, and a zoomed inset near ratio 1.
Usage:  python make_profiles.py            (run from the drivers folder)
Writes  profile_iter / profile_fev / profile_time  as .pdf, .png, .eps in ../results
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
res = os.path.join(here, "..", "results")
d = pd.read_csv(os.path.join(res, "results_monotone.csv"))
names = ["RSANME", "RSANME-Q", "Modified Alg. 2.1", "Liu Alg. 2.1", "DFSANME", "PAND-SR"]
style = {
    "RSANME":            dict(color="#0033FF", ls="-",  marker="o"),
    "RSANME-Q":          dict(color="#FF1A1A", ls="--", marker="s"),
    "Modified Alg. 2.1": dict(color="#FF9900", ls="-.", marker="^"),
    "Liu Alg. 2.1":      dict(color="#A800FF", ls=":",  marker="D"),
    "DFSANME":           dict(color="#00B21E", ls="-",  marker="v"),
    "PAND-SR":           dict(color="#FF00B4", ls="--", marker="P"),
}
metrics = [("Iter", "iterations", "profile_iter"),
           ("Fev", "function evaluations", "profile_fev"),
           ("Time", "computing time", "profile_time")]
key = ["Problem", "Dim", "Start"]
plt.rcParams.update({"font.size": 12, "axes.labelsize": 13, "legend.fontsize": 11, "font.family": "serif",
                     "font.serif": ["STIXGeneral", "DejaVu Serif"], "mathtext.fontset": "stix",
                     "ps.fonttype": 42, "pdf.fonttype": 42, "axes.linewidth": 1.1})

for m, label, out in metrics:
    P = d.pivot_table(index=key, columns="Solver", values=m)
    S = d.pivot_table(index=key, columns="Solver", values="Success")
    C = P.where(S > 0, np.inf)[names]
    best = C.min(axis=1)
    C = C[np.isfinite(best)]
    best = best[np.isfinite(best)]
    tmax = 1e3
    taus = np.unique(np.concatenate([[1.0], np.logspace(0, np.log10(tmax), 600)]))
    rho = {s: np.array([(C[s].values <= t * best.values * (1 + 1e-12)).mean() for t in taus]) for s in names}

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for k, s in enumerate(names):
        st = style[s]
        ax.step(taus, rho[s], where="post", color=st["color"], ls=st["ls"], lw=2.4, label=s)
        # markers at log-spaced ratios, shifted per solver so overlapping curves stay distinguishable
        mt = np.logspace(0, np.log10(tmax), 11) * (1 + 0.035 * k)
        mt = mt[mt <= tmax]
        mv = [rho[s][min(np.searchsorted(taus, t), len(taus) - 1)] for t in mt]
        ax.plot(mt, mv, ls="none", marker=st["marker"], ms=7.5, mfc="white" if k % 2 else st["color"],
                mec=st["color"], mew=1.8)
    ax.set_xscale("log")
    ax.set_xlim(1, tmax)
    ax.set_ylim(0, 1.03)
    ax.set_xlabel("performance ratio w.r.t. " + label)
    ax.set_ylabel(r"$P(\tau)$")
    ax.grid(True, which="major", color="#C8C8C8", lw=0.8)
    ax.grid(True, which="minor", color="#E6E6E6", lw=0.5)
    handles = [plt.Line2D([], [], color=style[s]["color"], ls=style[s]["ls"], lw=2.4, marker=style[s]["marker"],
                          ms=7.5, mfc="white" if k % 2 else style[s]["color"], mec=style[s]["color"], mew=1.8)
               for k, s in enumerate(names)]
    ax.legend(handles, names, loc="center right", bbox_to_anchor=(1.0, 0.42), framealpha=1.0, edgecolor="#555555", ncol=1)

    # zoomed inset: ratios between 1 and 2 on a linear axis
    ins = ax.inset_axes([0.34, 0.08, 0.30, 0.34])
    z = np.linspace(1, 2, 400)
    for k, s in enumerate(names):
        st = style[s]
        y = np.array([(C[s].values <= t * best.values * (1 + 1e-12)).mean() for t in z])
        ins.step(z, y, where="post", color=st["color"], ls=st["ls"], lw=1.9)
    top4 = [rho[s][0] for s in ("RSANME", "RSANME-Q", "DFSANME", "PAND-SR")]
    lo = max(0.0, np.floor(min(top4) * 20) / 20 - 0.05)
    ins.set_xlim(1, 2)
    ins.set_ylim(lo, 1.01)
    ins.set_title("zoom: ratio 1 to 2", fontsize=10, pad=2)
    ins.tick_params(labelsize=9)
    ins.grid(True, color="#C8C8C8", lw=0.6)
    fig.tight_layout()
    for ext, kw in (("pdf", {}), ("png", {"dpi": 300}), ("eps", {})):
        fig.savefig(os.path.join(res, f"{out}.{ext}"), **kw)
    plt.close(fig)
    print(out, "instances:", len(C), {s: round(float(rho[s][-1]), 3) for s in names})
