"""Plot measured vs theory L0..L3 as a function of delta for one (q, N) pair.

Reads the existing margin_sweep.csv and produces a 4-panel figure
(one panel per LKC) with theory as a solid line and numerical mean ± std
as points with error bars.

Usage:
    mamba run -n rf_explo python plot_lkc_vs_delta.py --q 16 --N 256
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RESULTS = Path(__file__).parent / "results"


def load(q: int, N: int) -> list[dict]:
    rows = []
    with (RESULTS / "margin_sweep.csv").open() as f:
        for row in csv.DictReader(f):
            if int(row["q"]) == q and int(row["N"]) == N:
                def _parse(k, v):
                    if k == "realization": return int(v)
                    if k == "spans_any": return v.lower() == "true"
                    return float(v)
                rows.append({k: _parse(k, v) for k, v in row.items()})
    return rows


def aggregate(rows: list[dict]) -> dict:
    """Group by delta, return mean ± std for measured quantities."""
    by_phi: dict[float, list] = {}
    for r in rows:
        by_phi.setdefault(r["target_phi"], []).append(r)

    out = {"delta": [], "theory_phi": [], "theory_euler": [], "theory_l1": [], "theory_l2": [],
           "phi_m": [], "phi_s": [], "euler_m": [], "euler_s": [],
           "l1_m": [], "l1_s": [], "l2_m": [], "l2_s": []}

    for phi in sorted(by_phi):
        cell = by_phi[phi]
        out["delta"].append(np.mean([r["delta"] for r in cell]))
        out["theory_phi"].append(cell[0]["theory_phi"])
        out["theory_euler"].append(cell[0]["theory_euler"])
        out["theory_l1"].append(cell[0]["theory_l1"])
        out["theory_l2"].append(cell[0]["theory_l2"])
        for key, col in [("phi", "measured_phi"), ("euler", "euler"),
                         ("l1", "l1_winner"), ("l2", "l2_winner")]:
            vals = [r[col] for r in cell]
            out[f"{key}_m"].append(np.mean(vals))
            out[f"{key}_s"].append(np.std(vals, ddof=1) if len(vals) > 1 else 0.0)

    return {k: np.array(v) for k, v in out.items()}


def plot(q: int, N: int) -> None:
    rows = load(q, N)
    if not rows:
        raise ValueError(f"No data for q={q}, N={N} in margin_sweep.csv")
    d = aggregate(rows)

    fig, axes = plt.subplots(1, 4, figsize=(16, 4), constrained_layout=True)

    panels = [
        ("L0  (Euler characteristic)", d["theory_euler"], d["euler_m"], d["euler_s"], None),
        ("L1  (mean breadth)",    d["theory_l1"],  d["l1_m"],    d["l1_s"],    None),
        ("L2  (surface area)",    d["theory_l2"],  d["l2_m"],    d["l2_s"],    None),
        ("L3  (volume fraction)", d["theory_phi"], d["phi_m"],   d["phi_s"],   None),
    ]

    for ax, (title, theory, meas, std, _) in zip(axes, panels):
        ax.plot(d["delta"], theory, "-", color="steelblue", lw=2, label="GKF theory")
        ax.errorbar(d["delta"], meas, yerr=std, fmt="o", color="tomato",
                    ms=5, capsize=3, lw=1.4, label=f"numerical  (N={N})")
        ax.set_xlabel(r"$\delta$  (margin threshold)", fontsize=11)
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)

    phi0 = next(iter(set(r["target_phi"] for r in rows)))
    n_real = sum(1 for r in rows if r["target_phi"] == phi0)
    fig.suptitle(f"Margin tessellation GKF validation — $q={q}$, $N={N}$, "
                 rf"$L=64$, $l_c=5$  ({n_real} realisations)",
                 fontsize=13)

    out = RESULTS / f"lkc_vs_delta_q{q}_n{N}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--q", type=int, default=16)
    parser.add_argument("--N", type=int, default=256)
    args = parser.parse_args()
    plot(args.q, args.N)
