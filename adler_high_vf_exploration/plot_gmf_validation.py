"""Plot closed GMF formulas against periodic numerical estimates."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from gaussian_minkowski import margin_gmfs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "csv", nargs="?", type=Path, default=Path("results/numerical_gmfs/numerical_gmfs.csv")
    )
    parser.add_argument("--output", type=Path, default=Path("figures/general_gmfs_validation.png"))
    args = parser.parse_args()

    rows = []
    with args.csv.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            rows.append({key: float(value) for key, value in row.items()})
    q_values = sorted({int(row["q"]) for row in rows})
    fig, axes = plt.subplots(1, len(q_values), figsize=(3.5 * len(q_values), 3.8), constrained_layout=True)
    axes = np.atleast_1d(axes)
    for axis, q in zip(axes, q_values):
        subset = sorted((row for row in rows if int(row["q"]) == q), key=lambda row: row["dimension"])
        delta = subset[0]["delta"]
        numerical = np.array([0.7] + [row["estimate"] for row in subset])
        sem = np.array([0.0] + [row["sem"] for row in subset])
        analytic = margin_gmfs(q, delta)
        x = np.arange(4)
        axis.errorbar(x - 0.06, numerical, yerr=sem, fmt="o", capsize=3, label="periodic numerical")
        axis.plot(x + 0.06, analytic, "s", label="closed formula")
        axis.axhline(0.0, color="black", lw=0.7)
        axis.set_xticks(x, [rf"$\mathcal{{M}}_{j}$" for j in x])
        axis.set_title(rf"$q={q}$")
        axis.grid(alpha=0.3)
        axis.legend(fontsize=7)
    fig.suptitle(r"General Gaussian Minkowski functionals at $\phi=0.70$")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
