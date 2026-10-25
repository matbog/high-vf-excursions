"""Aggregate margin-tessellation sweeps and create a compact figure/table."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("--figure", type=Path, default=Path("results/sweep_summary.png"))
    parser.add_argument("--markdown", type=Path, default=Path("results/sweep_summary.md"))
    args = parser.parse_args()

    grouped: dict[tuple[int, float], list[dict[str, float | bool]]] = defaultdict(list)
    with args.csv.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            parsed = {
                "phi": float(row["measured_phi"]),
                "largest": float(row["largest_component_fraction_phase"]),
                "spans": any(row[key] == "True" for key in ("spans_x", "spans_y", "spans_z")),
            }
            grouped[(int(row["q"]), float(row["theoretical_phi"]))].append(parsed)

    summary = []
    for (q, target), rows in sorted(grouped.items()):
        summary.append(
            {
                "q": q,
                "target": target,
                "mean_phi": float(np.mean([r["phi"] for r in rows])),
                "mean_largest": float(np.mean([r["largest"] for r in rows])),
                "span_probability": float(np.mean([r["spans"] for r in rows])),
                "n": len(rows),
            }
        )

    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "| q | target phi | measured phi | largest / phase | span probability | n |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary:
        lines.append(
            f"| {row['q']} | {row['target']:.2f} | {row['mean_phi']:.3f} | "
            f"{row['mean_largest']:.3f} | {row['span_probability']:.2f} | {row['n']} |"
        )
    args.markdown.write_text("\n".join(lines) + "\n", encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for q in sorted({row["q"] for row in summary}):
        subset = [row for row in summary if row["q"] == q]
        x = [row["target"] for row in subset]
        axes[0].plot(x, [row["span_probability"] for row in subset], "o-", label=f"q={q}")
        axes[1].plot(x, [row["mean_largest"] for row in subset], "o-", label=f"q={q}")
    axes[0].set(xlabel="target volume fraction", ylabel="spanning probability", ylim=(-0.05, 1.05))
    axes[1].set(xlabel="target volume fraction", ylabel="largest component / retained phase")
    for axis in axes:
        axis.grid(alpha=0.3)
        axis.legend()
    fig.suptitle("Gaussian margin tessellation pilot")
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.figure, dpi=200)


if __name__ == "__main__":
    main()

