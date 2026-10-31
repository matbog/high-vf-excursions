"""Plot identical physical slices of the PL volumes and their centre voxels."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap
import numpy as np
import pyvista as pv

from gstools_fields import center_samples
from model import margin_from_fields


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--realization", type=int, default=1)
    args = parser.parse_args()
    out = args.results.resolve()
    here = Path(__file__).resolve().parent
    if not out.is_relative_to(here) or out == here:
        raise ValueError("Select a result subdirectory inside experiment 07.")
    config = json.loads((out/"metadata.json").read_text())["config"]
    n, length = config["n"], config["length"]
    slice_index = n//2
    z = (slice_index+.5)*length/n
    with np.load(out/f"fields_r{args.realization}.npz") as saved:
        groups = [[saved[f"g_{k}_{j}"] for j in range(item["q"])] for k, item in enumerate(config["classes"])]
    margins = [margin_from_fields(fields)[0] for fields in center_samples(groups)]
    with (out/"theory.csv").open(newline="") as stream:
        theory = list(csv.DictReader(stream))
    fig, axes = plt.subplots(len(config["target_phis"]), 2, figsize=(9, 4*len(config["target_phis"])),
                             squeeze=False, constrained_layout=True)
    cmap = ListedColormap(["#f3f1ea", "#55758c"])
    for row, phi in enumerate(config["target_phis"]):
        deltas = [float(next(t["delta"] for t in theory if t["kind"] == "raw"
                            and int(t["class_index"]) == k and float(t["target_phi"]) == phi))
                  for k in range(len(groups))]
        mask = np.logical_or.reduce([m >= d for m, d in zip(margins, deltas)])
        axes[row, 0].imshow(mask[:, :, slice_index].T, origin="lower", cmap=cmap,
                            vmin=0, vmax=1, interpolation="nearest", extent=(0, length, 0, length))
        mesh_dir = out/f"meshes_r{args.realization}_vf_{phi:.4f}"
        for path in sorted(mesh_dir.glob("class_*_winner_*.vtu")):
            sliced = pv.read(path).slice(normal="z", origin=(length/2, length/2, z)).triangulate()
            if sliced.n_cells:
                faces = sliced.faces.reshape(-1, 4)[:, 1:]
                polygons = sliced.points[faces, :2]
                axes[row, 1].add_collection(PolyCollection(polygons, facecolors="#55758c", edgecolors="none",
                                                          antialiaseds=False))
        for column, label in enumerate(("Centre voxels of the same PL fields", "Clipped component interpolants")):
            ax = axes[row, column]
            ax.set(xlim=(0, length), ylim=(0, length), aspect="equal", facecolor="#f3f1ea",
                   title=f"Target VF {phi:.0%}: {label}", xlabel=config["length_unit"], ylabel=config["length_unit"])
    fig.suptitle(f"Experiment 07: N={n}, q={[c['q'] for c in config['classes']]}, z={z:.3f} {config['length_unit']}\n"
                 "Same fields and slice plane; all aggregate classes shown in blue")
    path = out/"geometry_slices.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(path)


if __name__ == "__main__":
    main()
