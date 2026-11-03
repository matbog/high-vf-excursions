"""High-resolution ceiling-push scene at N=512 to test physical disconnection.

Same (q, phi) as scene_ceiling.py but at N=512 instead of 128. The point:
at N=128 the film is 0.27 voxels (sub-voxel) and the physical binary
union merges most winner cells into one giant cluster despite the
per-winner construction being disconnected. At N=512 the film is 1.09
voxels (just resolved) — this is where we expect the physical binary
union to actually disconnect into separate particles.

Measures BOTH per-winner and physical binary-union spanning + component
counts, then exports a full ParaView scene set for direct comparison
with the N=128 scene.

Expected runtime: ~3 h on a laptop (dominated by field generation).
Peak RAM: ~6 GB.

Output in paraview/q32_phi75_n512/:
    field_g1.vti    first Gaussian competitor
    margin.vti      max - 2nd_max scalar
    phase.vti       binary inclusion mask
    winners.vti     winning field index 0..q-1
    components.vti  per-winner component id (unique int per aggregate)
    components_physical.vti  binary-union connected-component id
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
from skimage import measure

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[1]))

from src import adler, excursions, fields, io, measures  # noqa: E402
from src.fields import GridSpec  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("--N",   type=int,   default=512)
parser.add_argument("--q",   type=int,   default=32)
parser.add_argument("--phi", type=float, default=0.75)
args = parser.parse_args()

N      = args.N
LENGTH = 64.0
LC     = 5.0
Q      = args.q
PHI    = args.phi
SEED   = 1729 + 10000 * Q

OUT = HERE / "paraview" / f"q{Q}_phi{int(PHI*100):02d}_n{N}"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> None:
    grid = GridSpec(n=N, length=LENGTH)
    rng = np.random.default_rng(SEED)

    print(f"=== High-res scene: q={Q}, phi={PHI}, N={N} ===")
    print(f"grid: L={LENGTH}, dx={grid.dx:.4f}, lc={LC}, L/lc={LENGTH/LC:.1f}")

    # Rough per-field time scaling: 5.4s at N=128, 40s at N=256, 320s at N=512
    # (empirical, benchmarks from Phase B)
    per_field_s = 5.4 * (N / 128) ** 3
    t0 = time.time()
    print(f"\n[1/5] Generating {Q} Gaussian fields at N={N} "
          f"(~{Q * per_field_s / 60:.1f} min expected)...")
    fs = fields.generate_margin_fields(grid, LC, Q, rng)
    print(f"       done in {(time.time() - t0) / 60:.1f} min")
    print(f"       shape = {fs.shape}, mem = {fs.nbytes / 1e9:.1f} GB")

    delta = adler.delta_for_vf(Q, PHI)
    film_voxels = delta / grid.dx
    print(f"\n[2/5] delta for phi={PHI}: {delta:.5f} physical  "
          f"({film_voxels:.3f} voxels)")
    if film_voxels < 1:
        print(f"       WARNING: film still sub-voxel at N={N}")
    else:
        print(f"       film is resolved (>1 voxel)")

    t0 = time.time()
    print(f"\n[3/5] Computing margin + phase + winners...")
    top2   = np.partition(fs, kth=Q - 2, axis=0)[-2:]
    margin = (top2[1] - top2[0]).astype(np.float32)
    del top2
    phase, winners = excursions.margin_excursion(fs, delta)
    measured_vf = float(phase.mean())
    print(f"       measured VF: {measured_vf:.4f}  (target {PHI})")
    print(f"       done in {time.time() - t0:.1f}s")

    # Free field stack memory before heavy component labelling
    del fs

    # Per-winner components (same as scene_ceiling.py)
    t0 = time.time()
    print(f"\n[4/5] Labelling per-winner components...")
    components = np.zeros_like(winners, dtype=np.int32)
    cur = 1
    for w in range(Q):
        mask = phase & (winners == w)
        if not mask.any():
            continue
        lab = measure.label(mask, connectivity=1)
        n_comp = int(lab.max())
        components[lab > 0] = lab[lab > 0] + (cur - 1)
        cur += n_comp
    n_per_winner = cur - 1
    print(f"       per-winner components: {n_per_winner}")

    # Physical binary-union components (ignoring winner labels)
    components_phys = measure.label(phase, connectivity=1).astype(np.int32)
    n_physical = int(components_phys.max())
    print(f"       physical (binary-union) components: {n_physical}")

    # Largest physical component
    sizes = np.bincount(components_phys.ravel())
    sizes[0] = 0   # ignore matrix
    largest_id   = int(sizes.argmax())
    largest_frac = float(sizes[largest_id] / phase.sum())
    print(f"       largest physical component: id={largest_id}, "
          f"holds {largest_frac * 100:.1f}% of inclusion volume")
    print(f"       done in {time.time() - t0:.1f}s")

    # Spanning checks (both kinds)
    t0 = time.time()
    print(f"\n[5/5] Spanning checks...")
    perc_binary     = measures.measure_percolation(phase)
    spans_binary    = any(perc_binary.values())
    spans_per_winner = measures.measure_percolation_per_winner(phase, winners, Q)
    print(f"       per-winner spanning (any single winner): {spans_per_winner}")
    print(f"       binary-union spanning (phase alone): x={perc_binary['x']}, "
          f"y={perc_binary['y']}, z={perc_binary['z']} -> any={spans_binary}")
    print(f"       done in {time.time() - t0:.1f}s")

    # Save
    print(f"\nWriting VTI files to {OUT.relative_to(HERE.parent)}/...")
    t0 = time.time()
    # NOTE: g1 field not re-saved (would need to keep fs around) — skip for memory
    io.save_vti(OUT / "margin.vti",             margin,                      grid, "margin")
    io.save_vti(OUT / "phase.vti",              phase.astype(np.uint8),      grid, "phase")
    io.save_vti(OUT / "winners.vti",            winners.astype(np.int16),    grid, "winner_label")
    io.save_vti(OUT / "components.vti",         components,                  grid, "component_id")
    io.save_vti(OUT / "components_physical.vti", components_phys,            grid, "phys_component_id")
    print(f"       done in {time.time() - t0:.1f}s")

    # Dump a one-line summary
    with (OUT / "summary.txt").open("w") as f:
        f.write(f"Scene: q={Q}, phi={PHI} (measured {measured_vf:.4f}), "
                f"N={N}, L={LENGTH}, lc={LC}\n")
        f.write(f"delta = {delta:.5f} physical = {film_voxels:.3f} voxels\n")
        f.write(f"per-winner components: {n_per_winner}\n")
        f.write(f"physical binary-union components: {n_physical}\n")
        f.write(f"largest physical component holds {largest_frac*100:.2f}% "
                f"of inclusion volume\n")
        f.write(f"per-winner spanning: {spans_per_winner}\n")
        f.write(f"binary-union spanning: x={perc_binary['x']}, "
                f"y={perc_binary['y']}, z={perc_binary['z']}\n")
    print(f"\nSummary at {OUT.relative_to(HERE.parent) / 'summary.txt'}")
    print(f"Total files: {sum(1 for _ in OUT.iterdir())}")


if __name__ == "__main__":
    main()
