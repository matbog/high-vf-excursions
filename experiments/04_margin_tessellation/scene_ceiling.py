"""Export one representative ceiling-push realization to ParaView.

Picks q = 32 at target phi = 0.75 on an N = 128 grid — a point from the
ceiling push's "safe regime" (never spans across 10 realizations in
Phase A). Writes a complete set of VTI files for ParaView inspection:

    field_g1.vti    first Gaussian competitor G_1 (shows the underlying field)
    margin.vti      max_i G_i - second-max_i G_i (threshold at delta -> phase)
    phase.vti       binary inclusion mask after thresholding
    winners.vti     integer winning-field index, 0..q-1 (colour per cell)
    components.vti  per-winner connected-component ID (unique int per aggregate)

In ParaView:
    - Open field_g1.vti alone for a 3D volume render of the Gaussian field.
    - Open phase.vti and threshold at 0.5 to isolate the aggregate phase.
    - Open winners.vti and colour the phase by winner_label.
    - Open components.vti and colour by component_id to see each individual
      aggregate as a distinct object. The total count is printed here.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from skimage import measure

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[1]))

from src import adler, excursions, fields, io  # noqa: E402
from src.fields import GridSpec  # noqa: E402

N      = 128
LENGTH = 64.0
LC     = 5.0
Q      = 32
PHI    = 0.75
SEED   = 1729 + 10000 * Q   # same scheme as run_ceiling.py realization 0

OUT = HERE / "paraview" / f"q{Q}_phi{int(PHI*100):02d}"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> None:
    grid = GridSpec(n=N, length=LENGTH)
    rng = np.random.default_rng(SEED)

    print(f"Generating {Q} Gaussian competitors at N={N} (lc={LC})...")
    fs = fields.generate_margin_fields(grid, LC, Q, rng)
    print(f"  shape = {fs.shape}, mem = {fs.nbytes / 1e6:.0f} MB")

    delta = adler.delta_for_vf(Q, PHI)
    print(f"delta for phi={PHI}: {delta:.5f} physical  "
          f"({delta / grid.dx:.3f} voxels)")

    # Margin scalar
    top2   = np.partition(fs, kth=Q - 2, axis=0)[-2:]
    margin = (top2[1] - top2[0]).astype(np.float32)

    # Phase + per-voxel winner label
    phase, winners = excursions.margin_excursion(fs, delta)
    print(f"measured VF: {phase.mean():.4f}  (target {PHI})")

    # Per-winner connected components, with unique integer IDs across winners
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
    n_total = cur - 1
    print(f"per-winner connected components: {n_total} total across {Q} winners")

    io.save_vti(OUT / "field_g1.vti",   fs[0].astype(np.float32),  grid, "G_1")
    io.save_vti(OUT / "margin.vti",     margin,                     grid, "margin")
    io.save_vti(OUT / "phase.vti",      phase.astype(np.uint8),     grid, "phase")
    io.save_vti(OUT / "winners.vti",    winners.astype(np.int16),   grid, "winner_label")
    io.save_vti(OUT / "components.vti", components,                 grid, "component_id")
    print(f"\nSaved 5 VTI files to {OUT.relative_to(HERE.parent)}/")
    print(f"  (phase binary + {n_total} labelled aggregates)")


if __name__ == "__main__":
    main()
