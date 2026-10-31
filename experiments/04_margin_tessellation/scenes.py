"""Regenerate representative margin-tessellation realizations and export to
ParaView-readable formats (VTI + VTP).

For each (q, phi) scene we save:
  margin_scalar.vti   max_i G_i − 2nd-max_i G_i      (threshold at δ in ParaView to get phase)
  winners.vti         integer winner index 0..q-1    (colour per phase to see winners)
  phase.vti           binary mask of A_δ             (optional; = margin_scalar >= δ)
  surface_full.vtp    marching-cubes isosurface of the per-winner signed field (smooth)
                      — ONE file per winner, concatenated into a MultiBlock would be
                      cleaner but we keep per-winner files for ease of colouring.

Run:
    mamba run -n rf_explo python scenes.py

Output: experiments/04_margin_tessellation/paraview/q{Q}_phi{PHI}/
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2]))
from src import adler, excursions, fields, io
from src.fields import GridSpec

OUT_ROOT = Path(__file__).parent / "paraview"
OUT_ROOT.mkdir(exist_ok=True)

SEED   = 1729
N      = 128
LENGTH = 64.0
LC     = 5.0

# Representative scenes — chosen to show the morphology at the main
# validation resolution across a range of q and phi.
SCENES = [
    (2,  0.50),   # two chambers, loose film
    (8,  0.60),   # intermediate structure, well-resolved film (~1.7 vox)
    (8,  0.70),   # q=8 near percolation
    (16, 0.70),   # dense tessellation, film ~1.9 vox
    (16, 0.85),   # near sub-voxel limit
]


def export_scene(q: int, phi: float) -> None:
    grid = GridSpec(n=N, length=LENGTH)
    rng  = np.random.default_rng(SEED)
    fs   = fields.generate_margin_fields(grid, LC, q, rng)

    delta = adler.delta_for_vf(q, phi)
    phase, winners = excursions.margin_excursion(fs, delta)

    # Margin scalar field = max − 2nd_max. Threshold at δ in ParaView to recover phase.
    top2    = np.partition(fs, kth=q - 2, axis=0)[-2:]
    margin  = (top2[1] - top2[0]).astype(np.float32)

    out = OUT_ROOT / f"q{q}_phi{int(phi * 100):02d}"
    out.mkdir(exist_ok=True)

    io.save_vti(out / "margin_scalar.vti", margin,  grid, "margin")
    io.save_vti(out / "winners.vti",       winners.astype(np.int16), grid, "winner_label")
    io.save_vti(out / "phase.vti",         phase.astype(np.uint8),   grid, "phase")

    # Per-winner smooth isosurfaces: signed field s_i = G_i − max_{j≠i} G_j − δ.
    # {s_i >= 0} is winner i's cell. ParaView contour at level=0 gives the smooth
    # boundary (no voxel staircase).
    sec_field = top2[0]
    max_field = top2[1]
    for i in range(q):
        is_i_winner = (winners == i)
        max_other   = np.where(is_i_winner, sec_field, max_field)
        s_i         = (fs[i] - max_other - delta).astype(np.float32)
        io.save_vti(out / f"signed_winner_{i:02d}.vti", s_i, grid, "signed_margin")

    print(f"  q={q} phi={phi:.2f}  delta={delta:.4f}  →  {out.relative_to(OUT_ROOT.parent)}")


if __name__ == "__main__":
    print(f"Exporting ParaView scenes for exp04 to {OUT_ROOT.relative_to(Path.cwd())}")
    print(f"  (N={N}, LENGTH={LENGTH}, LC={LC}, SEED={SEED})\n")
    for q, phi in SCENES:
        export_scene(q, phi)
    print("\nIn ParaView:")
    print("  - Open margin_scalar.vti, Contour at level=δ (printed above) → the matrix film")
    print("  - Open winners.vti, Threshold between 0..q-1 → coloured winner cells")
    print("  - Open signed_winner_XX.vti, Contour at level=0 → smooth per-winner surface")
    print("  - Colour the signed-winner contours by the 'winner_label' field for grouping")
