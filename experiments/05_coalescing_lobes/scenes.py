"""Regenerate representative coalescing-lobes realizations and export to
ParaView-readable formats (VTI).

For each scene (one radius r/r_c), we save:
  distance.vti    scalar field d(x) = min_a ||G(x) − c_a||.
                  ParaView Contour at level=r → smooth isosurface of A_r.
  nearest.vti     int label a(x) = argmin_a ||G(x) − c_a||, 0..K-1.
                  Colour per label to see the K lobe families.
  phase.vti       binary mask (d(x) <= r). Threshold at 0.5 to isolate the phase.

Run:
    mamba run -n rf_explo python scenes.py

Output: experiments/05_coalescing_lobes/paraview/rrc{R:.2f}/
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2]))
from src import excursions, fields, io
from src.fields import GridSpec

OUT_ROOT = Path(__file__).parent / "paraview"
OUT_ROOT.mkdir(exist_ok=True)

SEED    = 9417
N       = 100
LENGTH  = 100.0
LC      = 10.0
PADDING = 20
RHO     = 1.5

# Representative scenes: pre-contact (disjoint balls), at and after contact
SCENES = [0.65, 0.85, 0.98, 1.15]


def export_scene(ratio: float, distances: np.ndarray, nearest: np.ndarray,
                 grid: GridSpec, rc: float) -> None:
    radius = rc * ratio
    phase  = distances <= radius

    out = OUT_ROOT / f"rrc{ratio:.2f}".replace(".", "p")
    out.mkdir(exist_ok=True)

    io.save_vti(out / "distance.vti", distances.astype(np.float32), grid, "distance")
    io.save_vti(out / "nearest.vti",  nearest.astype(np.int16),    grid, "lobe_label")
    io.save_vti(out / "phase.vti",    phase.astype(np.uint8),      grid, "phase")

    print(f"  r/r_c={ratio:.2f}  r={radius:.4f}  phi={phase.mean():.3f}  →  {out.name}")


if __name__ == "__main__":
    print(f"Exporting ParaView scenes for exp05 to {OUT_ROOT.relative_to(Path.cwd())}")
    print(f"  (N={N}, LENGTH={LENGTH}, LC={LC}, RHO={RHO}, SEED={SEED})\n")

    grid    = GridSpec(n=N, length=LENGTH)
    centres = excursions.tetrahedral_centres(RHO)
    rc      = float(np.linalg.norm(centres[0] - centres[1]) / 2.0)
    print(f"  r_c = {rc:.4f}  (contact radius)")

    # One field, all radii use the same realisation to show the morphology evolve.
    vf = fields.generate_vector_field(grid, LC, 3, np.random.default_rng(SEED),
                                      padding=PADDING)
    distances, nearest = excursions.lobe_distances(vf, centres)

    for ratio in SCENES:
        export_scene(ratio, distances, nearest, grid, rc)

    print("\nIn ParaView:")
    print("  - distance.vti → Contour at level=r (printed above) → smooth isosurface of A_r")
    print("  - nearest.vti   → Threshold or colour-by to show the 4 lobe families")
    print("  - For post-contact (r/r_c=1.15), you'll see merged lobes")
