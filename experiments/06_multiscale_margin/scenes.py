"""Convert exp06 sample_vf_*.npz checkpoints to ParaView-readable VTI files.

For each target VF, the main run.py already saves the realization-1 binary
phase and class ownership map. This script just re-exports them to VTI format
so you can open them in ParaView without re-running the sweep.

Output: experiments/06_multiscale_margin/paraview/vf{VF}/{phase,class_owner}.vti

In ParaView:
  - phase.vti        → Threshold at 0.5 to isolate the aggregate phase.
  - class_owner.vti  → Threshold / colour-by to see coarse/fine/sand ownership
                       (0 = matrix, 1 = coarse_gravel, 2 = fine_gravel, 3 = sand_proxy).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2]))
from src import io
from src.fields import GridSpec

HERE   = Path(__file__).parent
RESULTS = HERE / "results"
OUT     = HERE / "paraview"
OUT.mkdir(exist_ok=True)


def export(npz_path: Path) -> None:
    data   = np.load(npz_path, allow_pickle=True)
    phase  = data["phase"]
    owner  = data["class_owner"]
    length = float(data["length"])
    dx     = float(data["dx"])
    phi    = float(data["target_phi"])
    class_names = list(data["class_names"])
    deltas      = data["deltas"]

    n    = phase.shape[0]
    grid = GridSpec(n=n, length=length)

    out_dir = OUT / f"vf{int(phi * 100):02d}"
    out_dir.mkdir(exist_ok=True)

    io.save_vti(out_dir / "phase.vti",       phase.astype(np.uint8), grid, "phase")
    io.save_vti(out_dir / "class_owner.vti", owner.astype(np.int16),  grid, "class_label")

    label_map = {0: "matrix", **{i + 1: name for i, name in enumerate(class_names)}}
    print(f"  φ={phi:.2f}  N={n}  dx={dx:.3f}  deltas={[f'{d:.3f}' for d in deltas]}")
    print(f"    class labels: {label_map}")
    print(f"    → {out_dir.relative_to(HERE.parent)}")


if __name__ == "__main__":
    npz_files = sorted(RESULTS.glob("sample_vf_*.npz"))
    if not npz_files:
        raise SystemExit(f"No sample_vf_*.npz files found in {RESULTS}. "
                         "Run run.py with save_volumes=true in config first.")
    print(f"Exporting {len(npz_files)} scenes from {RESULTS.name}/ to {OUT.name}/\n")
    for p in npz_files:
        export(p)
    print("\nDone. Open any phase.vti + class_owner.vti together in ParaView;")
    print("colour the thresholded phase by class_label for the material classes.")
