"""HDF5 I/O for random fields and excursion sets.

Provides compact, compressed storage for scalar fields and binary excursion
sets.  All data is stored as float32 / uint8 with gzip compression (level 4),
giving 5–10× size reduction over uncompressed VTK for smooth Gaussian fields.

For ParaView visualisation, use io.save_vti / io.save_vtp.
For storing large grids or result archives, use these HDF5 functions.

XDMF note: to open an HDF5 field in ParaView, write a companion .xdmf file
(see write_xdmf_companion).
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

from .fields import GridSpec


# ---------------------------------------------------------------------------
# Scalar field
# ---------------------------------------------------------------------------

def save_field_h5(
    path: str | Path,
    field: np.ndarray,
    grid: GridSpec,
    name: str = "field",
    compression_level: int = 4,
) -> None:
    """Save a 3D scalar field to HDF5 (float32, gzip compressed).

    Metadata stored as attributes: n, length, dx.

    Args:
        path:              Output .h5 path.
        field:             3D array of shape (n, n, n).
        grid:              GridSpec (stores n, length, dx as attributes).
        name:              Dataset name inside the file.
        compression_level: gzip level 0–9 (4 is a good speed/size trade-off).
    """
    with h5py.File(path, "w") as f:
        f.create_dataset(
            name,
            data=field.astype(np.float32),
            compression="gzip",
            compression_opts=compression_level,
        )
        f.attrs["n"]      = grid.n
        f.attrs["length"] = grid.length
        f.attrs["dx"]     = grid.dx


def load_field_h5(path: str | Path, name: str = "field") -> tuple[np.ndarray, GridSpec]:
    """Load a scalar field and its GridSpec from an HDF5 file.

    Returns:
        (field array as float32, GridSpec).
    """
    with h5py.File(path, "r") as f:
        field  = f[name][:]
        n      = int(f.attrs["n"])
        length = float(f.attrs["length"])
    return field.astype(np.float32), GridSpec(n=n, length=length)


# ---------------------------------------------------------------------------
# Excursion set (binary phase + optional winner labels)
# ---------------------------------------------------------------------------

def save_excursion_h5(
    path: str | Path,
    phase: np.ndarray,
    grid: GridSpec,
    winners: np.ndarray | None = None,
    compression_level: int = 6,
) -> None:
    """Save a binary excursion set (and optional winner-type labels) to HDF5.

    Binary data compresses very well (run-length-encoding via gzip on smooth
    morphologies typically achieves 10–30× reduction).

    Args:
        path:    Output .h5 path.
        phase:   Boolean 3D array (True = included voxel).
        grid:    GridSpec.
        winners: Integer array of winner-type labels (0..q-1), same shape as
                 phase. Stored as int8 (supports up to 127 winner types).
    """
    with h5py.File(path, "w") as f:
        f.create_dataset(
            "phase",
            data=phase.astype(np.uint8),
            compression="gzip",
            compression_opts=compression_level,
        )
        if winners is not None:
            f.create_dataset(
                "winners",
                data=winners.astype(np.int8),
                compression="gzip",
                compression_opts=compression_level,
            )
        f.attrs["n"]      = grid.n
        f.attrs["length"] = grid.length
        f.attrs["dx"]     = grid.dx


def load_excursion_h5(
    path: str | Path,
) -> tuple[np.ndarray, np.ndarray | None, GridSpec]:
    """Load an excursion set from HDF5.

    Returns:
        (phase as bool array, winners as int8 array or None, GridSpec).
    """
    with h5py.File(path, "r") as f:
        phase   = f["phase"][:].astype(bool)
        winners = f["winners"][:] if "winners" in f else None
        n       = int(f.attrs["n"])
        length  = float(f.attrs["length"])
    return phase, winners, GridSpec(n=n, length=length)


# ---------------------------------------------------------------------------
# XDMF companion for ParaView
# ---------------------------------------------------------------------------

def write_xdmf_companion(h5_path: str | Path, dataset_name: str = "field") -> Path:
    """Write a minimal XDMF file that lets ParaView open an HDF5 scalar field.

    The .xdmf file is written alongside the .h5 file.  Open the .xdmf in
    ParaView; it will read the heavy data from the .h5.

    Returns:
        Path of the written .xdmf file.
    """
    h5_path = Path(h5_path)
    xdmf_path = h5_path.with_suffix(".xdmf")

    with h5py.File(h5_path, "r") as f:
        n      = int(f.attrs["n"])
        dx     = float(f.attrs["dx"])

    xdmf = f"""\
<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="3.0">
  <Domain>
    <Grid Name="mesh" GridType="Uniform">
      <Topology TopologyType="3DCoRectMesh" Dimensions="{n} {n} {n}"/>
      <Geometry GeometryType="ORIGIN_DXDYDZ">
        <DataItem Dimensions="3" NumberType="Float" Format="XML">0 0 0</DataItem>
        <DataItem Dimensions="3" NumberType="Float" Format="XML">{dx} {dx} {dx}</DataItem>
      </Geometry>
      <Attribute Name="{dataset_name}" AttributeType="Scalar" Center="Node">
        <DataItem Dimensions="{n} {n} {n}" NumberType="Float" Precision="4" Format="HDF">
          {h5_path.name}:/{dataset_name}
        </DataItem>
      </Attribute>
    </Grid>
  </Domain>
</Xdmf>
"""
    xdmf_path.write_text(xdmf)
    return xdmf_path


# ---------------------------------------------------------------------------
# Compression benchmark helper
# ---------------------------------------------------------------------------

def compression_stats(h5_path: str | Path) -> dict[str, float]:
    """Return file size and per-dataset compression ratios.

    Returns:
        Dict mapping dataset_name -> compression_ratio (uncompressed/compressed).
        Plus key 'file_mb' for total file size in MB.
    """
    h5_path = Path(h5_path)
    stats: dict[str, float] = {"file_mb": h5_path.stat().st_size / 1e6}
    with h5py.File(h5_path, "r") as f:
        for name, ds in f.items():
            if not isinstance(ds, h5py.Dataset):
                continue
            uncompressed = ds.size * ds.dtype.itemsize
            compressed   = ds.id.get_storage_size()
            stats[name]  = uncompressed / max(compressed, 1)
    return stats
