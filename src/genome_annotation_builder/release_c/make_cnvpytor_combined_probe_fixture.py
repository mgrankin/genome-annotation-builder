"""Create a tiny synthetic HDF5 fixture for the combined CNVpytor probe smoke."""

from __future__ import annotations

import argparse
from pathlib import Path

from .cnvpytor_combined_hdf5_probe import validate_combined_row


def make_fixture(output: Path) -> None:
    try:
        import h5py
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("h5py and numpy are required in the probe-control runtime") from exc
    row = validate_combined_row(
        [
            20.0,
            -1.0,
            100001.0,
            160000.0,
            60000.0,
            0.5,
            1e-8,
            0.8,
            0.1,
            0.1,
            0.1,
            60.0,
            0.2,
            1e-8,
            1e-6,
            0.0,
            20.0,
            5.0,
            0.0,
            0.0,
            0.9,
            1.0,
            1.0,
            0.0,
            0.9,
            1.0,
            2.0,
            1.0,
            1.0,
            0.1,
            1.0,
            *([0.0] * 40),
        ]
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output, "w") as destination:
        destination.attrs["Version"] = "1.3.2"
        destination.create_dataset("rd_chromosomes", data=np.array([b"1"]))
        destination.create_dataset("snp_chromosomes", data=np.array([b"1"]))
        for name in (
            "1_rd_p",
            "1_rd_u",
            "his_rd_p_1_1000_GC_mask",
            "his_rd_p_1_1000",
            "his_rd_u_1_1000",
            "snp_bafc_01_1_1000_mask",
            "snp_bafc_10_1_1000_mask",
            "snp_bafc_11_1_1000_mask",
            "snp_baf_1_1000_mask",
        ):
            destination.create_dataset(name, data=np.array([1.0], dtype="float64"))
        destination.create_dataset(
            "calls_2d_1_1000_mask_GC_mask",
            data=np.array([row], dtype="float64"),
        )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    make_fixture(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
