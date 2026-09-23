"""Convert WaveDump ROOT files to HDF5 for collaborators.

Run from the project root:
    python python/root_to_hdf5.py
"""

import os
import sys

import h5py

sys.path.insert(0, "/home/e-work/Software/Analysis-Utilities/python")
from analysis_utilities.io import load_tree_data

FILES = {
    "Cs137_raw": ("macros/root_files/Cs137_raw.root", "Data_R"),
    "Am241_raw": ("macros/root_files/Am241_raw.root", "Data_R"),
    "Cs137": ("macros/root_files/Cs137.root", "features"),
    "Am241": ("macros/root_files/Am241.root", "features"),
}


def convert_to_hdf5(root_path, tree_name, hdf5_path):
    """Convert a ROOT file to HDF5."""

    df, waveforms = load_tree_data(
        root_path,
        tree_name=tree_name,
        array_branch="Samples",
    )

    with h5py.File(hdf5_path, "w") as hf:
        for col in df.columns:
            hf.create_dataset(col, data=df[col].values, compression="gzip")
        if waveforms is not None:
            hf.create_dataset("Samples", data=waveforms, compression="gzip")

    n_events, n_samples = waveforms.shape if waveforms is not None else (len(df), 0)
    print(f"Converted {root_path} -> {hdf5_path}")
    print(f"  {n_events} events, {n_samples} samples/event")


if __name__ == "__main__":
    os.makedirs("hdf5_files", exist_ok=True)

    for name, (root_path, tree_name) in FILES.items():
        if not os.path.exists(root_path):
            print(f"Skipping {root_path} (not found)")
            continue
        convert_to_hdf5(root_path, tree_name, f"hdf5_files/{name}.hdf5")
