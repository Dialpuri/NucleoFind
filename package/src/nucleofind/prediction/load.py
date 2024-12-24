from pathlib import Path
import gemmi
import numpy as np
from typing import List
from .util import find_map_coefficients, check_density_path
import onnxruntime as rt
import sys
import logging


def load_mtz(
    path: Path | str,
    column_names: List[str] | None,
    resolution_cutoff: float | None,
) -> gemmi.FloatGrid:
    """Load MTZ file and transform to map with 0.7A grid spacing and with resolution cutoff if specified."""
    mtz = gemmi.read_mtz_file(str(path))
    if None in column_names:
        logging.warning(
            "No map coefficients were specified, NucleoFind will try and find some but they may be wrong."
        )
        column_names = find_map_coefficients(mtz)

    res = mtz.resolution_high()
    spacing = 0.7
    sample_rate = res / spacing
    grid = mtz.transform_f_phi_to_map(*column_names, sample_rate=sample_rate)
    grid.normalize()
    if resolution_cutoff:
        data = np.array(mtz, copy=False)
        mtz.set_data(data[mtz.make_d_array() >= resolution_cutoff])
    return grid


def load_map(path: Path | str) -> gemmi.FloatGrid:
    """Load map file and normalize"""
    map = gemmi.read_ccp4_map(str(path))
    grid = map.grid
    grid.normalize()
    return grid


def load_density(
    density_path: Path | str,
    column_names: List[str] | None,
    resolution_cutoff: float | None,
) -> gemmi.FloatGrid:
    """Load density from MTZ file, or map file"""
    density_path = check_density_path(density_path)

    if density_path.suffix == ".mtz":
        return load_mtz(density_path, column_names, resolution_cutoff)
    else:
        return load_map(density_path)


def load_onnx_model(
    model_path: Path | str, use_gpu: bool = True
) -> rt.InferenceSession:
    """Load ONNX model from model_path"""
    providers = ["CPUExecutionProvider"]
    if use_gpu:
        providers.insert(0, "CUDAExecutionProvider")
    sess_options = rt.SessionOptions()
    sess_options.intra_op_num_threads = 1
    try:
        return rt.InferenceSession(
            model_path, providers=providers, sess_options=sess_options
        )
    except OSError as e:
        logging.critical(
            "This model is corrupted, perhaps due to an incomplete download. Try downloading it again with "
            "nucleofind-install -m TYPE --reinstall"
        )
        sys.exit(1)
