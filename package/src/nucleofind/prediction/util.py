import logging
from pathlib import Path

import gemmi
from typing import Tuple
import sys


def find_map_coefficients(mtz: gemmi.Mtz) -> Tuple[str, str]:
    """Find F and P columns in MTZ file."""
    Fs = mtz.columns_with_type("F")
    Ps = mtz.columns_with_type("P")
    Fs = [F.label for F in Fs]
    Ps = [P.label for P in Ps]

    if not Fs or not Ps:
        logging.critical("No F and P columns found in MTZ file.")
        sys.exit(1)

    if "FWT" in Fs and "PHWT" in Ps:
        logging.warning(f"FWT and PHWT found, using them.")
        return "FWT", "PHWT"

    F, P = Fs[0], Ps[0]
    if len(Fs) != 1 or len(Ps) != 1:
        logging.warning(f"Multiple F and P columns found. Using first set. {F=}, {P=}")
    return F.label, P.label


def check_density_path(density_path):
    """Check that density path is a valid type and exists."""
    allowed_extensions = [".mtz", ".map", ".ccp4", ".mrc", ".gz"]
    density_path = Path(density_path)
    if not density_path.exists():
        logging.critical(f"Density file {density_path} does not exist.")
        sys.exit(1)

    if any(suffix not in allowed_extensions for suffix in density_path.suffixes):
        logging.critical(f"Density file must be one of {allowed_extensions}")
        sys.exit(1)
    return density_path
