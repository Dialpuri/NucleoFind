import gemmi
from pathlib import Path


def save_grid(grid: gemmi.FloatGrid, path: Path | str):
    """Save grid to CCP4 map file."""
    map = gemmi.Ccp4Map()
    map.grid = grid
    map.update_ccp4_header()
    map.write_ccp4_map(str(path))
