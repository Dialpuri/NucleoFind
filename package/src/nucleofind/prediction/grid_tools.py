import gemmi
import numpy as np
from typing import Tuple, List

from .config import Configuration


def interpolate_grid(
    grid: gemmi.FloatGrid, configuration: Configuration
) -> Tuple[np.ndarray, gemmi.Transform]:
    """Interpolate grid to 0.7A grid spacing surrounding the unit cell and return interpolated grid and transform."""
    if configuration.compute_entire_unit_cell:
        extent = gemmi.FractionalBox()
        extent.extend(gemmi.Fractional(0, 0, 0))
        extent.extend(gemmi.Fractional(1, 1, 1))
    else:
        extent = gemmi.find_asu_brick(grid.spacegroup).get_extent()

    box = grid.unit_cell.orthogonalize_box(extent)
    margin = configuration.spacing * (configuration.box_size // 2)
    box.add_margin(margin)
    size = box.get_size()
    numx = -(
        -int(size.x / configuration.spacing)
        // configuration.overlap
        * configuration.overlap
    )
    numy = -(
        -int(size.y / configuration.spacing)
        // configuration.overlap
        * configuration.overlap
    )
    numz = -(
        -int(size.z / configuration.spacing)
        // configuration.overlap
        * configuration.overlap
    )
    array = np.zeros((numx, numy, numz), dtype=np.float32)
    scale = gemmi.Mat33(configuration.spacing * np.eye(3))
    transform: gemmi.Transform = gemmi.Transform(scale, box.minimum)
    grid.interpolate_values(array, transform)
    return array, transform


def precompute_slices(grid_shape: np.ndarray, overlap: int = 16) -> List[List[int]]:
    """Precompute indices of slices to run inference on."""
    slices = []

    for i in range(0, grid_shape[0] - overlap, overlap):
        for j in range(0, grid_shape[1] - overlap, overlap):
            for k in range(0, grid_shape[2] - overlap, overlap):
                slices.append([i, j, k])
    return slices


def reinterpolate_grid(
    work_array: np.ndarray,
    transform: gemmi.Transform,
    template_grid: gemmi.FloatGrid,
    compute_entire_unit_cell: bool = True,
) -> gemmi.FloatGrid:
    """Reinterpolate grid to original unit cell."""

    output_grid = gemmi.FloatGrid()
    output_grid.spacegroup = template_grid.spacegroup
    output_grid.set_unit_cell(template_grid.unit_cell)

    grid_spacing = 0.7
    output_grid.set_size(*template_grid.shape)
    size_x = work_array.shape[0] * grid_spacing
    size_y = work_array.shape[1] * grid_spacing
    size_z = work_array.shape[2] * grid_spacing

    array_cell = gemmi.UnitCell(size_x, size_y, size_z, 90, 90, 90)
    array_grid = gemmi.FloatGrid(work_array, array_cell)

    if compute_entire_unit_cell:
        grid_iterable = output_grid
    else:
        grid_iterable = output_grid.masked_asu()

    for point in grid_iterable:
        position = output_grid.point_to_position(point) - gemmi.Position(transform.vec)
        point.value = array_grid.interpolate_value(position)

    output_grid.symmetrize_max()
    return output_grid
