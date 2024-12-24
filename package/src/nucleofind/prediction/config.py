import dataclasses
import enum


@dataclasses.dataclass
class Configuration:
    """Configuration for NucleoFind"""

    use_gpu: bool = False
    use_multiprocessing: bool = True
    disable_progress_bar: bool = True
    compute_entire_unit_cell: bool = True
    compute_variance: bool = False
    spacing: float = 0.7
    box_size = 32
    channels: int = 4
    overlap: int = 16


class MapType(enum.Enum):
    """Map types for NucleoFind, i.e. model will output 1 for phosphate..."""

    combined: int = 0
    phosphate: int = 1
    sugar: int = 2
    base: int = 3
    all: int = 4