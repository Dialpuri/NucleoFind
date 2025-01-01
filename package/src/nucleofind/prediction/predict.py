import functools
import logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List, Tuple
import numpy as np
from tqdm import tqdm

from .load import load_density, load_onnx_model
from .grid_tools import interpolate_grid, reinterpolate_grid, precompute_slices
from .arguments import parse_arguments
from .model import find_model, get_model_config
from .save import save_grid
from ..logs import setup_logging
from .config import Configuration, MapType


class NucleoFind:
    def __init__(self, model_path: Path | str, configuration: Configuration):
        self.model_path = model_path
        self.configuration = configuration

        self.model = None
        self.predicted_grids = {}

    def _process_sample(
        self,
        input_name: str,
        output_shape: Tuple[int],
        box_size: int,
        array: np.ndarray[np.float32],
        translation: Tuple[int, int, int],
    ) -> Tuple[np.ndarray, Tuple[int, int, int]]:
        """Perform inference on a single sample of shape (1, box_size, box_size, box_size, 1) and return an array of shape
        (box_size, box_size, box_size, output_channels) and the translation (for putting back into an array).
        """
        i, j, k = translation
        input_sub = array[i : i + box_size, j : j + box_size, k : k + box_size]
        input_sub = input_sub[np.newaxis, ..., np.newaxis].astype(np.float16)

        return np.array(self.model.run(None, {input_name: input_sub})).reshape(
            output_shape
        ), translation

    def _run_prediction(self, work_grid: np.ndarray) -> np.ndarray:
        """Run prediction on work_grid and calculate the average predicted grid"""
        work_grid_shape = np.array(work_grid.shape)
        slices = precompute_slices(work_grid_shape, overlap=self.configuration.overlap)
        box_size = self.configuration.box_size

        total_array = np.zeros((*work_grid_shape, 4), dtype=np.float32)
        count_array = np.zeros_like(total_array, dtype=np.float32)

        # Variance arrays for Welch's one pass variance method
        variance_mean = np.zeros_like(total_array, dtype=np.float32)
        variance_m2 = np.zeros_like(total_array, dtype=np.float32)

        channels = self.configuration.channels
        input_name = self.model.get_inputs()[0].name
        output_shape = (box_size, box_size, box_size, channels)
        process_sample_worker = functools.partial(
            self._process_sample,
            input_name,
            output_shape,
            box_size,
            work_grid,
        )

        miniters = 1_000 if len(slices) > 10_000 else 1
        max_workers = self.configuration.n_threads
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = list(
                tqdm(
                    executor.map(process_sample_worker, slices),
                    total=len(slices),
                    desc="Predicting",
                    miniters=miniters,
                    disable=self.configuration.disable_progress_bar,
                )
            )

        ones = np.ones(channels)
        for result in tqdm(results, desc="Processing results"):
            predicted_sub, (i, j, k) = result
            box_slice = (
                slice(i, i + box_size),
                slice(j, j + box_size),
                slice(k, k + box_size),
                slice(None),
            )
            total_array[box_slice] += predicted_sub
            count_array[box_slice] += ones

            if self.configuration.compute_variance:
                delta_variance = total_array[box_slice] - variance_mean[box_slice]
                variance_mean[box_slice] += delta_variance / count_array[box_slice]
                variance_m2[box_slice] += delta_variance * (
                    total_array[box_slice] - variance_mean[box_slice]
                )

        predicted_array = total_array / count_array
        if self.configuration.use_raw_values:
            return predicted_array.astype(np.float32)

        if self.configuration.compute_variance:
            variance_array = variance_m2 / (np.subtract(count_array, 1))
            return variance_array.astype(np.float32)

        argmax_array = np.argmax(predicted_array, axis=-1).squeeze()
        return argmax_array.astype(np.float32)

    def predict(
        self,
        density_path: Path | str,
        column_names: List[str] | List[None],
        resolution_cutoff: float | None = None,
    ):
        """Run a nucleofind prediction on specified density file. If density file is an MTZ, supply column names and an
        optional resolution cutoff. If density file is a MAP, these will be ignored."""
        self.model = load_onnx_model(self.model_path, self.configuration.use_gpu)
        input_grid = load_density(density_path, column_names, resolution_cutoff)

        work_grid, transform = interpolate_grid(input_grid, self.configuration)
        predicted_array = self._run_prediction(work_grid)

        rounded_array = np.round(predicted_array)
        for i in range(1, self.configuration.channels):
            if self.configuration.use_raw_values or self.configuration.compute_variance:
                index_array = predicted_array[:, :, :, i].astype(np.float32)
            else:
                index_array = (rounded_array == i).astype(np.float32)

            interpolated_index_array = reinterpolate_grid(
                index_array,
                transform,
                input_grid,
                self.configuration.compute_entire_unit_cell,
            )
            self.predicted_grids[MapType(i)] = interpolated_index_array

    def save_grid(self, type: MapType, output_path: Path | str):
        """Save the predicted grid to directory specified by output_path, with filename nucleofind-{type}.map."""
        output_path = Path(output_path)
        output_path.mkdir(exist_ok=True, parents=True)
        logging.info(f"Saving grid of type {type} to {output_path}")

        suffix = ".map"
        suffix = ".variance.map" if self.configuration.compute_variance else suffix
        suffix = ".raw.map" if self.configuration.use_raw_values else suffix

        if type == MapType.all:
            for k, v in self.predicted_grids.items():
                save_grid(v, output_path / f"nucleofind-{k.name}{suffix}")
            return
        elif type not in self.predicted_grids:
            raise ValueError(f"No grid of type {type} found.")

        save_grid(
            self.predicted_grids[type],
            output_path / f"nucleofind-{type.name}{suffix}",
        )


def run():
    """Run prediction from command line arguments"""
    setup_logging()
    args = parse_arguments()
    model_path = find_model(args.m)
    model_configuration = get_model_config(model_path, args.overlap)
    configuration = Configuration(
        use_gpu=args.gpu,
        disable_progress_bar=args.silent,
        compute_entire_unit_cell=False,
        use_raw_values=args.raw,
        compute_variance=args.variance,
        n_threads=args.n,
        **vars(model_configuration),
    )
    nucleofind = NucleoFind(model_path, configuration)
    nucleofind.predict(
        args.i,
        [args.amplitude, args.phase],
    )
    output_dir = Path(args.o)
    nucleofind.save_grid(MapType.phosphate, output_dir)
    nucleofind.save_grid(MapType.sugar, output_dir)
    nucleofind.save_grid(MapType.base, output_dir)


def predict_map(
    model: str,
    input: str,
    output: str,
    resolution: float = None,
    amplitude: str = "FWT",
    phase: str = "PHWT",
    overlap: int = None,
):
    """Run prediction from Python"""
    logging.info(
        f"Running prediction with model {model}, input {input}, output {output}, resolution {resolution}, amplitude {amplitude}, phase {phase}, overlap {overlap}"
    )
    model_path = find_model(model)
    model_configuration = get_model_config(model_path, overlap)
    configuration = Configuration(
        use_gpu=False,
        disable_progress_bar=False,
        compute_entire_unit_cell=False,
        n_threads=None,
        **vars(model_configuration),
    )
    prediction = NucleoFind(model_path, configuration=configuration)
    prediction.predict(input, [amplitude, phase], resolution_cutoff=resolution)
    prediction.save_grid(MapType.phosphate, output)
    prediction.save_grid(MapType.sugar, output)
    prediction.save_grid(MapType.base, output)
