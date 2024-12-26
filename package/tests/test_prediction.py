import os
import tempfile
import nucleofind.prediction.predict as p
import nucleofind.prediction.config as config
import pytest
import hashlib
from pathlib import Path
from types import SimpleNamespace

from nucleofind.prediction.config import MapType


@pytest.fixture(scope='session')
def data_base_path():
    return Path(__file__).parent / "test_data" / "5d5w"


@pytest.fixture(scope='session')
def parameters(data_base_path):
    # map_types = ["nano", "core"]
    model_name = "core"
    mtzin = data_base_path / "hklout.mtz"
    colinfc = "FWT,PHWT"

    with tempfile.TemporaryDirectory(delete=False) as tmp_directory:
        tmp_directory = Path(tmp_directory)
        yield SimpleNamespace(
            model_name=model_name,
            map_types=["phosphate", "sugar", "base"],
            mtzin=str(mtzin),
            colinfc=colinfc,
            output=tmp_directory / "prediction"
            )


@pytest.fixture(scope='session')
def md5sums():
    phosphate_map_md5sum = "478a977c95544a7e616b3148d822c80a"
    sugar_map_md5sum = "34b4ac22a3331b39783bb8013349a707"
    base_map_md5sum = "7976c788c4a76222898b02b1d58e8b9a"
    return SimpleNamespace(
        phosphate=phosphate_map_md5sum,
        sugar=sugar_map_md5sum,
        base=base_map_md5sum
    )


@pytest.fixture(scope='session')
def predictions_python(parameters) -> Path:
    """
    Run NucleoFind using the Python API.

    Args:
        parameters (SimpleNamespace): A dictionary containing the necessary input parameters.

    Returns:
        str: The path to the output file.
    """
    output = parameters.output.parent / ("python_" + parameters.output.stem)
    fwt, phwt = parameters.colinfc.split(",")
    p.predict_map(parameters.model_name, parameters.mtzin, output, amplitude=fwt, phase=phwt)
    return output


@pytest.fixture(scope='session')
def predictions_cmdline(parameters) -> Path:
    """
    Run NucleoFind using the command line interface.

    Args:
        parameters (SimpleNamespace): A dictionary containing the necessary input parameters.

    Returns:
        str: The path to the output file.
    """
    output = parameters.output.parent / ("cmd_" + parameters.output.stem)
    cmd = f'nucleofind -i "{parameters.mtzin}" -o "{output}" -m core -no-symmetry'
    os.system(cmd)
    return output


def test_python_prediction(predictions_python, parameters, md5sums):
    """
    This function is used to test the predictions made by a Python model. It compares the MD5 sums of the predictions
    to the known MD5 sums.

    Parameters:
    - predictions_python (path): Path to the output model base directory.
    - parameters (SimpleNamespace): A dictionary of parameters used in the predictions.
    - md5sums (SimpleNamespace): A dictionary of MD5 sums corresponding to the predictions.

    Raises:
        AssertionError: An error occurs if the calculated prediction MD5 sums do not equal the known MD5 sums.
    """
    compare_sums(md5sums, parameters, predictions_python)


def test_cmdline_prediction(predictions_cmdline, parameters, md5sums):
    """
    This function is used to test the predictions made by a command line interface. It compares the MD5 sums of the
    predictions to the known MD5 sums.

    Parameters:
    - predictions_python (path): Path to the output model base directory.
    - parameters (SimpleNamespace): A dictionary of parameters used in the predictions.
    - md5sums (SimpleNamespace): A dictionary of MD5 sums corresponding to the predictions.

    Raises:
        AssertionError: An error occurs if the calculated prediction MD5 sums do not equal the known MD5 sums.
       """
    compare_sums(md5sums, parameters, predictions_cmdline)


def compare_sums(md5sums, parameters, base_path):
    for map_type in parameters.map_types:
        output_map_path = base_path / f"nucleofind-{map_type}.map"
        assert output_map_path.exists()

        with open(output_map_path, "rb") as f:
            file_hash = hashlib.md5()
            while chunk := f.read(8192):
                file_hash.update(chunk)

            found_sum = file_hash.hexdigest()
        assert found_sum == getattr(md5sums, map_type)
