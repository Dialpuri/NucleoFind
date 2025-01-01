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
def expected_md5sums(model_type):
    data = {
        "nano": SimpleNamespace(
            phosphate = "38cbdd8aaf54dacfb37a06b45dd48b6b",
            sugar = "b30c9d392ad441ef80fb3914d155310f",
            base = "43be8bbbef2b691cb7556c6536b9baae"
        ),
        "core": SimpleNamespace(
            phosphate = "478a977c95544a7e616b3148d822c80a",
            sugar = "34b4ac22a3331b39783bb8013349a707",
            base = "dea227621d152fdb36280f986284421a"
        ),
        "ultra": SimpleNamespace(
            phosphate = "946df04305b3c46b43a6c1a00def6223",
            sugar = "15a7808c90210aff7f0707e0429a6db4",
            base = "bfeb128562d95f6804cc6c4cf424f2b6"
        )

    }
    return data[model_type]

@pytest.fixture(scope='session', params=["nano", "core", "ultra"])
def model_type(request):
    return request.param

@pytest.fixture(scope='session')
def parameters(data_base_path, model_type):
    mtzin = data_base_path / "hklout.mtz"
    colinfc = "FWT,PHWT"

    with tempfile.TemporaryDirectory() as tmp_directory:
        tmp_directory = Path(tmp_directory)
        yield SimpleNamespace(
            model_name=model_type,
            map_types=["phosphate", "sugar", "base"],
            mtzin=str(mtzin),
            colinfc=colinfc,
            output=tmp_directory / "prediction"
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
    cmd = f'nucleofind -i "{parameters.mtzin}" -o "{output}" -m {parameters.model_name} -no-symmetry'
    os.system(cmd)
    return output


def test_python_prediction(predictions_python, parameters, expected_md5sums):
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
    compare_sums(expected_md5sums, parameters, predictions_python)


def test_cmdline_prediction(predictions_cmdline, parameters, expected_md5sums):
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
    compare_sums(expected_md5sums, parameters, predictions_cmdline)


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
