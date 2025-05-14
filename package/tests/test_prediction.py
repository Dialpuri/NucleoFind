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
            phosphate = "615b7f3ee576dd19f2cd3ffc762013ff",
            sugar = "cf60e8d12058021046cad1a1f5e17e65",
            base = "a2737fd6cc3f08a12e468d888485af16"
        ),
        "core": SimpleNamespace(
            phosphate = "7351570fa904c36d806c4d1c20681c76",
            sugar = "fdc128322bed205b39eb7f9fc39ff8a5",
            base = "3ddd3206e340540e1471ef6f42132266"
        )
    }
    return data[model_type]

def pytest_generate_tests(metafunc):
    """Generate model matrix with runslow option"""
    if "model_type" in metafunc.fixturenames:
        runslow = metafunc.config.getoption("--runslow")
        params = ["nano", "core", "ultra"] if runslow else ["core"]
        metafunc.parametrize("model_type", params, ids=params, scope="session")

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
    cmd = f'nucleofind -i "{parameters.mtzin}" -o "{output}" -m {parameters.model_name} --no-use-symmetry'
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
