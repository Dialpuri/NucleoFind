import os
import tempfile
import nucleofind.predict as p
import pytest
import hashlib
from pathlib import Path
from types import SimpleNamespace


@pytest.fixture(scope='session')
def data_base_path():
    return Path(__file__).parent / "test_data" / "5d5w"


@pytest.fixture(scope='session')
def parameters(data_base_path):
    map_types = ["phosphate", "sugar", 'base']
    model_paths = [p.find_model(m) for m in map_types]
    mtzin = data_base_path / "hklout.mtz"
    colinfc = "FWT,PHWT"

    with tempfile.TemporaryDirectory() as tmp_directory:
        tmp_directory = Path(tmp_directory)
        yield SimpleNamespace(
            model_paths=model_paths,
            map_types=map_types,
            mtzin=str(mtzin),
            colinfc=colinfc,
            output=tmp_directory / "prediction"
        )


@pytest.fixture(scope='session')
def md5sums():
    phosphate_map_md5sum = "f9f3470859cac6c976b5146a635c58e7"
    sugar_map_md5sum = "ffa9190ab0f9d4c7316fa101a7ab3e10"
    base_map_md5sum = "3fa375bb2435def19020088dfd45b490"
    return SimpleNamespace(
        phosphate=phosphate_map_md5sum,
        sugar=sugar_map_md5sum,
        base=base_map_md5sum
    )


@pytest.fixture(scope='session')
def predictions_python(parameters):
    """
    Run NucleoFind using the Python API.

    Args:
        parameters (dict): A dictionary containing the necessary input parameters.

    Returns:
        str: The path to the output file.
    """
    output = parameters.output.parent / ("python_" + parameters.output.stem)
    prediction = p.Prediction(model_paths=parameters.model_paths)
    prediction.make_prediction(parameters.mtzin, parameters.colinfc.split(","), overlap=32)
    prediction.save_predicted_map(str(output))
    return output


@pytest.fixture(scope='session')
def predictions_cmdline(parameters):
    """
    Run NucleoFind using the command line interface.

    Args:
        parameters (dict): A dictionary containing the necessary input parameters.

    Returns:
        str: The path to the output file.
    """
    output = parameters.output.parent / ("cmd_" + parameters.output.stem)
    cmd = f"nucleofind -i {parameters.mtzin} -o {output} -m all -overlap 32"
    os.system(cmd)
    return output


def test_python_prediction(predictions_python, parameters, md5sums):
    """
    This function is used to test the predictions made by a Python model. It compares the MD5 sums of the predictions
    to the known MD5 sums.

    Parameters:
    - predictions_python (str): Path to the output model base directory.
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
       - predictions_python (str): Path to the output model base directory.
       - parameters (SimpleNamespace): A dictionary of parameters used in the predictions.
       - md5sums (SimpleNamespace): A dictionary of MD5 sums corresponding to the predictions.

       Raises:
           AssertionError: An error occurs if the calculated prediction MD5 sums do not equal the known MD5 sums.
       """
    compare_sums(md5sums, parameters, predictions_cmdline)


def compare_sums(md5sums, parameters, base_path):
    for map_type in parameters.map_types:
        output_map_path = base_path.parent / (base_path.name + f"_{map_type}.map")
        assert output_map_path.exists()

        with open(output_map_path, "rb") as f:
            file_hash = hashlib.md5()
            while chunk := f.read(8192):
                file_hash.update(chunk)

            found_sum = file_hash.hexdigest()
        assert found_sum == getattr(md5sums, map_type)
