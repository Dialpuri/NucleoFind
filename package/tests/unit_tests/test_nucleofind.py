import hashlib
import nucleofind.predict as p
import pytest
from pathlib import Path
from types import SimpleNamespace
import tempfile
import gemmi


@pytest.fixture(scope='session')
def data_base_path():
    return Path(__file__).parent.parent / "test_data" / "5d5w"


@pytest.fixture(scope='session')
def parameters(data_base_path):
    map_types = ["phosphate"]
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
            output=tmp_directory / "prediction",
            unit_test_dir=data_base_path / "unit_test_data",
        )


@pytest.fixture(scope="session")
def predictions(parameters):
    return p.Prediction(model_paths=parameters.model_paths)


@pytest.fixture(scope="session", params=[{"cutoff": None}, {"cutoff": 4}])
def resolution_cutoff(request):
    return request.param["cutoff"]


# Test kwargs
def test_class_kwargs():
    model_paths = ["test1.onnx", "test2.onnx"]
    test_case_1 = p.Prediction(model_paths=model_paths)
    assert test_case_1.model_paths == model_paths

    test_case_2 = p.Prediction(model_paths=model_paths, use_gpu=False)
    assert not test_case_2.use_gpu

    test_case_3 = p.Prediction(model_paths=model_paths, compute_variance=True)
    assert test_case_3.compute_variance

    test_case_4 = p.Prediction(model_paths=model_paths, disable_progress_bar=True)
    assert test_case_4.disable_progress_bar

# Test load_mtz
@pytest.fixture(scope="session")
def mtz(predictions, parameters, resolution_cutoff):
    predictions._load_mtz(parameters.mtzin, resolution_cutoff=resolution_cutoff)
    return predictions.raw_grid


@pytest.fixture(scope="session")
def computed_hash(resolution_cutoff):
    data = {
        None: b'\xab\xcfU\x08b\xc6k\x13#-"\xe2\xcclj`',
        4: b'a0\xf5\x95\xa1~\xda\x87\x1b\xd2\xdf\x0bf\x85{\x0c'
    }
    return data[resolution_cutoff]


def compare_hashes(computed_hash, mtz):
    array_hash = hashlib.md5(mtz.array.tobytes()).digest()
    assert array_hash == computed_hash


@pytest.mark.unit
def test_load_mtz(parameters, mtz, resolution_cutoff, computed_hash):
    compare_hashes(computed_hash, mtz)


# Test get_bounding_box
@pytest.fixture(scope="session")
def expected_box_parameters(resolution_cutoff):
    data = {
        None: {
            "maximum": gemmi.Position(39.771, 78.397, 90.084),
            "minimum": gemmi.Position(0, 0, 0)
        },
        4: {
            "maximum": gemmi.Position(39.771, 78.397, 90.084),
            "minimum": gemmi.Position(0, 0, 0)
        }
    }
    return data[resolution_cutoff]


def compare_positions(position_1: gemmi.Position, position_2: gemmi.Position):
    assert position_1.x == position_2.x
    assert position_1.y == position_2.y
    assert position_1.z == position_2.z


@pytest.mark.unit
def test_get_bounding_box(mtz, expected_box_parameters):
    box: gemmi.PositionBox = p.Prediction._get_bounding_box(mtz)
    compare_positions(box.maximum, expected_box_parameters["maximum"])
    compare_positions(box.minimum, expected_box_parameters["minimum"])


@pytest.fixture(scope="session")
def interpolated_grid(predictions, parameters):
    predictions._load_mtz(parameters.mtzin, resolution_cutoff=None)
    predictions._interpolate_grid()
    return predictions.interpolated_grid


# Test interpolate_grid
@pytest.mark.unit
def test_interpolate_grid(interpolated_grid):
    expected_interpolated_grid_hash = b'\x84\xb6\xf3\xd4(\xde\xb8\x00\xe3\xd1v^\xfd\x00\x87\x9a'
    compare_hashes(expected_interpolated_grid_hash, interpolated_grid)


# Test calculate_translations
def test_calculate_translations(predictions, interpolated_grid):
    predictions._calculate_translations(overlap=32)

    na, nb, nc = predictions.na, predictions.nb, predictions.nc
    assert na == 2
    assert nb == 4
    assert nc == 5
    assert len(predictions.translation_list) == 40

    predictions._calculate_translations(overlap=16)
    assert len(predictions.translation_list) == 292


#Test predict
# def test_predict(predictions):
#     pass
