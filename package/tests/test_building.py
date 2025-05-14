import os
import nucleofind.build as nucleofind_build
import pytest
import xml.etree.ElementTree as ET
import time
import tempfile
from pathlib import Path


@pytest.fixture(scope='session')
def expected_values(test_example):
    data = {
        "1hr2": {
            "time": 100,
            "fragments_built": 10,
            "residues_built": 250,
            "residues_sequenced": 150,
            "longest_fragment": 45
        },
        "5d5w": {
            "time": 20,
            "fragments_built": 1,
            "residues_built": 5,
            "residues_sequenced": 4,
            "longest_fragment": 5
        }
    }
    return data[test_example]


@pytest.fixture(scope='session')
def data_base_path():
    return Path(__file__).parent / "test_data"


@pytest.fixture(scope='session', params=[
    {"pdb": "1hr2"},
    {"pdb": "5d5w"},
])
def test_example(request):
    return request.param['pdb']


@pytest.fixture(scope='session')
def parameters(data_base_path, test_example):
    input = nucleofind_build.InputParameters()
    data_dir = data_base_path / test_example
    input.mtzin = data_dir / "hklout.mtz"
    input.pdbin = data_dir / "xyzout.pdb"

    if not input.pdbin.exists():
        input.pdbin = ""

    input.seqin = data_dir / f"{test_example}.fasta"
    input.phosin = data_dir / "nucleofind-phosphate.map"
    input.sugarin = data_dir / "nucleofind-sugar.map"
    input.basein = data_dir / "nucleofind-base.map"
    input.colinfo = "FP,SIGFP"
    input.colinfc = "FWT,PHWT"
    input.colinfree = "FREE"

    os.environ["CLIBD"] = str(data_base_path)
    os.environ["LD_LIBRARY_PATH"] = str(data_base_path)
    os.environ["CCP4"] = str(data_base_path)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        output = nucleofind_build.OutputParameters()
        output.pdbout = tmpdir / "built_model.pdb"
        output.xmlout = tmpdir / "built_model.xml"
        yield input, output


@pytest.fixture(scope='session')
def build(parameters):
    input, output = parameters
    start = time.time()
    nucleofind_build.build(input, output)
    end = time.time()
    return end - start


@pytest.fixture(scope='session')
def xml(build, parameters):
    _, output = parameters
    tree = ET.parse(output.xmlout)
    root = tree.getroot()
    return root


def test_pdbout(build, parameters):
    _, output = parameters
    assert output.pdbout.exists()


def test_xmlout(build, parameters):
    _, output = parameters
    assert os.path.exists(str(output.xmlout))


def test_buildtime(build,  expected_values):
    assert build < expected_values["time"]


def test_xml(xml, expected_values):
    fragments_built = int(xml.find("Final/FragmentsBuilt").text)
    residues_built = int(xml.find("Final/ResiduesBuilt").text)
    residues_sequenced = int(xml.find("Final/ResiduesSequenced").text)
    longest_fragment = int(xml.find("Final/ResiduesLongestFragment").text)

    assert fragments_built >= expected_values["fragments_built"]
    assert residues_built >= expected_values["residues_built"]
    assert residues_sequenced >= expected_values["residues_sequenced"]
    assert longest_fragment >= expected_values["longest_fragment"]
