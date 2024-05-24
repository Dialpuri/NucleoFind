import os
import nucleofind.build as nucleofind_build
import pytest
import xml.etree.ElementTree as ET
import time
import tempfile
from pathlib import Path

@pytest.fixture(scope='session')
def data_base_path():
    return Path(__file__).parent / "test_data"

@pytest.fixture(scope='session', params=[
    {"map_combination": [1, 0, 0]},
    {"map_combination": [1, 1, 0]},
    {"map_combination": [1, 0, 1]}
])
def permutation(request):
    return request.param['map_combination']

@pytest.fixture(scope='session')
def parameters(data_base_path, permutation):
    input = nucleofind_build.InputParameters()
    data_dir = data_base_path / "5d5w"
    input.mtzin = data_dir / "hklout.mtz"
    input.pdbin = data_dir / "xyzout.pdb"
    input.seqin = data_dir / "5d5w.fasta"
    input.phosin = data_dir / "phosphate.map" if permutation[0] else ""
    input.sugarin = data_dir / "sugar.map" if permutation[1] else ""
    input.basein = data_dir / "base.map" if permutation[2] else ""
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


def test_buildtime(build):
    assert build < 10

def test_xml(xml):
    fragments_built = int(xml.find("Final/FragmentsBuilt").text)
    residues_built = int(xml.find("Final/ResiduesBuilt").text)
    residues_sequenced = int(xml.find("Final/ResiduesSequenced").text)
    longest_fragment = int(xml.find("Final/ResiduesLongestFragment").text)

    assert fragments_built >= 1
    assert residues_built >= 5
    assert residues_sequenced >= 4
    assert longest_fragment >= 5
