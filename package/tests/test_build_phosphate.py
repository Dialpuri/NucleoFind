import os
import nucleofind.build as b
import pytest
import xml.etree.ElementTree as ET
import time

@pytest.fixture(scope='session')
def build():
    # Setup test environment
    mtzin = "tests/test_data/5d5w/hklout.mtz"
    pdbin = "tests/test_data/5d5w/xyzout.pdb"
    seqin = "tests/test_data/5d5w/5d5w.fasta"
    pdbout = "tests/5d5w.pdb"
    phosin = "tests/test_data/5d5w/phosphate.map"
    colinfo = "FP,SIGFP"
    colinfc = "FWT,PHWT"
    colinfree = "FREE"
    xmlout = "tests/5d5w.xml"
    cycles = 3
    data = {"mtzin": mtzin, "pdbin": pdbin,
            "seqin": seqin, "pdbout": pdbout,
            "phosin": phosin, "colinfo": colinfo,
            "colinfc": colinfc, "colinfree": colinfree,
            "xmlout": xmlout, "cycles": cycles}
    os.environ["CLIBD"]="tests/test_data"
    os.environ["LD_LIBRARY_PATH"]="tests/test_data"
    os.environ["CCP4"]="tests/test_data"

    if os.path.exists(data["xmlout"]):
        os.remove(data["xmlout"])

    if os.path.exists(data["pdbout"]):
        os.remove(data["pdbout"])

    t0 = time.time()
    b.build(data["mtzin"], data["seqin"], data["pdbin"], data["phosin"], "", "", data["colinfo"], data["colinfc"],
            data["colinfree"], data["pdbout"], data["xmlout"], data["cycles"])
    t1 = time.time()
    data["time"] = (t1-t0)
    tree = ET.parse(data["xmlout"])
    root = tree.getroot()
    data["xmlroot"] = root
    return data

def test_pdbout(build):
    assert os.path.exists(build["xmlout"])

def test_xmlout(build):
    assert os.path.exists(build["xmlout"])

def test_buildtime(build):
    assert build["time"] < 10

def test_xml(build):

    fragments_built=int(build["xmlroot"].find("Final/FragmentsBuilt").text)
    residues_built=int(build["xmlroot"].find("Final/ResiduesBuilt").text)
    residues_sequenced=int(build["xmlroot"].find("Final/ResiduesSequenced").text)
    longest_fragment=int(build["xmlroot"].find("Final/ResiduesLongestFragment").text)

    assert fragments_built >= 1
    assert residues_built >= 5
    assert residues_sequenced >= 4
    assert longest_fragment >= 5

