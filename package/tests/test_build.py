import os

import nucleofind.build as b
import pytest

@pytest.fixture
def data():
    mtzin = "tests/test_data/5d5w/hklout.mtz"
    pdbin = "tests/test_data/5d5w/xyzout.pdb"
    seqin = "tests/test_data/5d5w/5d5w.fasta"
    pdbout = "tests/5d5w.pdb"
    predin = "tests/test_data/5d5w/phosphate.map"
    colinfo = "FP,SIGFP"
    colinfc = "FWT,PHWT"
    colinfree = "FREE"
    return {"mtzin": mtzin, "pdbin": pdbin,
            "seqin": seqin, "pdbout": pdbout,
            "predin": predin, "colinfo": colinfo,
            "colinfc": colinfc, "colinfree": colinfree}


# mtzin: str, seqin: str, pdbin: str, predin: str, colin_fo: str, colin_fc: str, colin_free: str, pdbout: str
def test_build(data):
    os.environ["CLIBD"]="tests/test_data"
    os.environ["LD_LIBRARY_PATH"]="tests/test_data"
    os.environ["CCP4"]="tests/test_data"

    b.build(data["mtzin"], data["seqin"], data["pdbin"], data["predin"], data["colinfo"], data["colinfc"], data["colinfree"], data["pdbout"])
    assert os.path.exists(data["pdbout"])
