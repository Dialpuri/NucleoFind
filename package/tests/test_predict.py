import nucleofind.predict as p
import numpy as np
import pytest

@pytest.fixture
def data():
    model_path = "../models/phosphate.onnx"
    mtzin = "tests/test_data/5d5w/hklout.mtz"
    pdbin = "tests/test_data/5d5w/xyzout.pdb"
    seqin = "tests/test_data/5d5w/5d5w.fasta"
    pdbout = "tests/test_output/5d5w.pdb"
    predin = "tests/test_data/5d5w/phosphate.map"
    colinfo = "FP,SIGFP"
    colinfc = "FWT,PHWT"
    colinfree = "FREE"
    prediction = p.Prediction(model_dir=model_path, use_cache=False)
    return {"model_path": model_path, "mtzin": mtzin, "pdbin": pdbin,
            "seqin": seqin, "pdbout": pdbout, "predin": predin, "colinfo": colinfo,
            "colinfc": colinfc, "colinfree": colinfree,"prediction": prediction}

def test_prediction(data):
    data["prediction"].make_prediction(data["mtzin"], data["colinfc"].split(","), overlap=32)
    (i0, i1), (p0, p1) = np.unique(data["prediction"].predicted_map, return_counts=True)
    assert p1 > 2_000
    assert p0 > 1_000_000
    assert i0 == 0.0
    assert i1 == 1.0

