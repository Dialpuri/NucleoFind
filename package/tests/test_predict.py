import nucleofind.predict as p
import numpy as np
import pytest
import hashlib

@pytest.fixture
def data():
    model_path = p.find_model("phosphate")
    mtzin = "tests/test_data/5d5w/hklout.mtz"
    pdbin = "tests/test_data/5d5w/xyzout.pdb"
    seqin = "tests/test_data/5d5w/5d5w.fasta"
    pdbout = "tests/test_output/5d5w.pdb"
    predin = "tests/test_data/5d5w/phosphate.map"
    colinfo = "FP,SIGFP"
    colinfc = "FWT,PHWT"
    colinfree = "FREE"
    phosphate_map_md5sum = "f9f3470859cac6c976b5146a635c58e7"
    prediction = p.Prediction(model_dir=model_path)
    return {"model_path": model_path, "mtzin": mtzin, "pdbin": pdbin,
            "seqin": seqin, "pdbout": pdbout, "predin": predin, "colinfo": colinfo,
            "colinfc": colinfc, "colinfree": colinfree, "prediction": prediction,
            "phosphate_map_md5sum": phosphate_map_md5sum}

def test_prediction(data):
    data["prediction"].make_prediction(data["mtzin"], data["colinfc"].split(","), overlap=32)
    (i0, i1), (p0, p1) = np.unique(data["prediction"].predicted_map, return_counts=True)
    assert p1 > 2_000
    assert p0 > 1_000_000
    assert i0 == 0.0
    assert i1 == 1.0

    data["prediction"].save_predicted_map("phosphate.map")

    with open("phosphate.map", "rb") as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)

        found_sum = file_hash.hexdigest()
    assert found_sum == data["phosphate_map_md5sum"]
