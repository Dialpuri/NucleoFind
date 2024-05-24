import nucleofind.predict as p
import numpy as np
import pytest
import hashlib
import os

@pytest.fixture
def data():
    model_paths = [p.find_model(m) for m in ["phosphate", "sugar", "base"]]
    mtzin = "tests/test_data/5d5w/hklout.mtz"
    pdbin = "tests/test_data/5d5w/xyzout.pdb"
    seqin = "tests/test_data/5d5w/5d5w.fasta"
    pdbout = "tests/test_output/5d5w.pdb"
    predin = "tests/test_data/5d5w/base.map"
    colinfo = "FP,SIGFP"
    colinfc = "FWT,PHWT"
    colinfree = "FREE"

    phosphate_map_md5sum = "f9f3470859cac6c976b5146a635c58e7"
    sugar_map_md5sum = "ffa9190ab0f9d4c7316fa101a7ab3e10"
    base_map_md5sum = "3fa375bb2435def19020088dfd45b490"

    prediction = p.Prediction(model_paths=model_paths)
    return {"model_paths": model_paths, "mtzin": mtzin, "pdbin": pdbin,
            "seqin": seqin, "pdbout": pdbout, "predin": predin, "colinfo": colinfo,
            "colinfc": colinfc, "colinfree": colinfree, "prediction": prediction,
            "phosphate_map_md5sum": phosphate_map_md5sum, "sugar_map_md5sum": sugar_map_md5sum,
            "base_map_md5sum": base_map_md5sum}

def test_prediction(data):
    data["prediction"].make_prediction(data["mtzin"], data["colinfc"].split(","), overlap=32)
    data["prediction"].save_predicted_map("prediction")

    assert os.path.exists("prediction_phosphate.map")
    assert os.path.exists("prediction_sugar.map")
    assert os.path.exists("prediction_base.map")

    with open("prediction_phosphate.map", "rb") as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)

        found_sum = file_hash.hexdigest()
    assert found_sum == data["phosphate_map_md5sum"]

    with open("prediction_sugar.map", "rb") as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)

        found_sum = file_hash.hexdigest()
    assert found_sum == data["sugar_map_md5sum"]

    with open("prediction_base.map", "rb") as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)

        found_sum = file_hash.hexdigest()
    assert found_sum == data["base_map_md5sum"]

    os.remove("prediction_phosphate.map")
    os.remove("prediction_sugar.map")
    os.remove("prediction_base.map")