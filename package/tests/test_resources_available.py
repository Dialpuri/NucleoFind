import os
import hashlib
from pathlib import Path

import pytest

def get_md5_sums():
    core_md5sum = "a843baa5eb38e2ced9298177771d44d3"
    database_md5sum = "5e1e8b980b2359e74e1c9151780af926"

    return {
        "core": core_md5sum,
        "database": database_md5sum,
    }



def calculate_md5_hash(file_path: str) -> str:
    with open(file_path, "rb") as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)

        found_sum = file_hash.hexdigest()
    return found_sum

@pytest.mark.ccp4
def test_files_intact():
    data = get_md5_sums()

    clibd = os.environ.get("CLIBD", "")
    assert clibd != "", "Could not find CLIBD environment variable"

    nucleofind_model_folder_path = Path(clibd) / "nucleofind_models"
    assert nucleofind_model_folder_path.exists(), "Could not find nucleofind_models folder in CLIBD"

    onnx_model_path = nucleofind_model_folder_path / f"nucleofind-core.onnx"
    assert onnx_model_path.exists(), f"Could not find nucleofind-core.onnx model in nucleofind_models folder"

    onnx_model_sum = calculate_md5_hash(str(onnx_model_path))
    assert onnx_model_sum == data["core"], f"Model checksum mismatch for {onnx_model_path}"

    database_path = nucleofind_model_folder_path / "nucleofind-database.cif"
    assert database_path.exists(), f"Could not find nucleofind-database.cif in nucleofind_models folder"

    database_sum = calculate_md5_hash(str(database_path))
    assert database_sum == data["database"], f"Database checksum mismatch for {database_path}"

@pytest.mark.ccp4
def test_files_present():
    clibd = os.environ.get("CLIBD", "")
    assert clibd != "", "Could not find CLIBD environment variable"

    nucleofind_model_folder_path = Path(clibd) / "nucleofind_models"
    assert nucleofind_model_folder_path.exists(), "Could not find nucleofind_models folder in CLIBD"

    onnx_model_path = nucleofind_model_folder_path / f"nucleofind-core.onnx"
    assert onnx_model_path.exists(), f"Could not find nucleofind-core.onnx model in nucleofind_models folder"



