import subprocess
import site
import os
import hashlib
import pytest

def get_md5_sums():
    phosphate_md5_sum = "61c0f43b46ff9a686f23c63b4bac0583"
    sugar_md5_sum = "5e3f3a2024d55d8ad548560b00dbdbeb"
    base_md5_sum = "90d983747bc8b2efadb695c6cb802596"
    return {
        "phosphate": phosphate_md5_sum,
        "sugar": sugar_md5_sum,
        "base": base_md5_sum
    }


def perform_installation_test(model_type: str):
    command = f"nucleofind-install -m {model_type}"
    process = subprocess.run(command.split(" "))
    assert process.returncode == 0

    data = get_md5_sums()

    found_models = False
    found_sum = None
    for folder in site.getsitepackages():
        if 'nucleofind_models' in os.listdir(folder):
            possible_location = os.path.join(folder, "nucleofind_models", f"{model_type}.onnx")
            if os.path.exists(possible_location):
                found_models = True

                # Calculate file hash to ensure file was downloaded correctly
                with open(possible_location, "rb") as f:
                    file_hash = hashlib.md5()
                    while chunk := f.read(8192):
                        file_hash.update(chunk)

                    found_sum = file_hash.hexdigest()
                break

    assert found_models == True
    assert found_sum == data[model_type]


@pytest.mark.slow
def test_install():
    perform_installation_test("phosphate")
    perform_installation_test("sugar")
    perform_installation_test("base")
