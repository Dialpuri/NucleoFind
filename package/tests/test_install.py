import subprocess
import site
import os
import hashlib
from pathlib import Path

import pytest

def get_md5_sums():
    core_md5sum = "86a612ef06930355d7ac87b30e97f779"
    nano_md5sum = "bd88177bd06d593675d9eb5498c0048f"
    return {
        "core": core_md5sum,
        "nano": nano_md5sum,
    }


def perform_installation_test(model_type: str):
    command = f"nucleofind-install -m {model_type}"
    process = subprocess.run(command.split(" "))
    # assert process.returncode == 0

    data = get_md5_sums()

    found_models = False
    found_sum = None
    for folder in site.getsitepackages():
        folder = Path(folder)
        model_directory = list(folder.glob("nucleofind_models"))
        if model_directory:
            possible_location = folder / "nucleofind_models" / f"nucleofind-{model_type}.onnx"
            if possible_location.exists():
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
def test_core_install():
    perform_installation_test("core")


@pytest.mark.slow
def test_nano_install():
    perform_installation_test("nano")

@pytest.mark.slow
def test_ultra_install():
    perform_installation_test("ultra")


