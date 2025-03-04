import os
import site
import sys
import urllib
from pathlib import Path
import logging
from types import SimpleNamespace
from typing import Tuple, List

import requests
from enum import Enum
import re
import hashlib

from .errors import (
    show_missing_model_error,
    show_multiple_model_error,
    show_missing_specified_model_error,
)


class ModelType(Enum):
    """Types of NucleoFind Model available"""
    nano = 1
    core = 2



def calculate_sha256(file_path: Path):
    """Calculate SHA256 hash of file"""
    logging.debug("Calculating SHA256 hash for %s", file_path)
    with open(file_path, "rb") as f:
        file_hash = hashlib.sha256()
        while chunk := f.read(4096):
            file_hash.update(chunk)
        return file_hash.hexdigest()


def calculate_size(file_path: Path):
    """Calculate size of file"""
    logging.debug("Calculating size for %s", file_path)
    return file_path.stat().st_size


def get_latest_model_metadata(type: ModelType, latest_model: str) -> Tuple[str, str]:
    """Get latest model metadata from HuggingFace"""
    url = f"https://huggingface.co/dialpuri/NucleoFind-{type.name}/raw/main/{latest_model}"
    response = requests.get(url)
    text = response.text
    sha_match = re.search(r"sha256:([a-f0-9]+)", text)
    if not sha_match:
        raise RuntimeError("Failed to get SHA256 hash from model metadata")
    sha256 = sha_match.group(1)

    size_match = re.search(r"size ([0-9]+)", text)
    if not size_match:
        raise RuntimeError("Failed to get size from model metadata")
    size = size_match.group(1)
    return sha256, size


def is_model_valid(type: ModelType, model_path: Path, latest_model: str):
    """Compare current model hash with latest model hash"""
    current_model_hash = calculate_sha256(model_path)
    current_model_size = calculate_size(model_path)

    latest_model_hash, latest_model_size = get_latest_model_metadata(type, latest_model)
    if latest_model_hash != current_model_hash:
        logging.info("Latest model and current modal checksum failed")
        return False
    return True


def get_latest_model(type: ModelType) -> str:
    """Query the HuggingFace API to get URL for latest model"""
    base_url = "https://huggingface.co/api/models/Dialpuri/NucleoFind"
    url = f"{base_url}-{type.name}"
    logging.debug("Getting latest model for %s from %s", type.name, url)
    response = requests.get(url)
    json = response.json()
    if not json:
        raise RuntimeError("Failed to get model URL")

    siblings = json.get("siblings", None)
    if not siblings:
        raise RuntimeError("Failed to get siblings from model")

    possible_models = []
    for filename in siblings:
        file = filename.get("rfilename", "")
        if file.endswith(".onnx"):
            possible_models.append(file)

    # Get latest model out of list based on date
    possible_models = sorted(possible_models, reverse=True)
    latest_model = possible_models[0]
    logging.debug("Latest model for %s is %s", type.name, latest_model)
    return latest_model


def download_model(
    type: ModelType, folder: Path, reinstall: bool = False, dry_run: bool = False
):
    """Download model from HuggingFace"""
    latest_model = get_latest_model(type)
    nucleofind_model_dir = folder / "nucleofind_models"
    nucleofind_model_dir.mkdir(exist_ok=True)
    model_path = nucleofind_model_dir / f"nucleofind-{type.name}.onnx"

    # Check if model already exists and is the latest version.
    if model_path.exists() and not reinstall:
        status = is_model_valid(type, model_path, latest_model)
        if not status:
            logging.warning(
                "A model was found but did not pass the latest validation checks, it may be corrupted, or a newer version available. "
                "To update the model, run `nucleofind-install --update`"
            )
        else:
            logging.warning(
                "Model already exists at %s, skipping download.", model_path
            )
        return

    url = f"https://huggingface.co/dialpuri/NucleoFind-{type.name}/resolve/main/{latest_model}?download=true"
    logging.debug("Downloading model from %s", url)
    if not dry_run:
        urllib.request.urlretrieve(url, model_path)

        if not is_model_valid(type, model_path, latest_model):
            logging.error("Model verification failed, model may be corrupted.")


def find_all_potential_models():
    """Find all potential models in site-packages and CCP4"""
    model_extension = "*.onnx"

    potential_models = []

    for pkg in site.getsitepackages():
        model_directory = Path(pkg) / "nucleofind_models"
        models = list(model_directory.glob(model_extension))
        potential_models += models

    clibd = Path(os.environ.get("CLIBD", ""))
    if not clibd.exists() and not potential_models:
        logging.warning(
            """CCP4 Environment Variable - CLIBD is not found. 
            You can try sourcing it: 
            Ubuntu - source /opt/xtal/ccp4-X.X/bin/ccp4.setup-sh
            MacOS - source /Applications/ccp4-X.X/bin/ccp4.setup-sh
            """
        )
        return

    ccp4_model_path = clibd / "nucleofind_models"
    if not ccp4_model_path.exists() and not potential_models:
        show_missing_model_error()
        return

    potential_models += list(ccp4_model_path.glob(model_extension))

    if not potential_models:
        show_missing_model_error()
        sys.exit(1)

    return [Path(x) for x in potential_models]


def extract_model_names(models: List[Path]) -> List[str]:
    """Extract model names from model paths"""
    model_names = []
    for model in models:
        match = re.search(r"nucleofind-(\w+).onnx", model.name)
        if not match:
            raise RuntimeError(
                "Failed to extract model name from model path, have the models been renamed? Please report this issue on GitHub."
            )
        name = match.group(1)
        model_names.append(name)
    return model_names


def find_model(model: ModelType | str | None) -> Path | None:
    """Search through site-packages and CCP4/lib/data for a potential model"""
    potential_models = find_all_potential_models()
    if not potential_models:
        sys.exit(1)

    if not model and len(potential_models) == 1:
        return Path(potential_models[0])

    model_names = extract_model_names(potential_models)
    if not model:
        show_multiple_model_error(model_names)
        sys.exit(1)

    if isinstance(model, ModelType):
        specified_model_name = model.name
    else:
        specified_model_name = model

    for name in model_names:
        if name == specified_model_name:
            return Path(potential_models[model_names.index(name)])

    show_missing_specified_model_error(specified_model_name)
    sys.exit(1)


def get_model_config(model_path: Path, overlap: int | None) -> SimpleNamespace:
    """Get model configuration from model type"""
    model_type = model_path.stem.removeprefix("nucleofind-")
    if model_type not in ModelType.__members__:
        raise RuntimeError(f"Invalid model type - {model_type}")
    model_type = ModelType[model_type]
    match model_type:
        case ModelType.nano:
            return SimpleNamespace(box_size=128, overlap=64 if overlap is None else overlap)
        case ModelType.core:
            return SimpleNamespace(box_size=128, overlap=64 if overlap is None else overlap)
        case _:
            raise RuntimeError(f"Invalid model type - {model_type}")
