import logging
from pathlib import Path
from typing import List
import re


def show_missing_model_error():
    """Show error when no models are found"""
    logging.critical("""
    No models have been found in either site_packages or CCP4/lib/data.
    You can install models using the command:
    nucleofind-install -m {nano,core}
        """)


def show_missing_specified_model_error(model_name: str):
    """Show error when model with specified name is not found"""
    logging.critical(f"""
    No model with the name {model_name} has been found in either site_packages or CCP4/lib/data.""")


def show_multiple_model_error(model_names: List[str]):
    """Show warning when multiple models are found"""
    multiple_model_names = ""
    for model_name in model_names:
        multiple_model_names += f"\t-model {model_name}\n"

    logging.warning(f"""
    Multiple models have been found in either site_packages or CCP4/lib/data.
    Please specify either: 
    {multiple_model_names}""")
