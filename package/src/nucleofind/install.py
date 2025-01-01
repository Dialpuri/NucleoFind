#  Copyright (c) 2024 Jordan Dialpuri, Jon Agirre, Kathryn Cowtan, Paul Bond and University of York. All rights reserved

import site
import os
import urllib.request
import argparse
import enum
from pathlib import Path
from .__version__ import __version__
from .prediction import model
from .logs import setup_logging
import logging


class InstallLocation(enum.Enum):
    site_packages = 0
    ccp4 = 1


def clibd_error_msg():
    print("""CCP4 Environment Variable - CLIBD is not found.
                You can try sourcing it: 
                Ubuntu - source /opt/xtal/ccp4-X.X/bin/ccp4.setup-sh
                MacOS - source /Applications/ccp4-X.X/bin/ccp4.setup-sh
                """)


def install_model(type: model.ModelType, location: str, reinstall: bool) -> bool:
    logging.info(f"Installing {type.name} model to {location}")
    if InstallLocation[location] == InstallLocation.ccp4:
        clibd = os.environ.get("CLIBD", "")
        if not os.path.exists(clibd):
            clibd_error_msg()
            return False

        model.download_model(type=type, folder=clibd, reinstall=reinstall)
        return True

    if InstallLocation[location] == InstallLocation.site_packages:
        site_packages_dir = site.getsitepackages()
        if not site_packages_dir:
            raise RuntimeError(
                "Failed to get site packages directory, ensure you in a virtual environment"
            )
        first_sitepackages = Path(site_packages_dir[0])
        model.download_model(type=type, folder=first_sitepackages, reinstall=reinstall)
        return True


def run():
    setup_logging()
    output_list = ["ccp4", "site_packages"]

    parser = argparse.ArgumentParser(description="nucleofind Install")
    parser.add_argument(
        "-m", "--model", choices=[type.name for type in model.ModelType], required=False
    )
    parser.add_argument(
        "-o",
        "--output",
        choices=[location.name for location in InstallLocation],
        required=False,
        default=output_list[1],
    )
    parser.add_argument("--all", required=False, action="store_true")
    parser.add_argument("--update", required=False, action="store_true")
    parser.add_argument("-v", "--version", action="version", version=__version__)

    args = parser.parse_args()

    if not args.all and not args.model:
        print("Please specify a model you wish to download")
        return

    if args.all:
        status = [
            install_model(type=type, location=args.output, reinstall=args.update)
            for type in model.ModelType
        ]
        if not any(status):
            print("There was a problem with installation of one of the models.", status)
        return

    install_model(
        type=model.ModelType[args.model], location=args.output, reinstall=args.update
    )
