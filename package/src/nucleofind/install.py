import site
import os
import urllib.request
import argparse
import enum
from .__version__ import __version__


class ModelType(enum.Enum):
    phosphate = 0
    sugar = 1
    base = 2

class InstallLocation(enum.Enum):
    site_packages = 0
    ccp4 = 1


def clibd_error_msg():
    print("""CCP4 Environment Variable - CLIBD is not found.
                You can try sourcing it: 
                Ubuntu - source /opt/xtal/ccp4-X.X/bin/ccp4.setup-sh
                MacOS - source /Applications/ccp4-X.X/bin/ccp4.setup-sh
                """
          )


def download_model(type: ModelType, folder: str, reinstall: bool):
    model_name = f"{type.name}.onnx"
    url = f"http://www.ysbl.york.ac.uk/jsd523/{model_name}"

    nucleofind_model_dir = os.path.join(folder, "nucleofind_models")
    if not os.path.isdir(nucleofind_model_dir):
        os.makedirs(nucleofind_model_dir)

    model_path = os.path.join(nucleofind_model_dir, model_name)

    if os.path.exists(model_path) and not reinstall:
        print(f"Found {type.name} model in {model_path}, skipping.")
        return

    urllib.request.urlretrieve(url, model_path)

    if not os.path.exists(model_path):
        print("Something has gone wrong with the download. Aborting.")


def install_model(type: ModelType, location: str, reinstall: bool) -> bool:
    print(f"Installing {type.name} model to {location}")
    if InstallLocation[location] == InstallLocation.ccp4:
        clibd = os.environ.get("CLIBD", "")
        if not os.path.exists(clibd):
            clibd_error_msg()
            return False

        download_model(type=type, folder=clibd, reinstall=reinstall)
        return True

    if InstallLocation[location] == InstallLocation.site_packages:
        site_packages_dir = site.getsitepackages()
        if len(site_packages_dir) == 0:
            print("Site packages could not be found, ensure you are in a python virtualenvironment. Aborting.")
            return False

        download_model(type=type, folder=site_packages_dir[0], reinstall=reinstall)
        return True


def run():
    model_list = ["phosphate", "sugar", "base"]
    output_list = ["ccp4", "site_packages"]

    parser = argparse.ArgumentParser(description='nucleofind Install')
    parser.add_argument("-m", "--model", choices=model_list, required=False)
    parser.add_argument("-o", "--output", choices=[location.name for location in InstallLocation], required=False,
                        default=output_list[1])
    parser.add_argument("--all", required=False, action='store_true')
    parser.add_argument("--reinstall", required=False, action='store_true')
    parser.add_argument("-v", "--version", action="version", version=__version__)

    args = parser.parse_args()

    if not args.all and not args.model:
        print("Please specify a model you wish to download")
        return

    if args.all:
        status = [install_model(type=type, location=args.output, reinstall=args.reinstall) for type in ModelType]
        if not any(status):
            print("There was a problem with installation of one of the models.", status)
        return

    install_model(type=ModelType[args.model], location=args.output, reinstall=args.reinstall)
