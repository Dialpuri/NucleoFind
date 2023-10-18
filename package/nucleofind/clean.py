import site 
import os
import urllib.request
import argparse
import enum


class ModelType(enum.Enum):
    phosphate = 0
    sugar = 1
    base = 2

class InstallLocation(enum.Enum):
    site_packages = 0
    ccp4 = 1


def clibd_error_msg():
    print(   """CCP4 Environment Variable - CLIBD is not found.
                You can try sourcing it: 
                Ubuntu - source /opt/xtal/ccp4-X.X/bin/ccp4.setup-sh
                MacOS - source /Applications/ccp4-X.X/bin/ccp4.setup-sh
                """
                )

def download_model(type: ModelType, folder: str, reinstall: bool):
    model_name = f"{type.name}.hdf5"
    url = f"http://www.ysbl.york.ac.uk/jsd523/{model_name}"

    nucleofind_model_dir = os.path.join(folder, "nucleofind_models")
    if not os.path.isdir(nucleofind_model_dir):
        os.makedirs(nucleofind_model_dir)

    model_path = os.path.join(nucleofind_model_dir, model_name )

    if os.path.exists(model_path) and not reinstall:
        print(f"Found {type.name} model in {model_path}, skipping.")
        return

    urllib.request.urlretrieve(url, model_path)  

    if not os.path.exists(model_path):
        print("Something has gone wrong with the download. Aborting.")

def install_model(type: ModelType, location: str, reinstall: bool) -> bool: 
    print("Installing model of type", type.name, "to", location)
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

def clean_models(): 
    site_packages_dir = site.getsitepackages()

    found_models = []

    for folder in site_packages_dir:
        nucleofind_model_dir = os.path.join(folder, "nucleofind_models")
        for model in os.scandir(nucleofind_model_dir):
            found_models.append(model)

    if not found_models:
        print("No models were found in site-packages. Finishing.")
        return
    
    print("Pick an option to remove: ")
    if len(found_models)>1:
        print(f"1) All")

    for index, model in enumerate(found_models):
        print(f"{index+2}) {model.name}")


    option_selected = False
    while not option_selected:
        option = input("Number: ")

        try:
            choice = int(option)
            if choice <= 0 or choice > len(found_models)+1:
                raise ValueError()
            
            model_to_remove = found_models[choice-2]

            if choice == 1:
                print("Do you want to remove all the models?")
            else:
                print("You want to remove ")

            y_no_selected = False
            confirm = False
            while not y_no_selected:
                y_or_n = input("Y/N ").lower()
                if y_or_n not in ["y", "yes", "n", "no"]:
                    continue
                
                y_no_selected = True

                if y_or_n == "y" or y_or_n == "yes":
                    confirm = True

            if confirm:
                if choice == 1:
                    for model in found_models:
                        os.remove(model.path)
                        print("Removed", model.name)
                else:
                    os.remove(model_to_remove.path)
                    print("Removed", model_to_remove.name)
            

            option_selected = True
        except ValueError :
            print("Invalid choice")





def run():
    clean_models()
