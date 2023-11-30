import site
import os


def clean_models():
    site_packages_dir = site.getsitepackages()

    found_models = []

    for folder in site_packages_dir:
        nucleofind_model_dir = os.path.join(folder, "nucleofind_models")
        if os.path.exists(nucleofind_model_dir):
            for model in os.scandir(nucleofind_model_dir):
                found_models.append(model)

    if not found_models:
        print("No models were found in site-packages. Finishing.")
        return

    print("Pick an option to remove: ")
    if len(found_models) > 1:
        print(f"1) All")

    for index, model in enumerate(found_models):
        print(f"{index + 2}) {model.name}")

    option_selected = False
    while not option_selected:
        option = input("Number: ")

        try:
            choice = int(option)
            if choice <= 0 or choice > len(found_models) + 1:
                raise ValueError()

            model_to_remove = found_models[choice - 2]

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
        except ValueError:
            print("Invalid choice")


def run():
    clean_models()
