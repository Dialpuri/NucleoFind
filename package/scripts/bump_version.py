import re
import argparse


def bump_version(version, args):
    major, minor, patch = map(int, version.split("."))
    if args.major:
        major += 1
        minor = 0
        patch = 0

    if args.minor:
        minor += 1

    if args.patch:
        patch += 1

    return f"{major}.{minor}.{patch}"


def update_version_in_file(file_path, old_version, new_version):
    with open(file_path, "r") as file:
        content = file.read()

    updated_content = content.replace(old_version, new_version)

    with open(file_path, "w") as file:
        file.write(updated_content)


def main():
    parser = argparse.ArgumentParser(description="Bump project version.")
    parser.add_argument("-major", action=argparse.BooleanOptionalAction)
    parser.add_argument("-minor", action=argparse.BooleanOptionalAction)
    parser.add_argument("-patch", action=argparse.BooleanOptionalAction)
    args = parser.parse_args()

    # Update __version__.py
    version_file = "src/nucleofind/__version__.py"
    version_regex = r"(?:__version__\s*=\s*[\"'])(\d+\.\d+\.\d+)([\"'])"
    with open(version_file, "r") as file:
        match = re.search(version_regex, file.read())
        if not match:
            raise ValueError("Current version not found in __version__.py")
        current_version = match.group(1)

    # Update writerside.cfg
    cfg_file = "../docs/Writerside/writerside.cfg"
    cfg_version_regex = r'(<instance [^>]*version=")(\d+\.\d+\.\d+)(")'

    with open(cfg_file, "r") as file:
        match = re.search(cfg_version_regex, file.read())
        if not match:
            raise ValueError("Current version not found in writerside.cfg")
        cfg_current_version = match.group(2)

    if current_version != cfg_current_version:
        raise RuntimeError(f"Docs and Python version mismatch, {current_version=} != {cfg_current_version=}")

    new_version = bump_version(current_version, args)

    print(f"Bumping version {current_version} -> {new_version}")
    update_version_in_file(version_file, current_version, new_version)
    update_version_in_file(cfg_file, current_version, new_version)


if __name__ == "__main__":
    main()
