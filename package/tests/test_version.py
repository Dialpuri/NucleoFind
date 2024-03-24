import argparse
import nucleofind as nf
import nucleofind.build as nfb
import nucleofind.install as nfi
import re

def main(args):

    pattern = r'v(\d+\.\d+\.\d+)'
    match = re.search(pattern, args.version)
    version = args.version
    if match:
        version = match.group(1)
    else:
        raise RuntimeError(f'Version {args.version} does not match pattern v0.0.0 - Name')
    print(version)
    assert nf.__version__ == version
    assert nfb.__version__ == version
    assert nfi.__version__ == version
    assert nfb.__version__ == nf.__version__

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-version', required=True)
    args = parser.parse_args()
    main(args)