import argparse
from types import SimpleNamespace

from nucleofind.__version__ import __version__
from .model import ModelType


def parse_arguments() -> SimpleNamespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "-model",
        help="Model selection",
        choices=[type.name for type in ModelType],
        required=False,
    )
    parser.add_argument("-i", "-input", help="Input mtz", required=True)
    parser.add_argument(
        "-o",
        "-output",
        help="Output directory, if does not exist it will be created model",
        default="nucleofind-output",
        required=False,
    )
    parser.add_argument("-r", "-resolution", nargs="?", help="Resolution cutoff")
    parser.add_argument(
        "-n", "-nthreads", nargs="?", default=None, type=int, help="Number of threads to use"
    )
    parser.add_argument(
        "-amplitude", "-f", nargs="?", help="Name of amplitude column in MTZ, e.g. FWT"
    )
    parser.add_argument(
        "-phase", "-phi", nargs="?", help="Name of phase column in MTZ, e.g. PHWT"
    )
    parser.add_argument(
        "-overlap",
        nargs="?",
        help="Amount of overlap to use",
        default=None,
        type=int,
    )
    parser.add_argument(
        "-no-symmetry",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Compute predictions for the entire unit cell",
    )
    parser.add_argument(
        "-variance", action=argparse.BooleanOptionalAction, help="Output variance map"
    )
    parser.add_argument(
        "-raw", action=argparse.BooleanOptionalAction, help="Output raw map (no argmax)"
    )
    parser.add_argument(
        "-gpu", action=argparse.BooleanOptionalAction, help="Use GPU (experimental)"
    )
    parser.add_argument(
        "-debug", action=argparse.BooleanOptionalAction, help="Turn on debug logging"
    )
    parser.add_argument(
        "-silent",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Turn off progress bar",
    )
    parser.add_argument("-model_path", nargs="?", help="Path to model (development)")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    args = vars(parser.parse_args())
    return SimpleNamespace(**args)
