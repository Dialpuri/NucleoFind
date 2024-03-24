from .nautilus_module import Input, Output, run
import argparse
from .__version__ import __version__


def main():
    parser = argparse.ArgumentParser(description='nucleofind build')
    parser.add_argument("-mtzin", required=True)
    parser.add_argument("-seqin", required=True)
    parser.add_argument("-pdbin", required=False, default="")
    parser.add_argument("-pdbout", required=True)
    parser.add_argument("-predin", required=True)
    parser.add_argument("-colin-fo", required=True)
    parser.add_argument("-colin-fc", required=True)
    parser.add_argument("-colin-free", required=False, default="")
    parser.add_argument("-cycles", required=False, default=3)
    parser.add_argument("-v", "--version", action="version", version=__version__)

    args = parser.parse_args()

    input = Input(args.mtzin, args.seqin, args.pdbin, args.predin, "", "", args.colin_fo, "", "", args.colin_fc,
                  args.colin_free)
    output = Output(args.pdbout)

    run(input, output, args.cycles)


def build(mtzin: str, seqin: str, pdbin: str, predin: str, colin_fo: str, colin_fc: str, colin_free: str, pdbout: str):
    input = Input(mtzin, seqin, pdbin, predin, "", "", colin_fo, "", "", colin_fc,
                  colin_free)
    output = Output(pdbout)
    run(input, output, 1)
