from .nautilus_module import Input, Output, run
import argparse
from .__version__ import __version__
from .predict import predict_map
import traceback

def main():
    parser = argparse.ArgumentParser(description='nucleofind build')
    parser.add_argument("-mtzin", required=True)
    parser.add_argument("-seqin", required=True)
    parser.add_argument("-pdbin", required=False, default="")
    parser.add_argument("-pdbout", required=True)
    parser.add_argument("-phosin", required=True)
    parser.add_argument("-sugarin", required=False, default="")
    parser.add_argument("-basein", required=False, default="")
    parser.add_argument("-colin-fo", required=True)
    parser.add_argument("-colin-fc", required=True)
    parser.add_argument("-colin-free", required=False, default="")
    parser.add_argument("-xmlout", required=False, default="")
    parser.add_argument("-cycles", required=False, default=3)
    parser.add_argument("-v", "--version", action="version", version=__version__)

    args = parser.parse_args()

    if args.phosin == "auto":
        print("Before building, NucleoFind will predict a phosphate map and output it into the current working directory")
        split_calculated_sf = args.colin_fc.split(",")

        if not split_calculated_sf:
            print("No Calculated SFs provided, attempting to use FWT and PHWT")
            split_calculated_sf = ["FWT", "PHWT"]

        predict_map(model="phosphate",
                    input=args.mtzin,
                    output="phosphate.map",
                    intensity=split_calculated_sf[0],
                    phase=split_calculated_sf[1])

        args.predin = "phosphate.map"

    input = Input(args.mtzin,
                  args.seqin,
                  args.pdbin,
                  args.phosin,
                  args.sugarin,
                  args.basein,
                  args.colin_fo,
                  "",
                  "",
                  args.colin_fc,
                  args.colin_free)

    output = Output(args.pdbout, args.xmlout)

    try:
        run(input, output, args.cycles)
    except:
        print(traceback.format_exc())



def build(mtzin: str, seqin: str, pdbin: str, predin: str, colin_fo: str, colin_fc: str, colin_free: str, pdbout: str,
          xmlout: str, cycles: int):
    input = Input(mtzin, seqin, pdbin, predin, "", "", colin_fo, "", "", colin_fc,
                  colin_free)
    output = Output(pdbout, xmlout)
    run(input, output, cycles)
