#  Copyright (c) 2024 Jordan Dialpuri, Jon Agirre, Kathryn Cowtan, Paul Bond and University of York. All rights reserved

import dataclasses
import pathlib
import site
import traceback
from importlib.metadata import PackageNotFoundError


from .nautilus_module import Input, Output, run
import argparse
from .__version__ import __version__
from .prediction.predict import predict_map


@dataclasses.dataclass
class InputParameters:
    mtzin: str | pathlib.Path = ""
    pdbin: str | pathlib.Path = ""
    seqin: str | pathlib.Path = ""
    phosin: str | pathlib.Path = ""
    sugarin: str | pathlib.Path = ""
    basein: str | pathlib.Path = ""
    colinfo: str | pathlib.Path = ""
    colinphifom: str | pathlib.Path = ""
    colinfc: str | pathlib.Path = ""
    colinfree: str | pathlib.Path = ""
    em: bool = False
    cycles: int = 3


@dataclasses.dataclass
class OutputParameters:
    pdbout: str | pathlib.Path = ""
    xmlout: str | pathlib.Path = ""


def find_database() -> str:
    for pkg in site.getsitepackages():
        database_path = pathlib.Path(pkg) / "nucleofind_models" / "nucleofind-database.cif"
        if database_path.exists():
            return str(database_path)
    return ""

def main():
    parser = argparse.ArgumentParser(description="nucleofind build")
    parser.add_argument("--mtzin", required=True)
    parser.add_argument("--seqin", required=False, default="")
    parser.add_argument("--pdbin", required=False, default="")
    parser.add_argument("--pdbout", required=True)
    parser.add_argument("--phosin", required=False)
    parser.add_argument("--sugarin", required=False, default="")
    parser.add_argument("--basein", required=False, default="")
    parser.add_argument("--preddirin", required=False)
    parser.add_argument("--colin-fo", required=False, default="")
    parser.add_argument("--colin-fc", required=False, default="")
    parser.add_argument("--colin-phifom", required=False, default="")
    parser.add_argument("--colin-free", required=False, default="")
    parser.add_argument("--xmlout", required=False, default="")
    parser.add_argument("--cycles", required=False, default=3)
    parser.add_argument("--em", required=False, default=False, action='store_true')
    parser.add_argument("-v", "--version", action="version", version=__version__)

    args = parser.parse_args()

    if not args.phosin and not args.preddirin:
        parser.error("Please specify a predicted map, either -phosin, -sugarin, -basein or -preddirin")

    if args.phosin == "auto":
        print(
            "Before building, NucleoFind will predict a phosphate map and output it into the current working directory"
        )
        split_calculated_sf = []
        if args.colin_fc:
            split_calculated_sf = args.colin_fc.split(",")
        if args.colin_fo and args.colin_phifom:
            F = args.colin_fo.split(",")[0]
            PHI = args.colin_phifom.split(",")[0]
            split_calculated_sf = [F, PHI]

        if not split_calculated_sf:
            print("No Calculated SFs provided, attempting to use FWT and PHWT")
            split_calculated_sf = ["FWT", "PHWT"]

        predict_map(
            model="core",
            input=args.mtzin,
            output="nucleofind-run",
            amplitude=split_calculated_sf[0],
            phase=split_calculated_sf[1],
        )

        args.phosin = "nucleofind-run/nucleofind-phosphate.map"
        args.sugarin = "nucleofind-run/nucleofind-sugar.map"
        args.basein = "nucleofind-run/nucleofind-base.map"

    if args.preddirin:
        args.phosin = f"{args.preddirin}/nucleofind-phosphate.map"
        args.sugarin = f"{args.preddirin}/nucleofind-sugar.map"
        args.basein = f"{args.preddirin}/nucleofind-base.map"

    if not args.colin_fc and not args.colin_phifom:
        raise ValueError("Please specify columns for the phases, either FWT,PHWT or PHI,FOM")

    database = find_database()

    input = Input(
        args.mtzin,
        args.seqin,
        args.pdbin,
        args.phosin,
        args.sugarin,
        args.basein,
        args.colin_fo,
        "",
        args.colin_phifom,
        args.colin_fc,
        args.colin_free,
        args.em,
        database
    )

    output = Output(args.pdbout, args.xmlout)

    try:
        run(input, output, args.cycles)
    except:
        print(traceback.format_exc())


def build(input_parameters: InputParameters, output_parameters: OutputParameters):
    database = find_database()

    input = Input(
        str(input_parameters.mtzin),
        str(input_parameters.seqin),
        str(input_parameters.pdbin),
        str(input_parameters.phosin),
        str(input_parameters.sugarin),
        str(input_parameters.basein),
        str(input_parameters.colinfo),
        str(""),
        str(input_parameters.colinphifom),
        str(input_parameters.colinfc),
        str(input_parameters.colinfree),
        input_parameters.em,
        database
    )
    output = Output(str(output_parameters.pdbout), str(output_parameters.xmlout))

    if not all([input_parameters.phosin, input_parameters.sugarin, input_parameters.basein]):
        raise RuntimeError("The phosphate, sugar and base predictions are required.")

    try:
        run(input, output, input_parameters.cycles)
    except:
        print(traceback.format_exc())



# def build(mtzin: str, seqin: str, pdbin: str, phosin: str, sugarin: str, basein: str, colin_fo: str, colin_fc: str,
#           colin_free: str, pdbout: str,
#           xmlout: str, cycles: int):
#     input = Input(mtzin, seqin, pdbin, phosin, sugarin, basein, colin_fo, "", "", colin_fc,
#                   colin_free)
#     output = Output(pdbout, xmlout)
#     run(input, output, cycles)
