# import vtk
import pandas as pd

from paraview.simple import (
    ExtractBlock,
    ExtractSurface,
    GetParaViewVersion,
    Delete,
    SaveData,
)


from .method import load, info, getbounds, resultinfo
from .view import deformed


pd.options.mode.copy_on_write = True

import argparse, argcomplete
import os
import sys

from pint import UnitRegistry, Quantity

# Ignore warning for pint
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Quantity([])


def options(description: str, epilog: str):
    """
    define options
    """
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    subparsers = parser.add_subparsers(
        title="dimmension", dest="dimmension", help="sub-dimmension help"
    )
    parser_3D = subparsers.add_parser("3D", help="3D model")
    parser_2D = subparsers.add_parser("2D", help="2D model")
    parser_Axi = subparsers.add_parser("Axi", help="Axi model")
    for allparsers in [parser_3D, parser_2D, parser_Axi]:
        allparsers.add_argument(
            "file", type=str, help="input case file (ex. Export.case)"
        )
        allparsers.add_argument(
            "--field", type=str, help="select field to display", default=""
        )
        allparsers.add_argument(
            "--stats", help="activate stats calculations", action="store_true"
        )
        allparsers.add_argument(
            "--histos", help="activate histograms calculations", action="store_true"
        )
        allparsers.add_argument(
            "--bins", type=int, help="set bins number (default 10)", default=10
        )
        allparsers.add_argument(
            "--plots", help="activate plots calculations", action="store_true"
        )
        allparsers.add_argument(
            "--views", help="activate views calculations", action="store_true"
        )
        allparsers.add_argument(
            "--r", nargs="*", type=float, help="select r in m to display"
        )
        if allparsers != parser_Axi:
            allparsers.add_argument(
                "--theta", nargs="*", type=float, help="select theta in deg to display"
            )
        if allparsers != parser_2D:
            allparsers.add_argument(
                "--channels", help="activate views calculations", action="store_true"
            )
            allparsers.add_argument(
                "--z", nargs="*", type=float, help="select z in m to display"
            )

        allparsers.add_argument(
            "--save",
            help="save graphs",
            action="store_true",
        )
        allparsers.add_argument(
            "--verbose",
            help="activate verbose mode",
            action="store_true",
        )

    # TODO get Exports section from json model file
    # data['PostProcess'][method_params[0]]['Exports']['expr']?
    #
    # provides
    # * field: symbol, unit, support, ...
    # * list of operations to perform (to be implemented?)
    #

    return parser


def main():

    parser = options("", "")
    argcomplete.autocomplete(parser)
    args = parser.parse_args()
    print(f"args: {args}")

    # get current working directory
    cwd = os.getcwd()

    basedir = f"{os.path.dirname(args.file)}/paraview.exports"
    # basedir = os.path.dirname(args.file).replace(f"{toolbox}.export", "paraview.export")
    print("Results are stored in: ", basedir)
    os.makedirs(basedir, exist_ok=True)

    # Pint configuration
    ureg = UnitRegistry()
    ureg.define("percent = 0.01 = %")
    ureg.define("ppm = 1e-6")
    ureg.default_system = "SI"
    ureg.autoconvert_offset_to_baseunit = True

    # set default output unit to millimeter
    distance_unit = "millimeter"  # or "meter"

    # fieldunits: dict( Quantity: symbol: str, units: [ in_unit, out_unit ]
    # !!! watchout matplotlib setup for latex support and packages required !!!

    # build fieldunits when creating setup and store as a json
    # load json
    match args.dimmension:
        case "3D":
            from .meshinfo import meshinfo
            from .case3D.plot import plotTheta, plotOr, plotOz
            from .case3D.display3D import makeview
            from .case3D.method3D import create_dicts

            dim = 3
            axis = False
        case "2D":
            from .meshinfo import meshinfo
            from .case2D.plot import plotTheta, plotOr
            from .case2D.display2D import makeview
            from .case2D.method2D import create_dicts

            dim = 2
            axis = False
        case "Axi":
            from .meshinfoAxi import meshinfo
            from .caseAxi.plot import plotOz, plotOr
            from .case2D.display2D import makeview
            from .caseAxi.methodAxi import create_dicts

            dim = 2
            axis = True
        case _:
            pass

    fieldunits, ignored_keys = create_dicts(ureg, distance_unit, basedir)

    # check paraview version
    version = GetParaViewVersion()
    print(f"Paraview version: {version}")

    # args.file = "../../HL-31/test/hybride-Bh27.7T-Bb9.15T-Bs9.05T_HPfixed_BPfree/bmap/np_32/thermo-electric.exports/Export.case"
    reader = load(args.file)
    # print(f"help(reader) = {dir(reader)}")
    bounds = getbounds(reader)
    print(f"bounds={bounds}")  # , type={type(bounds)}")
    info(reader)
    color = None
    CellToData = False
    if args.field:
        if args.field in list(reader.CellData.keys()):
            field = reader.CellData[args.field]
            # if field.GetNumberOfComponents() == 1:
            color = ["CELLS", args.field]
            # if field.GetNumberOfComponents() == dim:
            #    color = ["CELLS", args.field, "Magnitude"]
        if args.field in list(reader.PointData.keys()):
            field = reader.PointData[args.field]
            # if field.GetNumberOfComponents() == 1:
            color = ["POINTS", args.field]
            # if field.GetNumberOfComponents() == dim:
            #    color = ["POINTS", args.field, "Magnitude"]

    # get Block info
    cellsize, blockdata, statsdict = meshinfo(
        reader,
        dim,
        fieldunits,
        ignored_keys,
        basedir,
        ureg,
        ComputeStats=args.stats,
        ComputeHisto=args.histos,
        BinCount=args.bins,
        show=(not args.save),
        verbose=args.verbose,
    )

    # Plots
    if args.plots:
        if args.r:
            if axis:
                print(f"plots: r={args.r}, z={args.z}")
                # plotOr(reader, r, z, show=(not args.save))# with r=[r1, r2], z: float
                # plotOz(reader, r, z, show=(not args.save)) # with r: float, z=[z1,z2]
            elif dim == 3 and args.z:
                for z in args.z:
                    for r in args.r:
                        plotTheta(cellsize, r, z, show=(not args.save))
            elif dim == 2:
                for r in args.r:
                    plotTheta(cellsize, r, basedir, show=(not args.save))
                    for theta in args.theta:
                        plotOr(cellsize, r, theta, basedir, show=(not args.save))
        # add plotOr
        # add plotOz

    # When dealing with elasticity
    suffix = ""
    cellsize_deformed = None
    datadict = resultinfo(cellsize, ignored_keys)
    found = False
    for field in list(reader.PointData.keys()):
        if field.endswith("displacement"):
            found = True
            break
    print(f"displacement found={found} in {list(reader.PointData.keys())}")

    if found and (dim == 3 or axis):
        # make3Dview(cellsize, blockdata, key, color, addruler=True)
        if args.channels:
            print("Save stl for original geometries:")
            for i, block in enumerate(blockdata.keys()):
                name = blockdata[block]["name"]
                actual_name = name.replace("/root/", "")
                print(f"\t{name}: actual_name={actual_name}", end="")
                if not actual_name.endswith("Isolant") and not "Air" in actual_name:
                    print(" saved", end="", flush=True)
                    extractBlock1 = ExtractBlock(registrationName=name, Input=cellsize)
                    extractBlock1.Selectors = [block]
                    extractBlock1.UpdatePipeline()
                    extractSurface1 = ExtractSurface(
                        registrationName="ExtractSurface1", Input=extractBlock1
                    )

                    print(f" file={basedir}/{actual_name}.stl", flush=True)
                    SaveData(f"{basedir}/{actual_name}.stl", proxy=extractSurface1)
                    Delete(extractBlock1)
                    del extractBlock1
                else:
                    print(" ignored", flush=True)

        geometry = deformed(cellsize, factor=1)

        # compute channel deformation
        # use MeshLib see test-meshlib example
        if args.channels:
            print("Save stl for deformed geometries:")
            for i, block in enumerate(blockdata.keys()):
                name = blockdata[block]["name"]
                actual_name = name.replace("/root/", "")
                print(f"\t{name}: actual_name={actual_name}", end="")
                if not actual_name.endswith("Isolant") and not "Air" in actual_name:
                    print(" saved", flush=True)
                    extractBlock1 = ExtractBlock(registrationName=name, Input=geometry)
                    extractBlock1.Selectors = [block]
                    extractBlock1.UpdatePipeline()
                    extractSurface1 = ExtractSurface(
                        registrationName="ExtractSurface1", Input=extractBlock1
                    )

                    SaveData(
                        f"{basedir}/{actual_name}-deformed.stl", proxy=extractSurface1
                    )
                    Delete(extractBlock1)
                    del extractBlock1
                else:
                    print(" ignored", flush=True)

        suffix = "-deformed"
        cellsize_deformed = geometry
    elif found and dim == 2 and not axis:
        cellsize_deformed = deformed(cellsize, factor=1)

        # suffix = "-deformed"
        # cellsize = geometry

    # Views
    if args.views:
        if not args.field:
            for vkey in list(cellsize.PointData.keys()) + list(
                cellsize.CellData.keys()
            ):

                if not vkey in ignored_keys:
                    if vkey in list(cellsize.CellData.keys()):
                        color = ["CELLS", vkey]
                    if vkey in list(cellsize.PointData.keys()):
                        color = ["POINTS", vkey]

                    makeview(
                        args,
                        cellsize,
                        blockdata,
                        vkey,
                        fieldunits,
                        color,
                        basedir,
                        suffix="",
                        addruler=False,
                    )

            if found:
                for vkey in list(cellsize_deformed.PointData.keys()) + list(
                    cellsize_deformed.CellData.keys()
                ):

                    if not vkey in ignored_keys:
                        if vkey in list(cellsize_deformed.CellData.keys()):
                            color = ["CELLS", vkey]
                        if vkey in list(cellsize_deformed.PointData.keys()):
                            color = ["POINTS", vkey]

                        makeview(
                            args,
                            cellsize_deformed,
                            blockdata,
                            vkey,
                            fieldunits,
                            color,
                            basedir,
                            suffix=suffix,
                            addruler=False,
                        )

        else:
            if args.field in list(cellsize.PointData.keys()) + list(
                cellsize.CellData.keys()
            ):
                makeview(
                    args,
                    cellsize,
                    blockdata,
                    args.field,
                    fieldunits,
                    color,
                    basedir,
                    suffix="",
                    addruler=False,
                )

            if found and (
                args.field
                in list(cellsize_deformed.PointData.keys())
                + list(cellsize_deformed.CellData.keys())
            ):
                makeview(
                    args,
                    cellsize_deformed,
                    blockdata,
                    args.field,
                    fieldunits,
                    color,
                    basedir,
                    suffix=suffix,
                    addruler=False,
                )

    # for magnetfield:
    #   - view contour for magnetic potential (see pv-contours.py)
    #   - view glyph for MagneticField
    #   - compute spherical expansion coefficients

    # # save vtkjs

    # # export view
    # ExportView('/home/LNCMI-G/trophime/export-scene.vtkjs', view=renderView1, ParaViewGlanceHTML='')


if __name__ == "__main__":
    sys.exit(main())
