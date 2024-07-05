import numpy as np
import pandas as pd

from paraview.simple import *

import sys
import os

import argparse

######################## README ########################
#
# Create plots if paraview for 'display_results_v0.1.py'
#
# Options:
#   dir: add the directory where you want the results (must contain subdirectory like cooling/heatcorrelation/friction/np_...)
#   --np: specify nb of cores for results
#   --coolings: specify all coolings you want to display
#   --heatcorrelations: specify all heat correlations you want to display
#   --frictions: specify all frictions you want to display
#   --rint: specify inner radius (if 2D give 1, if axi give 2, one for each bitter)
#   --rint: specify outer radius (if 2D give 1, if axi give 2, one for each bitter)
#   --toolbox: specify toolbox (cfpdes or thermo_electric)
#
# axi: profiles for axi coord
#   --measures: give measures to add to plot (can be multiple)
#               -for "display_results_v0.1.py plot_profiles", need "Jth","temperature"-
#               -can also add B for profile on Z axis-
#   --Z: give one Z for each Bitter for 2ndtop & bottom profiles
#   --Zmax: give Zmax for B profile on Z axis
#
# 2D: profiles for 2D coord
#   --measures: give measures to add to plot (can be multiple)
#               -for "display_results_v0.1.py plot_profiles", need "temperature"-
#   --yaml_file: input yaml file of bitter for more precise boxplots
#   --nb: specify number of boxplots (only used if yaml_file is not specified)
#   --theta: specify angle of the slice (default pi/16)
#
########################################################
# python profiles.py axi bitters_old/M9Bitters-frans --rint 0.2 0.343 --rext 0.34 0.5 --Z 0.22292 0.250615 --Zmax 0.15 --coolings gradH meanH --frictions Colebrook --measures Jth temperature B
# python profiles.py 2D bitters_old/M9Bi-frans-nonlinear-ansys/ --yaml_file yaml/M9_Bi.yaml --rint 0.2 --rext 0.34 --coolings gradH meanH --frictions Colebrook
########################################################


def options(description: str, epilog: str):
    """
    define options
    """

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    subparsers = parser.add_subparsers(
        title="commands", dest="command", help="sub-command help"
    )
    parser_axi = subparsers.add_parser("axi", help="profiles for axi")
    parser_2D = subparsers.add_parser("2D", help="profiles for 2D")

    ## For all:
    for allparser in [parser_axi, parser_2D]:
        allparser.add_argument(
            "dir", help="specify measure directory", type=str, default=""
        )
        allparser.add_argument(
            "--coolings",
            help="what coolings do you want",
            nargs="+",
            metavar="coolings",
            type=str,
            choices=["mean", "meanH", "grad", "gradH", "gradHZ"],
            default=["mean", "meanH", "grad", "gradH"],
        )
        allparser.add_argument(
            "--heatcorrelations",
            help="what heat correlations do you want",
            nargs="+",
            metavar="heatcorrelations",
            type=str,
            choices=["Montgomery", "Dittus", "Colburn", "Silverberg"],
            default=["Montgomery"],
        )
        allparser.add_argument(
            "--frictions",
            help="what frictions do you want",
            nargs="+",
            metavar="frictions",
            type=str,
            choices=["Constant", "Blasius", "Filonenko", "Colebrook", "Swanee"],
            default=["Constant"],
        )
        allparser.add_argument(
            "--np", help="specify nb of cores for results", type=int, default=1
        )
        allparser.add_argument(
            "--rint",
            help="specify inner radius",
            nargs="+",
            metavar="rint",
            type=float,
            default=[0.2],
        )
        allparser.add_argument(
            "--rext",
            help="specify outer radius",
            nargs="+",
            metavar="rext",
            type=float,
            default=[0.34],
        )
        allparser.add_argument(
            "--toolbox",
            help="specify toolbox",
            choices=["cfpdes", "thermo_electric"],
            type=str,
            default="cfpdes",
        )

    parser_axi.add_argument(
        "--measures",
        help="what measures do you want",
        nargs="+",
        metavar="measures",
        type=str,
        default=["Jth", "temperature", "Vonmises"],
    )
    parser_axi.add_argument(
        "--Z",
        help="specify Z for 2nd profiles",
        nargs="+",
        metavar="Z",
        type=float,
        default=[0.22292],
    )
    parser_axi.add_argument("--Zmax", help="specify Zmax", type=float, default=None)

    parser_2D.add_argument(
        "--measures",
        help="what measures do you want",
        nargs="+",
        metavar="measures",
        type=str,
        default=["temperature", "Vonmises"],
    )
    parser_2D.add_argument(
        "--yaml_file", help="input yaml file", type=str, default=None
    )
    parser_2D.add_argument("--nb", help="specify nb of boxplots", type=int, default=10)
    parser_2D.add_argument(
        "--theta", help="specify angle of the slice", type=str, default="pi/16"
    )

    return parser


def profiles_2D(args, pwd: str, directory: str, exportcase):

    r = np.linspace(args.rint[0] * 1e3, args.rext[0] * 1e3, args.nb)
    if args.yaml_file:
        rslits = [args.rint[0] * 1e3]
        save_r = False
        odd_angle = 0
        save_angle = True
        with open(f"{pwd}/{args.yaml_file}", "r") as f:
            Lines = f.readlines()
            for line in Lines:
                if "!<Slit>" in line:
                    save_r = True
                if "angle:" in line and save_angle and save_r:
                    angle = float(line.replace("angle:", ""))
                    if angle != 0.0:
                        save_angle = False
                if "r:" in line and save_r:
                    rslits.append(float(line.replace("r:", "")))
                    save_r = False

        rslits.append(args.rext[0] * 1e3)

        # print(rslits)
        r = []
        for i in range(len(rslits) - 1):
            r.append(rslits[i])
            r.append(rslits[i] + (rslits[i + 1] - rslits[i]) / 2)

        r.append(rslits[-1])

    if "Vonmises" not in args.measures:
        print("r=", r)
        pts_array = []
        for measure in args.measures:
            if measure == "temperature":
                pts_array.append(f"{args.toolbox}.heat.temperature")
            else:
                pts_array.append(f"{args.toolbox}.expr.{measure}")
        exportcase.PointArrays = pts_array
        pts_array.insert(0, "arc_length")

        boxplots = pd.DataFrame()
        for x in r:
            # create a new 'Plot On Intersection Curves'
            plotOnIntersectionCurves1 = PlotOnIntersectionCurves(
                registrationName="PlotOnIntersectionCurves1", Input=exportcase
            )
            # Properties modified on plotOnIntersectionCurves1
            plotOnIntersectionCurves1.SliceType = "Cylinder"

            # Properties modified on plotOnIntersectionCurves1.SliceType
            plotOnIntersectionCurves1.SliceType.Axis = [0.0, 0.0, 1.0]
            plotOnIntersectionCurves1.SliceType.Radius = x * 1e-3

            # save data
            SaveData(
                f"{pwd}/{directory}/r_profile.csv",
                proxy=plotOnIntersectionCurves1,
                PointDataArrays=pts_array,
            )

            df = pd.read_csv(f"{pwd}/{directory}/r_profile.csv")
            df["r"] = [x] * len(df)
            boxplots = pd.concat([boxplots, df], axis=0)

            os.remove(f"{pwd}/{directory}/r_profile.csv")

            # set active source
            SetActiveSource(exportcase)
            # destroy plotOnIntersectionCurves1
            Delete(plotOnIntersectionCurves1)
            del plotOnIntersectionCurves1

        boxplots.to_csv(f"{pwd}/{directory}/r_profiles.csv")

        def rotate(xa, ya, xb, yb, angle):
            x1 = (np.cos(np.radians(angle)) * xa) + (np.sin(np.radians(angle)) * ya)
            y1 = (1 * np.sin(np.radians(angle)) * xa) - (np.cos(np.radians(angle)) * ya)

            x2 = (np.cos(np.radians(angle)) * xb) + (np.sin(np.radians(angle)) * yb)
            y2 = (1 * np.sin(np.radians(angle)) * xb) - (np.cos(np.radians(angle)) * yb)

            return [[x1, y1], [x2, y2]]

        print()

        plotOverLine1 = PlotOverLine(registrationName="PlotOverLine1", Input=exportcase)

        # init the 'Line' selected
        plotOverLine1.Point1 = [args.rint[0], 0.0, 0.0]
        plotOverLine1.Point2 = [args.rext[0], 0.0, 0.0]

        print(f"plotOverLine1.Point1 = [{args.rint[0]}, 0.0, 0.0]")
        print(f"plotOverLine1.Point2 = [{args.rext[0]},  0.0, 0.0]")

        # save data
        SaveData(
            f"{pwd}/{directory}/center_profile.csv",
            proxy=plotOverLine1,
            PointDataArrays=pts_array,
        )

        # set active source
        SetActiveSource(exportcase)
        # destroy plotOverLine1
        Delete(plotOverLine1)
        del plotOverLine1

        if "pi" in args.theta:
            theta = np.pi / float(args.theta.replace("pi/", "")) / 2
        else:
            theta = float(args.theta.replace("pi/", "")) / 2

        plotOverLine1 = PlotOverLine(registrationName="PlotOverLine1", Input=exportcase)

        # init the 'Line' selected

        print(
            f"plotOverLine1.Point1 = [{args.rint[0]* np.cos(theta / 2.0)}, {args.rint[0]* np.sin(theta / 2.0)}, 0.0]"
        )
        print(
            f"plotOverLine1.Point2 = [{args.rext[0] * np.cos(theta / 2.0)},  {args.rext[0]* np.sin(theta / 2.0)}, 0.0]"
        )

        plotOverLine1.Point1 = [
            args.rint[0] * np.cos(theta / 2.0),
            args.rint[0] * np.sin(theta / 2.0),
            0.0,
        ]
        plotOverLine1.Point2 = [
            args.rext[0] * np.cos(theta / 2.0),
            args.rext[0] * np.sin(theta / 2.0),
            0.0,
        ]

        # save data
        SaveData(
            f"{pwd}/{directory}/border_profile.csv",
            proxy=plotOverLine1,
            PointDataArrays=pts_array,
        )

        # set active source
        SetActiveSource(exportcase)
        # destroy plotOverLine1
        Delete(plotOverLine1)
        del plotOverLine1

        [[x1, y1], [x2, y2]] = rotate(args.rint[0], 0.0, args.rext[0], 0.0, 3 * angle)

        plotOverLine1 = PlotOverLine(registrationName="PlotOverLine1", Input=exportcase)

        # init the 'Line' selected

        print(f"plotOverLine1.Point1 = [{x1}, {y1}, 0.0]")
        print(f"plotOverLine1.Point2 = [{x2},  {y2}, 0.0]")

        plotOverLine1.Point1 = [
            x1,
            y1,
            0.0,
        ]
        plotOverLine1.Point2 = [
            x2,
            y2,
            0.0,
        ]

        # save data
        SaveData(
            f"{pwd}/{directory}/odd_profile.csv",
            proxy=plotOverLine1,
            PointDataArrays=pts_array,
        )
        # set active source
        SetActiveSource(exportcase)
        # destroy plotOverLine1
        Delete(plotOverLine1)
        del plotOverLine1
    else:
        print("r=", r)
        pts_array = []
        for measure in args.measures:
            if measure == "temperature":
                pts_array.append(f"{args.toolbox}.heat.temperature")
            else:
                pts_array.append(f"{args.toolbox}.expr.{measure}")
        exportcase.PointArrays = pts_array
        pts_array.insert(0, "arc_length")

        def rotate(xa, ya, xb, yb, angle):
            x1 = (np.cos(np.radians(angle)) * xa) + (np.sin(np.radians(angle)) * ya)
            y1 = (1 * np.sin(np.radians(angle)) * xa) - (np.cos(np.radians(angle)) * ya)

            x2 = (np.cos(np.radians(angle)) * xb) + (np.sin(np.radians(angle)) * yb)
            y2 = (1 * np.sin(np.radians(angle)) * xb) - (np.cos(np.radians(angle)) * yb)

            return [[x1, y1], [x2, y2]]

        if "pi" in args.theta:
            theta = np.pi / float(args.theta.replace("pi/", "")) / 2
        else:
            theta = float(args.theta.replace("pi/", "")) / 2

        plotOverLine1 = PlotOverLine(registrationName="PlotOverLine1", Input=exportcase)

        [[x1, y1], [x2, y2]] = rotate(args.rint[0], 0.0, args.rext[0], 0.0, 4 * angle)

        print(f"plotOverLine1.Point1 = [{x1}, {y1}, 0.0]")
        print(f"plotOverLine1.Point2 = [{x2},  {y2}, 0.0]")

        plotOverLine1.Point1 = [
            x1,
            y1,
            0.0,
        ]
        plotOverLine1.Point2 = [
            x2,
            y2,
            0.0,
        ]

        # save data
        SaveData(
            f"{pwd}/{directory}/center_profile.csv",
            proxy=plotOverLine1,
            PointDataArrays=pts_array,
        )
        # set active source
        SetActiveSource(exportcase)
        # destroy plotOverLine1
        Delete(plotOverLine1)
        del plotOverLine1

        plotOverLine1 = PlotOverLine(registrationName="PlotOverLine1", Input=exportcase)

        [[x1, y1], [x2, y2]] = rotate(args.rint[0], 0.0, args.rext[0], 0.0, 8 * angle)

        print(f"plotOverLine1.Point1 = [{x1}, {y1}, 0.0]")
        print(f"plotOverLine1.Point2 = [{x2},  {y2}, 0.0]")

        plotOverLine1.Point1 = [
            x1,
            y1,
            0.0,
        ]
        plotOverLine1.Point2 = [
            x2,
            y2,
            0.0,
        ]

        # save data
        SaveData(
            f"{pwd}/{directory}/border_profile.csv",
            proxy=plotOverLine1,
            PointDataArrays=pts_array,
        )

        # set active source
        SetActiveSource(exportcase)
        # destroy plotOverLine1
        Delete(plotOverLine1)
        del plotOverLine1

        [[x1, y1], [x2, y2]] = rotate(args.rint[0], 0.0, args.rext[0], 0.0, 7 * angle)

        plotOverLine1 = PlotOverLine(registrationName="PlotOverLine1", Input=exportcase)

        # init the 'Line' selected

        print(f"plotOverLine1.Point1 = [{x1}, {y1}, 0.0]")
        print(f"plotOverLine1.Point2 = [{x2},  {y2}, 0.0]")

        plotOverLine1.Point1 = [
            x1,
            y1,
            0.0,
        ]
        plotOverLine1.Point2 = [
            x2,
            y2,
            0.0,
        ]

        # save data
        SaveData(
            f"{pwd}/{directory}/odd_profile.csv",
            proxy=plotOverLine1,
            PointDataArrays=pts_array,
        )
        # set active source
        SetActiveSource(exportcase)
        # destroy plotOverLine1
        Delete(plotOverLine1)
        del plotOverLine1

    return 0


def profiles_axi(args, pwd: str, directory: str, exportcase):
    pts_array = []
    cell_array = []
    for measure in args.measures:
        if measure == "temperature":
            pts_array.append(f"{args.toolbox}.heat.temperature")
        if measure == "B":
            cell_array.append("cfpdes.expr.B")
            exportcase.CellArrays = cell_array
        else:
            pts_array.append(f"{args.toolbox}.expr.{measure}")
    exportcase.PointArrays = pts_array
    pts_array.insert(0, "arc_length")

    if pts_array:

        pts = [
            [
                [args.rint[0], -args.Z[0], 0.0],
                [args.rext[0], -args.Z[0], 0.0],
                "inner_2ndtop",
            ],
            [[args.rint[0], 0.0, 0.0], [args.rext[0], 0.0, 0.0], "inner_core"],
            [
                [args.rint[0], args.Z[0], 0.0],
                [args.rext[0], args.Z[0], 0.0],
                "inner_2ndbot",
            ],
            [
                [args.rint[1], -args.Z[1], 0.0],
                [args.rext[1], -args.Z[1], 0.0],
                "outer_2ndtop",
            ],
            [[args.rint[1], -0.0, 0.0], [args.rext[1], -0.0, 0.0], "outer_core"],
            [
                [args.rint[1], args.Z[1], 0.0],
                [args.rext[1], args.Z[1], 0.0],
                "outer_2ndbot",
            ],
        ]

        for p in pts:
            # create a new 'Plot Over Line'
            plotOverLine1 = PlotOverLine(
                registrationName="PlotOverLine1", Input=exportcase
            )

            print(f"Line {p[2]}: p1={p[0]}; p2={p[1]}")

            # init the 'Line' selected
            plotOverLine1.Point1 = p[0]
            plotOverLine1.Point2 = p[1]

            # save data
            SaveData(
                f"{pwd}/{directory}/{p[2]}_profile.csv",
                proxy=plotOverLine1,
                PointDataArrays=pts_array,
            )

            # set active source
            SetActiveSource(exportcase)

            # destroy plotOverLine1
            Delete(plotOverLine1)
            del plotOverLine1

    if cell_array and args.Zmax:
        cell_array.insert(0, "arc_length")
        # create a new 'Plot Over Line'
        plotOverLine1 = PlotOverLine(registrationName="PlotOverLine1", Input=exportcase)

        print(f"Line B: p1={[0, -args.Zmax, 0.0]}; p2={[0, args.Zmax, 0.0]}")
        # init the 'Line' selected for 'Source'
        plotOverLine1.Point1 = [0, -args.Zmax, 0.0]
        plotOverLine1.Point2 = [0, args.Zmax, 0.0]

        # save data
        SaveData(
            f"{pwd}/{directory}/B_profile.csv",
            proxy=plotOverLine1,
            PointDataArrays=cell_array,
        )

        # set active source
        SetActiveSource(exportcase)

        # destroy plotOverLine1
        Delete(plotOverLine1)
        del plotOverLine1

    return 0


def main():

    description = "export.case"
    epilog = "\n"
    parser = options(description, epilog)
    args = parser.parse_args()

    pwd = os.getcwd()

    for cooling in args.coolings:
        heatcorrelations = args.heatcorrelations
        frictions = args.frictions
        if cooling in ["mean", "grad"]:
            heatcorrelations = ["Montgomery"]
            frictions = ["Constant"]

        for heatcorrelation in heatcorrelations:
            for friction in frictions:

                print(f"\n## {cooling} {heatcorrelation} {friction}")

                if args.toolbox == "thermo_electric":
                    export_dir = "thermoelectric.exports/Export.case"
                else:
                    export_dir = f"{args.toolbox}.exports/Export.case"

                directory = (
                    f"{args.dir}/{cooling}/{heatcorrelation}/{friction}/np_{args.np}"
                )

                # create a new 'EnSight Reader'
                exportcase = EnSightReader(
                    registrationName="Export.case",
                    CaseFileName=f"{pwd}/{directory}/{export_dir}",
                )

                if args.command == "2D":
                    profiles_2D(args, pwd, directory, exportcase)

                if args.command == "axi":
                    profiles_axi(args, pwd, directory, exportcase)

                # destroy exportcase
                Delete(exportcase)
                del exportcase

    return 0


if __name__ == "__main__":
    sys.exit(main())
