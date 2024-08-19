import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json

from tabulate import tabulate
from scipy.optimize import curve_fit
from natsort import natsorted

import sys
import os

import argparse

######################## README ########################
#
# Compare feelpp results with Ansys or SimpleModel results
#
# Options:
#   --mtype: bitter or helix
#   --filter: magnet filter if measures from MSite
#   --dir: add the directory where you want the results (must contain subdirectory like cooling/heatcorrelation/friction/np_...)
#   --np: specify nb of cores for results
#   --pupitrefile: txt file of exp
#   --coolings: specify all coolings you want to display
#   --heatcorrelations: specify all heat correlations you want to display
#   --frictions: specify all frictions you want to display
#   --oneI: IF RESULTS ARE NOT COMMISSIONING, ONLY ONE I
#   --save: save plots & err tables
#   --show: show plots
#
# plot_measure: plot measures from measures.csv vs I & give relative error
#   --measures: choose from ["Ucoil","B0","Power","flow","R(I)","MinTH","MaxTH","MeanTH","TH","Tout"] (can be multiple)
#
# plot_values
#   --measures: choose from ["Uw","DT","Flux","HeatCoeff","PowerH","statsTH","cf"] (can be multiple)
#   --current: set specific currents to plot out of the hole commissioning
#   --modeldir: specify external model directory with : [..].measures/values.csv
#
# plot_profiles:
#   files: specify simple model files containing temperature profiles
#   --measures: choose from ["T","J"] (can be multiple)
#   --parts: choose on which parts of the magnet the profile will be plot; choices=["2ndtop", "core", "2ndbot"] (can be multiple)
#   --toolbox: specify toolbox used between ["cfpdes", "thermo_electric"]
#   --hideSimpleModel: hide the SimpleModel profiles
#   --boxplotdir: specify 2D measures directory
#   --show2Dlines: show 2D profiles
#   --Ansys: specify Ansys file containing temperature profiles
#
########################################################
# python display_results_v0.1.py plot_measures --mtype bitter --dir bitters_new/M9Bitters_18MWcommi9june/ --coolings meanH gradH gradHZ --frictions Colebrook  --pupitrefile pupitre/M9_2023.06.09---13\:48\:53_from13\:54\:09to14\:00\:05.txt --save --show --measures "Ucoil" "B0" "Power" "flow" "R(I)" "MinTH" "MaxTH" "MeanTH" "TH" Tout
# python display_results_v0.1.py plot_values --mtype bitter --dir bitters_new/M9Bitters_18MW25oct/ --coolings meanH gradH gradHZ --frictions Colebrook --oneI --pupitrefile pupitre/M9_2023.10.25---15\:50\:45_to16\:08\:42.txt --save --show --measures "Uw" "DT" "Flux" "HeatCoeff" "PowerH" "statsTH" "cf" --modeldir Exported_Results_Modified_Inner_Bitter/Exported_Results/
# python display_results_v0.1.py plot_profiles '{"inner":"Exported_Results/20231127-LNCMI_Inner_Original-Simple_Model_Data.csv","outer":"Exported_Results/20231127-LNCMI_Outer_Original-Simple_Model_Data.csv"}' --coolings mean grad --frictions Constant --dir bitters_old/M9Bitters-frans --save --show --measures T J
# python display_results_v0.1.py plot_profiles '{"inner":"Exported_Results/20231127-LNCMI_Inner_Original-Simple_Model_Data.csv"}' --part core --coolings meanH --frictions Colebrook --dir bitters_old/M9Bitters-frans --boxplotsdir bitters_old/M9Bi-frans-nonlinear-ansys/ --save --show --show2Dlines --Ansys Exported_Results_fix/20240129-LNCMI_Inner_Original-Cross-sectional_data.csv
########################################################


def options(description: str, epilog: str):
    """
    define options
    """

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    subparsers = parser.add_subparsers(
        title="commands", dest="command", help="sub-command help"
    )
    parser_plotm = subparsers.add_parser("plot_measures", help="plot measures.csv vs I")
    parser_plotv = subparsers.add_parser("plot_values", help="plot values.csv")
    parser_plotp = subparsers.add_parser(
        "plot_profiles", help="plot temperature profile"
    )

    ## For all:
    for allparser in [parser_plotm, parser_plotv, parser_plotp]:
        allparser.add_argument(
            "--mtype",
            help="input magnet type",
            type=str,
            choices=["helix", "bitter"],  # , "site"],
            default="bitter",
        )
        allparser.add_argument(
            "--filter",
            help="input magnet filter name",
            type=str,
            default="",
        )
        allparser.add_argument(
            "--dir", help="specify measure directory", type=str, default=""
        )
        allparser.add_argument(
            "--np", help="specify nb of cores for results", type=int, default=1
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
            "--FF", help="input Field Factor", type=float, default=None
        )
        allparser.add_argument("--Tin", help="specify Tinit", type=str, default=None)
        allparser.add_argument("--oneI", help="measures for 1 I", action="store_true")

        allparser.add_argument(
            "--magfile",
            help="magfile from control command (ex. MAGFile.conf)",
            type=str,
            default=None,
        )
        allparser.add_argument(
            "--pupitrefile",
            help="input pupitre file",
            type=str,
            default=None,
        )
        allparser.add_argument("--debug", help="activate debug", action="store_true")
        allparser.add_argument("--save", help="save Tables", action="store_true")
        allparser.add_argument("--show", help="show graph", action="store_true")

    ## For plot_measures
    parser_plotm.add_argument(
        "--measures",
        help="what measures do you want",
        nargs="+",
        metavar="measures",
        type=str,
        choices=[
            "Ucoil",
            "B0",
            "Power",
            "flow",
            "R(I)",
            "MinTH",
            "MaxTH",
            "MeanTH",
            "TH",
            "Tout",
        ],
        default=[
            "Ucoil"
        ],  # "Ucoil" "B0" "Power" "flow" "R(I)" "MinTH" "MaxTH" "MeanTH" "TH"
    )

    ## For plot values:
    parser_plotv.add_argument(
        "--measures",
        help="what measures do you want",
        nargs="+",
        metavar="measures",
        type=str,
        choices=[
            "Uw",
            "DT",
            "Flux",
            "HeatCoeff",
            "PowerH",
            "statsTH",
            "cf",
        ],
        default=["Uw"],  # "Uw" "DT" "Flux" "HeatCoeff" "PowerH" "statsTH" "cf"
    )
    parser_plotv.add_argument(
        "--current",
        help="set a specific current",
        nargs="+",
        metavar="currents",
        type=int,
        default=None,
    )
    parser_plotv.add_argument(
        "--modeldir", help="specify model dir directory", type=str, default=None
    )

    ## For plot_profiles
    parser_plotp.add_argument(
        "files",
        help="specify cooling files",
        type=json.loads,  #'{"mean":{"file":"values.csv"}}'
    )
    parser_plotp.add_argument(
        "--measures",
        help="what measures do you want",
        nargs="+",
        metavar="measures",
        type=str,
        choices=["T", "J", "Vonmises", "HoopStress"],
        default=["T"],
    )
    parser_plotp.add_argument(
        "--parts",
        help="give parts top plot",
        nargs="+",
        metavar="parts",
        type=str,
        choices=["2ndtop", "core", "2ndbot"],
        default=["2ndtop", "core", "2ndbot"],
    )
    parser_plotp.add_argument(
        "--toolbox",
        help="specify toolbox",
        choices=["cfpdes", "thermo_electric"],
        type=str,
        default="cfpdes",
    )
    parser_plotp.add_argument(
        "--hideSimpleModel", help="show Simple model profiles", action="store_true"
    )
    parser_plotp.add_argument(
        "--hideAxi", help="show Axi profiles", action="store_true"
    )
    parser_plotp.add_argument(
        "--boxplotsdir", help="specify 2D measure directory", type=str, default=None
    )
    parser_plotp.add_argument(
        "--boxplotsdir2", help="specify 2D measure directory", type=str, default=None
    )
    parser_plotp.add_argument(
        "--show2Dlines", help="show 2D profiles", action="store_true"
    )
    parser_plotp.add_argument(
        "--Ansys", help="specify ansys profiles file", type=str, default=None
    )
    parser_plotp.add_argument(
        "--Ansysvm", help="specify ansys vonmises profiles file", type=str, default=None
    )

    return parser


def polyfit(x, a0, a1, a2):
    return a0 + a1 * x + a2 * x**2


def plot_greySlits(df: pd.DataFrame, cx: str, cy: str, ax, legend):
    nan_positions = df[cy].isnull()

    # Fill areas before and after NaN values
    x1 = None
    x2 = None
    for i in range(len(nan_positions)):
        if nan_positions[i] and i > 0 and not x1:
            x1 = i - 1
        elif not nan_positions[i] and i > 0:
            if nan_positions[i - 1]:
                x2 = i

        if x1 and x2:
            ax.axvspan(
                df[cx][x1],
                df[cx][x2],
                alpha=0.3,
                color="grey",
            )
            x1 = None
            x2 = None
            legend.append("_nolegend_")
    return legend


def compare(
    df1: pd.DataFrame,
    cx1: str,
    cy1: str,
    df2: pd.DataFrame,
    cx2: str,
    cy2: str,
):
    # blank_rows = df1.isnull().any(axis=1)
    # rows_to_remove = blank_rows | blank_rows.shift(-5) | blank_rows.shift(-6)
    # df1 = df1[~rows_to_remove].reset_index(drop=True)

    x1 = df1[cx1].to_numpy()
    bx1 = df1[cy1].to_numpy()
    y1 = pd.Series(bx1, index=x1)

    # blank_rows = df2.isnull().any(axis=1)
    # rows_to_remove = blank_rows | blank_rows.shift(-1) | blank_rows.shift(-2)
    # df2 = df2[~rows_to_remove].reset_index(drop=True)
    x2 = df2[cx2].to_numpy()
    bx2 = df2[cy2].to_numpy()
    y2 = pd.Series(bx2, index=x2)

    df = pd.DataFrame({f"calc1": y1, f"calc2": y2})

    df["calc1"].interpolate(method="linear", direction="forward", inplace=True)
    df["calc2"].interpolate(method="linear", direction="forward", inplace=True)
    df["calc1-calc2"] = np.sqrt(np.square(df["calc1"] - df["calc2"]))
    L2ErrNorm = np.linalg.norm(df["calc1"] - df["calc2"])
    L2ErrRelaNorm = np.linalg.norm(df["calc1"] - df["calc2"]) / np.linalg.norm(
        df["calc1"]
    )
    maxErr = np.max(df["calc1-calc2"])
    max_error_rows = df[df["calc1-calc2"] == maxErr]
    # print(max_error_rows)
    max_error_x_values = max_error_rows.index.tolist()

    return (L2ErrNorm, L2ErrRelaNorm, maxErr, max_error_x_values[0])


def plot_measures(args, plotm_show):
    color_cycler = [
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#F0E442",
        "c",
        "m",
        "y",
        "k",
    ]
    linestyle_cycler = [
        "-",
        "--",
        "-.",
        ":",
        "-",
        "--",
        "-.",
        "-",
        "--",
        "-.",
        ":",
        "-",
        "--",
        "-.",
    ]
    pointstyle_cycler = [
        "o",
        "^",
        "s",
        "h",
        "v",
        "p",
        "*",
        "o",
        "^",
        "s",
        "h",
        "v",
        "p",
        "*",
    ]

    filter = args.filter

    Be = 16
    if args.pupitrefile:
        ## read exp files
        pupitre = pd.read_csv(args.pupitrefile, sep=r"\s+")

        # pupitre["Time"] = pd.to_datetime(pupitre["Date"] + " " + pupitre["Time"])
        # pupitre["Time"] = pd.to_datetime(pupitre["Time"])

        ## Temperature C to K
        pupitre["Tout"] = pupitre["Tout"] + 273.2
        col = [int(pup[4:]) for pup in pupitre.columns if "Tcal" in pup]
        Be = max(col)
        for i in range(1, Be + 1):
            pupitre[f"Tcal{i}"] = pupitre[f"Tcal{i}"] + 273.2
            ## Add R(I)= Ucoil/Icoil ??
            pupitre[f"Rcoil{i}"] = pupitre[f"Ucoil{i}"] / pupitre[f"Icoil{i}"]
        pupitre["Tin1"] = pupitre["Tin1"] + 273.2
        pupitre["Tin2"] = pupitre["Tin2"] + 273.2
        ## Calcul Power for helix & bitter
        if args.mtype == "bitter":
            pupitre["Power_bitter"] = (
                (pupitre[f"Ucoil{Be-1}"] + pupitre[f"Ucoil{Be}"])
                * pupitre[f"Icoil{Be}"]
                / 1e6
            )
        elif args.mtype == "helix":
            pupitre["Power_helix"] = (
                sum(pupitre[f"Ucoil{x}"] for x in range(1, 8)) * pupitre["Icoil7"] / 1e6
            )

        if args.FF:
            pupitre["Field_bitter"] = args.FF * pupitre[f"Icoil{Be}"]
            pupitre["Field_helix"] = args.FF * pupitre["Icoil7"]
        else:
            pupitre["Field_bitter"] = pupitre["Field"]
            pupitre["Field_helix"] = pupitre["Field"]

    # ## read simu dict and create a dict for df
    simu = {}
    configs = []
    for cooling in args.coolings:
        heatcorrelations = args.heatcorrelations
        frictions = args.frictions
        if cooling in ["mean", "grad"]:
            heatcorrelations = ["Montgomery"]
            frictions = ["Constant"]
        for heatcorrelation in heatcorrelations:
            for friction in frictions:
                config = f"{cooling}_{heatcorrelation}_{friction}"
                configs.append(config)
                if not args.oneI:
                    simu[config] = pd.read_csv(
                        f"{args.dir}/{cooling}/{heatcorrelation}/{friction}/np_{args.np}/measures.csv",
                        sep=r",",
                    )
                else:
                    simu[config] = pd.read_csv(
                        f"{args.dir}/{cooling}/{heatcorrelation}/{friction}/np_{args.np}/measures.csv",
                        sep=r",",
                        index_col=0,
                    )
                    simu[config] = simu[config].T
                simu[config]["B0[T]"] = abs(simu[config]["B0[T]"])

    ## stock a df for column names
    last_config = simu[config]

    color_cooling = {
        "mean": color_cycler[0],
        "meanH": color_cycler[1],
        "grad": color_cycler[2],
        "gradH": color_cycler[3],
        "gradHZ": color_cycler[4],
    }
    line_style = {
        "MinTH": "^--",
        "MaxTH": "s--",
        "MeanTH": "o-",
    }

    if args.pupitrefile:
        ## Create a translation dict for simu <-> exp
        Icoils = {"helix": "Icoil7", "bitter": f"Icoil{Be}"}
        dict_trad = {
            "Ucoil": {"unit": "V", "exp": "Ucoil"},
            "R(I)": {"unit": "ohm", "exp": ""},
            "MinTH": {"unit": "K", "exp": "Tcal", "line": "."},
            "MaxTH": {"unit": "K", "exp": "Tcal", "line": ".-"},
            "MeanTH": {"unit": "K", "exp": "Tcal", "line": "--"},
            "B0": {"unit": "T"},
            "flow": {"unit": "l/s"},
            "Tout": {"unit": "K"},
            "Power": {"unit": "MW"},
        }

        # dict_measures = {
        #     "B0": {
        #         "helix": {"B0[T]": "Field"},
        #         "bitter": {"B0[T]": "Field"},
        #     }
        # }
        dict_measures = {
            "B0": {
                "helix": {"B0[T]": "Field_helix"},
                "bitter": {"B0[T]": "Field_bitter"},
            }
        }
        for i in ["Ucoil", "R(I)", "MinTH", "MaxTH", "MeanTH"]:  # "R(I)",
            dict_measures[i] = {
                "helix": {
                    f"{filter}{i}_H1H2[{dict_trad[i]['unit']}]": f"{dict_trad[i]['exp']}1",
                    f"{filter}{i}_H3H4[{dict_trad[i]['unit']}]": f"{dict_trad[i]['exp']}2",
                    f"{filter}{i}_H5H6[{dict_trad[i]['unit']}]": f"{dict_trad[i]['exp']}3",
                    f"{filter}{i}_H7H8[{dict_trad[i]['unit']}]": f"{dict_trad[i]['exp']}4",
                    f"{filter}{i}_H9H10[{dict_trad[i]['unit']}]": f"{dict_trad[i]['exp']}5",
                    f"{filter}{i}_H11H12[{dict_trad[i]['unit']}]": f"{dict_trad[i]['exp']}6",
                    f"{filter}{i}_H13H14[{dict_trad[i]['unit']}]": f"{dict_trad[i]['exp']}7",
                },
                "bitter": {
                    f"{filter}{i}_{filter}M8_Bi[{dict_trad[i]['unit']}]": f"{dict_trad[i]['exp']}{Be-1}",
                    f"{filter}{i}_{filter}M8_Be[{dict_trad[i]['unit']}]": f"{dict_trad[i]['exp']}{Be}",
                    f"{filter}{i}_{filter}M9_Bi[{dict_trad[i]['unit']}]": f"{dict_trad[i]['exp']}{Be-1}",
                    f"{filter}{i}_{filter}M9_Be[{dict_trad[i]['unit']}]": f"{dict_trad[i]['exp']}{Be}",
                    f"{filter}{i}_{filter}M10_Bi[{dict_trad[i]['unit']}]": f"{dict_trad[i]['exp']}{Be-1}",
                    f"{filter}{i}_{filter}M10_Be[{dict_trad[i]['unit']}]": f"{dict_trad[i]['exp']}{Be}",
                },
            }
        dict_measures["TH"] = dict_measures["MinTH"]
        dict_measures["Power"] = {
            "helix": {f"{filter}PowerM[MW]": "Power_helix"},
            "bitter": {f"{filter}PowerM[MW]": "Power_bitter"},
        }
        if "M9" in args.pupitrefile:
            dict_measures["flow"] = {
                "helix": {f"{filter}flow[l/s]": "Flow1"},
                "bitter": {f"{filter}flow[l/s]": "Flow2"},
            }
        else:
            dict_measures["flow"] = {
                "helix": {f"{filter}flow[l/s]": "Flow2"},
                "bitter": {f"{filter}flow[l/s]": "Flow1"},
            }
        dict_measures["Tout"] = {
            "helix": {f"{filter}Tout[K]": "Tout"},
            "bitter": {f"{filter}Tout[K]": "Tout"},
        }

    Table = {}
    for meas in args.measures:  ## go through each measure
        Table[meas] = pd.DataFrame()  ## 1 df per measure
        Table[meas][f"{filter}I[A]"] = last_config[
            f"{filter}I[A]"
        ]  ## put all I from simu

        ## show all "TH" at the same time
        if meas == "TH":
            if args.pupitrefile:
                ## translation dict for measure
                given_meas = dict_measures[meas][args.mtype]
            for column in last_config:
                if f"{filter}MinTH" in column:
                    legend = []

                    if plotm_show:
                        if args.pupitrefile:
                            ## pupitre plot measure vs I
                            ax = pupitre.plot(
                                x=Icoils[args.mtype],
                                y=given_meas[column],
                                ylabel=column.replace("MinTH", "TH").replace(
                                    filter, ""
                                ),
                                xlabel=f"{filter}I[A]",
                                color=color_cycler[0],
                            )
                            legend.append(given_meas[column])
                        else:
                            fig, ax = plt.subplots()

                    ## browse all config df
                    for config in simu:
                        ## simu plot measure vs I
                        for m in ["MinTH", "MeanTH", "MaxTH"]:
                            selected_columns = simu[config][
                                [f"{filter}I[A]", column.replace("MinTH", m)]
                            ]
                            selected_columns = selected_columns.rename(
                                columns={
                                    column.replace(
                                        "MinTH", m
                                    ): f"{config}_{column.replace('MinTH', m)}"
                                }
                            )
                            Table[meas] = pd.merge(
                                Table[meas], selected_columns, on=f"{filter}I[A]"
                            )
                            if plotm_show:
                                simu[config].plot(
                                    x=f"{filter}I[A]",
                                    y=column.replace("MinTH", m),
                                    ax=ax,
                                    ylabel=column.replace("MinTH", "TH").replace(
                                        filter, ""
                                    ),
                                    title=column.replace("MinTH", "TH").replace(
                                        filter, ""
                                    ),
                                    style=line_style[m],
                                    c=color_cooling[config.split("_")[0]],
                                    markerfacecolor="none",
                                )
                                legend.append(f"{config}_{m}")
                    if plotm_show:
                        ax.legend(legend)
                        ax.grid(True)

        elif meas == "R(I)":
            filter_meas = f"{filter}R_"
            Table[meas] = pd.DataFrame()
            Table[meas]["configs"] = configs
            for column in last_config:
                if filter_meas in column:
                    if plotm_show:
                        legend = []
                        fig, ax = plt.subplots()

                        if "MeanTH" in args.measures:
                            legend2 = []
                            fig2, ax2 = plt.subplots()

                            meanth_col = column.replace("R", "MeanTH").replace(
                                "[ohm]", "[K]"
                            )

                    # U_col = column.replace("R(I)", "Ucoil").replace("[ohm]", "[V]")
                    a0 = []
                    a1 = []
                    a2 = []

                    ## browse all cooling df
                    color_index = 0
                    for config in simu:
                        # simu[cooling][column] = (
                        #     simu[cooling][U_col] / simu[cooling][f"{filter}I[A]"]
                        # )
                        if not args.oneI:
                            popt, pcov = curve_fit(
                                polyfit,
                                simu[config][f"{filter}I[A]"].values,
                                simu[config][column].values,
                            )

                            fitlegend = "Num: a0=%5.3e, a1=%5.3e, a2=%5.3e" % tuple(
                                popt
                            )
                            print(tuple(popt))
                            print("fit[", column, " =", fitlegend, "]")
                            a0.append(tuple(popt)[0])
                            a1.append(tuple(popt)[1])
                            a2.append(tuple(popt)[2])

                        if plotm_show:
                            # simu plot measure vs I
                            simu[config].plot(
                                x=f"{filter}I[A]",
                                y=column,
                                ax=ax,
                                ylabel=column.replace(filter, ""),
                                title=column.replace(filter, ""),
                                style=pointstyle_cycler[color_index],
                                markerfacecolor="none",
                                color=color_cycler[color_index + 2],
                            )
                            legend.append(f"{config}_R")

                            if not args.oneI:
                                I_data = np.linspace(
                                    simu[config][f"{filter}I[A]"].min(),
                                    simu[config][f"{filter}I[A]"].max(),
                                    100,
                                )
                                ax.plot(
                                    I_data,
                                    polyfit(I_data, *popt),
                                    color=color_cycler[color_index],
                                    linestyle=linestyle_cycler[color_index],
                                )
                                legend.append(f"{config}_polyfit")

                            if "MeanTH" in args.measures:
                                # simu plot measure vs meanTH
                                simu[config].plot(
                                    x=meanth_col,
                                    y=column,
                                    ax=ax2,
                                    ylabel=column.replace(filter, ""),
                                    style=pointstyle_cycler[color_index],
                                    markerfacecolor="none",
                                    color=color_cycler[color_index],
                                )
                                legend2.append(f"{config}_R")

                            color_index += 1

                    if not args.oneI:
                        Table[meas][f"a0_{column}"] = a0
                        Table[meas][f"a1_{column}"] = a1
                        Table[meas][f"a2_{column}"] = a2

                    if plotm_show:
                        ax.legend(legend)
                        ax.grid(True)
                        if "MeanTH" in args.measures:
                            ax2.legend(legend2)
                            ax2.grid(True)

        else:
            if args.pupitrefile:
                ## translation dict for measure
                given_meas = dict_measures[meas][args.mtype]

                ## find the closest I value in pupitre file for each I of the simu
                closest_values = [
                    pupitre[Icoils[args.mtype]].iloc[
                        (pupitre[Icoils[args.mtype]] - value).abs().idxmin()
                    ]
                    for value in Table[meas][f"{filter}I[A]"]
                ]

                ## select rows in pupitre file where Icoil column is in the list of closest values (only one per closest value)
                selected_rows = (
                    pupitre[pupitre[Icoils[args.mtype]].isin(closest_values)]
                    .drop_duplicates(subset=Icoils[args.mtype], keep="first")
                    .copy()
                )

                ## add a new column 'I[A]' with the corresponding I closest to the closest value
                selected_rows[f"{filter}I[A]"] = [
                    min(
                        Table[meas][f"{filter}I[A]"],
                        key=lambda x: abs(x - closest_value),
                    )
                    for closest_value in selected_rows[Icoils[args.mtype]]
                ]

            filter_meas = f"{filter}{meas}"
            if meas == "B0":
                filter_meas = meas

            ## find columns corresponding to measure
            for column in last_config:
                if filter_meas in column:
                    legend = []

                    if args.pupitrefile:
                        ## select columns corresponding to measure
                        selected_columns = selected_rows[
                            [f"{filter}I[A]", given_meas[column]]
                        ]
                        selected_columns = selected_columns.rename(
                            columns={
                                given_meas[
                                    column
                                ]: f"{given_meas[column]}[{dict_trad[meas]['unit']}]"
                            }
                        )
                        ## merge it with meaure table on I
                        Table[meas] = pd.merge(
                            Table[meas], selected_columns, on=f"{filter}I[A]"
                        )
                        if plotm_show:
                            ## pupitre plot measure vs I
                            ax = pupitre.plot(
                                x=Icoils[args.mtype],
                                y=given_meas[column],
                                ylabel=column.replace(filter, ""),
                                title=column.replace(filter, ""),
                                xlabel=f"{filter}I[A]",
                                color=color_cycler[0],
                            )
                            legend.append(given_meas[column])
                    elif plotm_show:
                        fig, ax = plt.subplots()

                    ## browse all cooling df
                    color_index = 1
                    for config in simu:
                        ## same thing as the pupitr file merge it to measure table
                        selected_columns = simu[config][[f"{filter}I[A]", column]]
                        selected_columns = selected_columns.rename(
                            columns={column: f"{config}_{column}"}
                        )
                        Table[meas] = pd.merge(
                            Table[meas], selected_columns, on=f"{filter}I[A]"
                        )
                        ## add relative error simu exp
                        if args.pupitrefile:
                            Table[meas][f"err_{config}_{given_meas[column]}(%)"] = (
                                round(
                                    abs(
                                        Table[meas][f"{config}_{column}"]
                                        - Table[meas][
                                            f"{given_meas[column]}[{dict_trad[meas]['unit']}]"
                                        ]
                                    )
                                    / Table[meas][
                                        f"{given_meas[column]}[{dict_trad[meas]['unit']}]"
                                    ]
                                    * 100,
                                    2,
                                )
                            )
                        if plotm_show:
                            ## simu plot measure vs I
                            simu[config].plot(
                                x=f"{filter}I[A]",
                                y=column,
                                ax=ax,
                                ylabel=column.replace(filter, ""),
                                title=column.replace(filter, ""),
                                style=pointstyle_cycler[color_index],
                                markerfacecolor="none",
                                color=color_cycler[color_index],
                            )
                            legend.append(config)
                            color_index += 1

                    if args.show:
                        ax.legend(legend)
                        ax.grid(True)

    for meas in args.measures:  ## go through each measure
        Table[meas] = Table[meas].reset_index(drop=True)
        for column in Table[meas]:
            Table[meas] = Table[meas].rename(
                columns={
                    column: column.replace(filter, "")
                    .replace("_bitter", "")
                    .replace("_helix", "")
                    .replace("[ohm]", "")
                }
            )
        if meas != "R(I)":
            Table[meas] = Table[meas].round(3)
        else:
            Table[meas] = Table[meas].rename(
                columns={column: column.replace("[ohm]", "")}
            )

    return Table, Be


def plot_values(args):
    color_cycler = [
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#F0E442",
    ]
    color_cooling = {
        "mean": color_cycler[0],
        "meanH": color_cycler[1],
        "grad": color_cycler[2],
        "gradH": color_cycler[3],
        "gradHZ": color_cycler[4],
    }
    measuresUnit = {
        "Uw": "Uw [m/s]",
        "DT": "DT [C]",
        "Flux": "Flux [l/s]",
        "HeatCoeff": "Heat Coeff [W/m²/K]",
        "PowerH": "PowerH [W]",
        "statsTH": "statsTH [K]",
        "cf": "cf",
    }
    measuresfeel = {
        "Uw": "Uw",
        "DT": "dTw",
        "Flux": "Flux",
        "HeatCoeff": "hw",
        "PowerH": "Power",
        "statsTH": "statsTH",
        "cf": "cf",
    }
    if args.mtype == "bitter":
        inout = {"Be": "Outer Bitter", "Bi": "Inner Bitter"}
    for meas in args.measures:
        if args.mtype == "bitter":
            figBi, axBi = plt.subplots(figsize=(12, 8))
            figBe, axBe = plt.subplots(figsize=(12, 8))
            figs = {"Be": figBe, "Bi": figBi}
            axs = {"Be": axBe, "Bi": axBi}
            legends = {"Be": [], "Bi": []}
            if args.modeldir and meas in ["cf", "DT", "HeatCoeff", "Uw"]:
                df = pd.read_csv(
                    f"{args.modeldir}/{meas}.measures/values.csv",
                )
                for b in ["Be", "Bi"]:

                    # Rename columns
                    selected_rows = df[df["Slit"].str.contains(b)]
                    selected_rows.index = selected_rows.index.map(
                        lambda x: str(x).replace(
                            f"{measuresfeel[meas]}_M9_{b}_Slit", ""
                        )
                    )
                    # selected_rows.plot(
                    #     x="Slit",
                    #     y=[meas],
                    #     kind="line",
                    #     ax=axs[b],
                    #     color="black",
                    #     linewidth=3,
                    # )
                    selected_rows.plot.bar(
                        alpha=0.7,
                        # xlabel="Cooling Slits",
                        # ylabel=measures[meas],
                        # title=f"{args.filter}{meas} {inout[b]} ",
                        rot=0,
                        ax=axs[b],
                        color=color_cycler[0],
                        position=0,
                        width=0.3,
                        edgecolor="black",
                    )
                    legends[b].append("Simple model")
        else:
            fig, ax = plt.subplots(figsize=(12, 8))
            legend = []
        pos = 0
        # maxpos = len(args.coolings) * len(args.heatcorrelations) * len(args.frictions)
        for cooling in args.coolings:
            heatcorrelations = args.heatcorrelations
            frictions = args.frictions
            if cooling in ["mean", "grad"]:
                heatcorrelations = ["Montgomery"]
                frictions = ["Constant"]

            for heatcorrelation in heatcorrelations:
                for friction in frictions:
                    if (
                        len(args.coolings) == 1
                        and len(heatcorrelations) == 1
                        and len(frictions) == 1
                    ):
                        config = "Feel++ Axi"
                    else:
                        config = f"{cooling}_{heatcorrelation}_{friction}"

                    df = pd.read_csv(
                        f"{args.dir}/{cooling}/{heatcorrelation}/{friction}/np_{args.np}/{args.filter}{meas}.measures/values.csv",
                        index_col=0,
                    )
                    df = df.reindex(index=natsorted(df.index))
                    if args.mtype == "bitter":
                        for b in ["Be", "Bi"]:

                            selected_rows = df[df.index.str.contains(b)]
                            # print(selected_rows)
                            selected_rows.index = selected_rows.index.map(
                                lambda x: str(x).replace(
                                    f"{measuresfeel[meas]}_M9_{b}_Slit", ""
                                )
                            )

                            # if meas == "DT":
                            #     selected_rows["DTT"] = selected_rows["I=33000.0A"]

                            selected_rows.plot.bar(
                                alpha=0.7,
                                xlabel="",
                                ylabel=measuresUnit[meas],
                                title=f"{args.filter}{meas} {inout[b]} ",
                                rot=0,
                                ax=axs[b],
                                color=color_cooling[cooling],  # "r"],
                                position=1,  # pos,
                                width=0.3,
                                edgecolor="black",
                            )
                            legends[b].append(config)
                    else:
                        df.plot.bar(
                            alpha=0.7,
                            xlabel="",
                            ylabel=measuresUnit[meas],
                            title=f"{args.filter}{meas}",
                            rot=45,
                            ax=axs[b],
                            color=color_cooling[cooling],
                            position=pos,
                            width=0.1,
                            edgecolor="black",
                        )
                        legend.append(config)
                    pos += 1  # / maxpos
        if args.mtype == "bitter":
            for b in ["Be", "Bi"]:
                axs[b].legend(legends[b], fontsize=18, loc="lower right")
                axs[b].set_title(f"{args.filter}{meas}  {inout[b]} ", fontsize=20)
                axs[b].grid(True, linestyle="--")
                axs[b].set_ylabel(measuresUnit[meas], fontsize=18)
                axs[b].set_xlabel("Cooling Slits", fontsize=18)
                axs[b].tick_params(axis="x", which="major", labelsize=15)
                axs[b].tick_params(axis="y", which="major", labelsize=15)
        else:
            ax.legend(legend, fontsize=18, loc="best")
            ax.grid(True, linestyle="--")

        directory = ""
        if args.dir:
            directory = f"{args.dir}/"
        if args.save:
            name = f"{args.coolings}_{args.heatcorrelations}_{args.frictions}"
            name = name.replace("[", "").replace("]", "").replace("'", "")
            os.makedirs(f"{directory}{args.filter}res_plot", exist_ok=True)
            if args.mtype == "bitter":
                for b in ["Be", "Bi"]:
                    figs[b].tight_layout()
                    figs[b].savefig(
                        f"{directory}{args.filter}res_plot/{name}_{inout[b]}_{args.filter}{meas}.png"
                    )
            else:
                fig.tight_layout()
                fig.savefig(
                    f"{directory}{args.filter}res_plot/{name}_{args.filter}{meas}.png"
                )


def plot_profiles(args):
    color_cycler = [
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#F0E442",
    ]
    dict_measures = {
        "Feel": {"T": "Celsius", "J": "J", "Vonmises": "vm", "HoopStress": "hs"},
        "SimpleModel": {"T": "temperature[C]", "J": "currentDensity[A/m^2]"},
        "Ansys": {
            "T": "T[C]",
            "J": "Current_Density[A/m2]",
            "Vonmises": "sigma[MPa]",
            # "VonMises": "VonMises[MPa]",
        },
    }
    dict_ylabel = {
        "T": "T [C]",
        "J": "Jth [A/m²]",
        "Vonmises": "Von-Mises Stress [MPa]",
        "HoopStress": "HoopStress [MPa]",
    }

    files = {}
    bitters = []
    for key in args.files:
        files[key] = pd.read_csv(args.files[key])
        bitters.append(key)

    fig2, ax2 = plt.subplots()

    for cooling in args.coolings:
        heatcorrelations = ["Montgomery"]
        frictions = ["Constant"]
        if "H" in cooling:
            heatcorrelations = args.heatcorrelations
            frictions = args.frictions

        for heatcorrelation in heatcorrelations:
            for friction in frictions:
                if args.debug:
                    print(
                        f"#################### {cooling}_{heatcorrelation}_{friction}"
                    )

                for meas in args.measures:
                    fig, ax = plt.subplots(
                        nrows=len(bitters), ncols=len(args.parts), figsize=(14, 9)
                    )

                    axindex1 = 0
                    for i in bitters:
                        # print(i)
                        file = files[i]
                        # print(file)
                        axindex2 = 0
                        for j in args.parts:
                            # print(j)

                            if len(bitters) != 1:
                                if len(args.parts) != 1:
                                    col = ax[axindex1, axindex2]
                                    axindex2 += 1
                                else:
                                    col = ax[axindex1]
                            else:
                                if len(args.parts) != 1:
                                    col = ax[axindex2]
                                    axindex2 += 1
                                else:
                                    col = ax

                            legend = []
                            legend2 = []
                            feel = pd.read_csv(
                                f"{args.dir}/{cooling}/{heatcorrelation}/{friction}/np_1/{i}_{j}_profile.csv"
                            )
                            # x1 = (feel["Points:0"] * 1e3).to_numpy()
                            feel["mm"] = feel["Points:0"] * 1e3
                            if meas == "T":
                                feel["Celsius"] = (
                                    feel["cfpdes.heat.temperature"] - 273.15
                                )
                                feel.plot(
                                    x="mm",
                                    y=["Celsius"],
                                    title="Bitters profile",
                                    style=".",
                                    xlabel="r[mm]",
                                    ylabel="T[C]",
                                    ax=ax2,
                                )
                                legend2.append(f"Feel++_{cooling}_{i}_{j}_T[C]")
                            elif meas == "J":
                                feel["J"] = -feel["cfpdes.expr.Jth"]
                                if not args.hideSimpleModel:
                                    file["core_currentDensity[A/m^2]"] = (
                                        file["core_currentDensity[A/mm^2]"] * 1e6
                                    )
                            elif meas == "Vonmises":
                                feel["vm"] = feel["cfpdes.expr.Vonmises"] * 1e-6
                            elif meas == "HoopStress":
                                feel["hs"] = feel["cfpdes.expr.HoopStress"] * 1e-6

                            if not args.hideSimpleModel and meas in ["T", "J"]:
                                file.plot(
                                    x="r",
                                    y=[f"{j}_{dict_measures['SimpleModel'][meas]}"],
                                    title=f"{i} Bitter {j}_profile",
                                    # xlabel="r[mm]",
                                    # ylabel="T[C]",
                                    ax=col,
                                    linewidth=3,
                                    color=color_cycler[0],
                                )
                                # legend.append("HMFL")
                                legend.append("Simple Model")

                            if not args.hideAxi:
                                feel.plot(
                                    x="mm",
                                    y=[dict_measures["Feel"][meas]],
                                    # title=f"{i} Bitter {j}_profile",
                                    # style=".",
                                    # xlabel="r[mm]",
                                    # ylabel="T[C]",
                                    ax=col,
                                    linewidth=3,
                                    color=color_cycler[1],
                                )

                                legend.append("Feel++ Axi")
                                # legend.append("LNCMI")

                            if args.Ansys:
                                ansys = pd.read_csv(f"{args.Ansys}")
                                if "Center_R[m]" in ansys.columns:
                                    ansys_center = ansys.drop_duplicates(
                                        subset="Center_R[m]"
                                    )
                                    ansys_even = ansys.drop_duplicates(
                                        subset="Even_Rings_R[m]"
                                    )
                                    ansys_odd = ansys.drop_duplicates(
                                        subset="Odd_Rings_R[m]"
                                    )
                                    ansys_center["Center_R[mm]"] = (
                                        ansys_center["Center_R[m]"].copy() * 1e3
                                    )
                                    ansys_even["Even_Rings_R[mm]"] = (
                                        ansys_even["Even_Rings_R[m]"] * 1e3
                                    )
                                    ansys_odd["Odd_Rings_R[mm]"] = (
                                        ansys_odd["Odd_Rings_R[m]"] * 1e3
                                    )
                                else:
                                    ansys_center = ansys.drop_duplicates(
                                        subset="Radial_Position[mm]"
                                    )
                                    ansys_even = ansys.drop_duplicates(
                                        subset="Radial_Position[mm]"
                                    )
                                    ansys_odd = ansys.drop_duplicates(
                                        subset="Radial_Position[mm]"
                                    )
                                    ansys_center["Center_R[mm]"] = ansys_center[
                                        "Radial_Position[mm]"
                                    ].copy()
                                    ansys_even["Even_Rings_R[mm]"] = ansys_even[
                                        "Radial_Position[mm]"
                                    ].copy()
                                    ansys_odd["Odd_Rings_R[mm]"] = ansys_odd[
                                        "Radial_Position[mm]"
                                    ].copy()

                                prefix = ""
                                if meas == "Vonmises":
                                    if (
                                        "laplace" in args.dir
                                        and "dilatation" not in args.dir
                                    ):
                                        prefix = "LorentzForces_"
                                    elif (
                                        "dilatation" in args.dir
                                        and "laplace" not in args.dir
                                    ):
                                        prefix = "ThermalLoad_"

                                ansys_center.plot(
                                    x="Center_R[mm]",
                                    y=f"{prefix}Center_{dict_measures['Ansys'][meas]}",
                                    style="--",
                                    ax=col,
                                    # linewidth=3,
                                    color=color_cycler[4],
                                )
                                legend.append("Ansys Center")

                                ansys_odd.plot(
                                    x="Odd_Rings_R[mm]",
                                    y=f"{prefix}Odd_Rings_{dict_measures['Ansys'][meas]}",
                                    style="--",
                                    ax=col,
                                    # linewidth=3,
                                    color=color_cycler[2],  # "green",
                                )
                                legend.append("Ansys Odd Rings")
                                legend = plot_greySlits(
                                    ansys_odd,
                                    "Odd_Rings_R[mm]",
                                    f"{prefix}Odd_Rings_{dict_measures['Ansys'][meas]}",
                                    col,
                                    legend,
                                )

                                ansys_even.plot(
                                    x="Even_Rings_R[mm]",
                                    y=f"{prefix}Even_Rings_{dict_measures['Ansys'][meas]}",
                                    style="--",
                                    ax=col,
                                    # linewidth=3,
                                    color=color_cycler[3],  # "yellow",
                                )
                                legend.append("Ansys Even Rings")
                                legend = plot_greySlits(
                                    ansys_even,
                                    "Even_Rings_R[mm]",
                                    f"{prefix}Even_Rings_{dict_measures['Ansys'][meas]}",
                                    col,
                                    legend,
                                )

                            if args.boxplotsdir:
                                if args.show2Dlines and meas in ["T", "J", "Vonmises"]:
                                    center = pd.read_csv(
                                        f"{args.boxplotsdir}/{cooling}/{heatcorrelation}/{friction}/np_1/center_profile.csv"
                                    )
                                    center["mm"] = center["Points:0"] * 1e3
                                    border = pd.read_csv(
                                        f"{args.boxplotsdir}/{cooling}/{heatcorrelation}/{friction}/np_1/border_profile.csv"
                                    )

                                    border["mm"] = (
                                        np.sqrt(
                                            border["Points:0"] ** 2
                                            + border["Points:1"] ** 2
                                        )
                                        * 1e3
                                    )
                                    odd = None
                                    try:
                                        odd = pd.read_csv(
                                            f"{args.boxplotsdir}/{cooling}/{heatcorrelation}/{friction}/np_1/odd_profile.csv"
                                        )
                                    except:
                                        odd = None

                                    if type(odd) != None:
                                        odd["mm"] = (
                                            np.sqrt(
                                                odd["Points:0"] ** 2
                                                + odd["Points:1"] ** 2
                                            )
                                            * 1e3
                                        )

                                    if meas == "T":
                                        center["Celsius"] = (
                                            center[f"{args.toolbox}.heat.temperature"]
                                            - 273.15
                                        )
                                        border["Celsius"] = (
                                            border[f"{args.toolbox}.heat.temperature"]
                                            - 273.15
                                        )
                                        if type(odd) != None:
                                            odd["Celsius"] = (
                                                odd[f"{args.toolbox}.heat.temperature"]
                                                - 273.15
                                            )
                                    if meas == "J":
                                        center["J"] = np.sqrt(
                                            center["cfpdes.expr.J:0"] ** 2
                                            + center["cfpdes.expr.J:1"] ** 2
                                        )
                                        border["J"] = np.sqrt(
                                            border["cfpdes.expr.J:0"] ** 2
                                            + border["cfpdes.expr.J:1"] ** 2
                                        )
                                        if type(odd) != None:
                                            odd["J"] = np.sqrt(
                                                odd["cfpdes.expr.J:0"] ** 2
                                                + odd["cfpdes.expr.J:1"] ** 2
                                            )
                                    if meas == "Vonmises":
                                        center["vm"] = (
                                            center[f"{args.toolbox}.expr.Vonmises"]
                                            * 1e-6
                                        )
                                        border["vm"] = (
                                            border[f"{args.toolbox}.expr.Vonmises"]
                                            * 1e-6
                                        )
                                        if type(odd) != None:
                                            odd["vm"] = (
                                                odd[f"{args.toolbox}.expr.Vonmises"]
                                                * 1e-6
                                            )

                                    center.plot(
                                        x="mm",
                                        y=dict_measures["Feel"][meas],
                                        style="-.",
                                        ax=col,
                                        color=color_cycler[2],
                                    )
                                    legend.append("Feel++ 2D Center")

                                    border.plot(
                                        x="mm",
                                        y=dict_measures["Feel"][meas],
                                        style="-.",
                                        ax=col,
                                        color=color_cycler[3],
                                    )
                                    # legend.append("Feel++ 2D Even Rings")
                                    if type(odd) != None:
                                        legend.append("Feel++ 2D Even Rings")
                                        odd.plot(
                                            x="mm",
                                            y=dict_measures["Feel"][meas],
                                            style="-.",
                                            ax=col,
                                            color=color_cycler[1],
                                        )
                                        legend.append("Feel++ 2D Odd Rings")
                                    else:
                                        legend.append("Feel++ 2D Border")

                                    if args.Ansys:
                                        figA, axA = plt.subplots(figsize=(15, 8))
                                        legendA = []

                                        # ansys_center.plot(
                                        #     x="Center_R[mm]",
                                        #     y=f"{prefix}Center_{dict_measures['Ansys'][meas]}",
                                        #     style="--",
                                        #     linewidth=2,
                                        #     ax=axA,
                                        #     color=color_cycler[0],
                                        # )
                                        # legendA.append("Ansys Center")
                                        # center.plot(
                                        #     x="mm",
                                        #     y=dict_measures["Feel"][meas],
                                        #     style="-",
                                        #     linewidth=2,
                                        #     ax=axA,
                                        #     color=color_cycler[1],
                                        #     alpha=0.7,
                                        # )
                                        # legendA.append("Feel++ 2D Center")
                                        # legendA = plot_greySlits(
                                        #     ansys_center,
                                        #     "Center_R[mm]",
                                        #     f"{prefix}Center_{dict_measures['Ansys'][meas]}",
                                        #     axA,
                                        #     legendA,
                                        # )
                                        # err = compare(
                                        #     ansys_center,
                                        #     "Center_R[mm]",
                                        #     f"{prefix}Center_{dict_measures['Ansys'][meas]}",
                                        #     center,
                                        #     "mm",
                                        #     dict_measures["Feel"][meas],
                                        # )
                                        # print(
                                        #     f"Center Feel++2D vs Ansys: max error={np.round(err[2],3)} at r={np.round(err[3],2)}"
                                        # )
                                        ansys_odd.plot(
                                            x="Odd_Rings_R[mm]",
                                            y=f"{prefix}Odd_Rings_{dict_measures['Ansys'][meas]}",
                                            style="--",
                                            linewidth=2,
                                            ax=axA,
                                            color=color_cycler[2],
                                        )
                                        legendA.append("Ansys Odd Rings")
                                        odd.plot(
                                            x="mm",
                                            y=dict_measures["Feel"][meas],
                                            style="-",
                                            linewidth=2,
                                            ax=axA,
                                            color=color_cycler[3],
                                            alpha=0.7,
                                        )
                                        legendA.append("Feel++ 2D Odd Rings")
                                        legendA = plot_greySlits(
                                            ansys_odd,
                                            "Odd_Rings_R[mm]",
                                            f"{prefix}Odd_Rings_{dict_measures['Ansys'][meas]}",
                                            axA,
                                            legendA,
                                        )
                                        err = compare(
                                            ansys_odd,
                                            "Odd_Rings_R[mm]",
                                            f"{prefix}Odd_Rings_{dict_measures['Ansys'][meas]}",
                                            odd,
                                            "mm",
                                            dict_measures["Feel"][meas],
                                        )
                                        print(
                                            f"Odd Rings Feel++2D vs Ansys: max error={np.round(err[2],3)} at r={err[3]}"
                                        )
                                        ansys_even.plot(
                                            x="Even_Rings_R[mm]",
                                            y=f"{prefix}Even_Rings_{dict_measures['Ansys'][meas]}",
                                            style="--",
                                            linewidth=2,
                                            ax=axA,
                                            color=color_cycler[4],
                                        )
                                        legendA.append("Ansys Even Rings")
                                        border.plot(
                                            x="mm",
                                            y=dict_measures["Feel"][meas],
                                            style="-",
                                            linewidth=2,
                                            ax=axA,
                                            color=color_cycler[5],
                                            alpha=0.7,
                                        )
                                        legendA.append("Feel++ 2D Even Rings")
                                        legendA = plot_greySlits(
                                            ansys_even,
                                            "Even_Rings_R[mm]",
                                            f"{prefix}Even_Rings_{dict_measures['Ansys'][meas]}",
                                            axA,
                                            legendA,
                                        )
                                        err = compare(
                                            ansys_even,
                                            "Even_Rings_R[mm]",
                                            f"{prefix}Even_Rings_{dict_measures['Ansys'][meas]}",
                                            border,
                                            "mm",
                                            dict_measures["Feel"][meas],
                                        )
                                        print(
                                            f"Even Rings Feel++2D vs Ansys: max error={np.round(err[2],3)} at r={err[3]}"
                                        )

                                        axA.legend(legendA, fontsize=18, loc="best")
                                        axA.set_xlabel("r[mm]", fontsize=18)
                                        axA.set_ylabel(dict_ylabel[meas], fontsize=18)
                                        axA.tick_params(
                                            axis="both", which="major", labelsize=15
                                        )
                                        axA.grid(linestyle="--")
                                        axA.margins(x=0)
                                        if args.save:
                                            directory = f"{args.boxplotsdir}/"
                                            os.makedirs(
                                                f"{directory}res_profile", exist_ok=True
                                            )
                                            figA.tight_layout()
                                            figA.savefig(
                                                f"{directory}res_profile/{cooling}_{heatcorrelation}_{friction}_2D_Ansys_{meas}_profiles.png"
                                            )

                                    if args.boxplotsdir2:
                                        center2 = pd.read_csv(
                                            f"{args.boxplotsdir2}/{cooling}/{heatcorrelation}/{friction}/np_1/center_profile.csv"
                                        )
                                        center2["mm"] = center2["Points:0"] * 1e3
                                        border2 = pd.read_csv(
                                            f"{args.boxplotsdir2}/{cooling}/{heatcorrelation}/{friction}/np_1/border_profile.csv"
                                        )

                                        border2["mm"] = (
                                            np.sqrt(
                                                border2["Points:0"] ** 2
                                                + border2["Points:1"] ** 2
                                            )
                                            * 1e3
                                        )
                                        odd2 = None
                                        try:
                                            odd2 = pd.read_csv(
                                                f"{args.boxplotsdir2}/{cooling}/{heatcorrelation}/{friction}/np_1/odd_profile.csv"
                                            )
                                        except:
                                            odd2 = None

                                        if type(odd2) != None:
                                            odd2["mm"] = (
                                                np.sqrt(
                                                    odd2["Points:0"] ** 2
                                                    + odd2["Points:1"] ** 2
                                                )
                                                * 1e3
                                            )

                                        if meas == "T":
                                            center2["Celsius"] = (
                                                center2[
                                                    f"{args.toolbox}.heat.temperature"
                                                ]
                                                - 273.15
                                            )
                                            border2["Celsius"] = (
                                                border2[
                                                    f"{args.toolbox}.heat.temperature"
                                                ]
                                                - 273.15
                                            )
                                            if type(odd2) != None:
                                                odd2["Celsius"] = (
                                                    odd2[
                                                        f"{args.toolbox}.heat.temperature"
                                                    ]
                                                    - 273.15
                                                )
                                        if meas == "J":
                                            center2["J"] = np.sqrt(
                                                center2["cfpdes.expr.J:0"] ** 2
                                                + center2["cfpdes.expr.J:1"] ** 2
                                            )
                                            border2["J"] = np.sqrt(
                                                border2["cfpdes.expr.J:0"] ** 2
                                                + border2["cfpdes.expr.J:1"] ** 2
                                            )
                                            if type(odd2) != None:
                                                odd2["J"] = np.sqrt(
                                                    odd2["cfpdes.expr.J:0"] ** 2
                                                    + odd2["cfpdes.expr.J:1"] ** 2
                                                )
                                        if meas == "Vonmises":
                                            center2["vm"] = (
                                                center2[f"{args.toolbox}.expr.Vonmises"]
                                                * 1e-6
                                            )
                                            border2["vm"] = (
                                                border2[f"{args.toolbox}.expr.Vonmises"]
                                                * 1e-6
                                            )
                                            if type(odd2) != None:
                                                odd2["vm"] = (
                                                    odd2[
                                                        f"{args.toolbox}.expr.Vonmises"
                                                    ]
                                                    * 1e-6
                                                )

                                        figA, axA = plt.subplots(figsize=(15, 8))
                                        legendA = []

                                        odd2.plot(
                                            x="mm",
                                            y=dict_measures["Feel"][meas],
                                            style="--",
                                            linewidth=2,
                                            ax=axA,
                                            color=color_cycler[3],
                                        )
                                        legendA.append(
                                            "Feel++ 2D Odd Rings w fixed Tierods"
                                        )
                                        odd.plot(
                                            x="mm",
                                            y=dict_measures["Feel"][meas],
                                            style="-",
                                            linewidth=2,
                                            ax=axA,
                                            color=color_cycler[3],
                                            alpha=0.7,
                                        )
                                        legendA.append("Feel++ 2D Odd Rings")
                                        legendA = plot_greySlits(
                                            odd2,
                                            "mm",
                                            dict_measures["Feel"][meas],
                                            axA,
                                            legendA,
                                        )

                                        border2.plot(
                                            x="mm",
                                            y=dict_measures["Feel"][meas],
                                            style="--",
                                            linewidth=2,
                                            ax=axA,
                                            color=color_cycler[5],
                                        )
                                        legendA.append(
                                            "Feel++ 2D Even Rings w fixed Tierods"
                                        )
                                        border.plot(
                                            x="mm",
                                            y=dict_measures["Feel"][meas],
                                            style="-",
                                            linewidth=2,
                                            ax=axA,
                                            color=color_cycler[5],
                                            alpha=0.7,
                                        )
                                        legendA.append("Feel++ 2D Even Rings")
                                        legendA = plot_greySlits(
                                            border2,
                                            "mm",
                                            dict_measures["Feel"][meas],
                                            axA,
                                            legendA,
                                        )

                                        axA.legend(legendA, fontsize=18, loc="best")
                                        axA.set_xlabel("r[mm]", fontsize=18)
                                        axA.set_ylabel(dict_ylabel[meas], fontsize=18)
                                        axA.tick_params(
                                            axis="both", which="major", labelsize=15
                                        )
                                        axA.grid(linestyle="--")
                                        axA.margins(x=0)
                                        if args.save:
                                            directory = f"{args.boxplotsdir}/"
                                            os.makedirs(
                                                f"{directory}res_profile", exist_ok=True
                                            )
                                            figA.tight_layout()
                                            figA.savefig(
                                                f"{directory}res_profile/{cooling}_{heatcorrelation}_{friction}_2D_vs_{meas}_profiles.png"
                                            )

                                if meas == "T":
                                    boxplots = pd.read_csv(
                                        f"{args.boxplotsdir}/{cooling}/{heatcorrelation}/{friction}/np_1/r_profiles.csv"
                                    )

                                    boxplots[f"{args.toolbox}.heat.temperature"] = (
                                        boxplots[f"{args.toolbox}.heat.temperature"]
                                        - 273.15
                                    )

                                    rs = np.unique(boxplots["r"])

                                    boxplots.boxplot(
                                        column=[f"{args.toolbox}.heat.temperature"],
                                        by=["r"],
                                        positions=rs,
                                        widths=3,
                                        boxprops=dict(
                                            linestyle="-", linewidth=1.5, color="r"
                                        ),
                                        flierprops=dict(
                                            linestyle="-", linewidth=1.5, color="r"
                                        ),
                                        medianprops=dict(
                                            linestyle="-", linewidth=1.5, color="r"
                                        ),
                                        whiskerprops=dict(
                                            linestyle="-", linewidth=1.5, color="r"
                                        ),
                                        capprops=dict(
                                            linestyle="-", linewidth=1.5, color="r"
                                        ),
                                        showfliers=False,
                                        ax=col,
                                    )
                                    if "BitterC" in args.dir:
                                        col.set_xticks(
                                            [
                                                120,
                                                125,
                                                130,
                                                135,
                                                140,
                                                145,
                                                150,
                                                155,
                                                160,
                                                165,
                                                170,
                                                175,
                                                180,
                                            ],
                                            [
                                                "120",
                                                "125",
                                                "130",
                                                "135",
                                                "140",
                                                "145",
                                                "150",
                                                "155",
                                                "160",
                                                "165",
                                                "170",
                                                "175",
                                                "180",
                                            ],
                                        )
                                    else:
                                        col.set_xticks(
                                            np.arange(
                                                int(rs[0] * 1e-1) * 10, rs[-1] + 10, 20
                                            ),
                                            np.arange(
                                                int(rs[0] * 1e-1) * 10, rs[-1] + 10, 20
                                            ),
                                        )
                                    legend.append("Feel++ 2D")

                        col.legend(legend, fontsize=18, loc="best")
                        col.set_xlabel("r[mm]", fontsize=18)
                        col.set_ylabel(dict_ylabel[meas], fontsize=18)
                        col.tick_params(axis="both", which="major", labelsize=15)
                        col.margins(x=0)
                        # col.set_ylim([29, 46])
                        col.grid(linestyle="--")
                        col.set_title(f"{i} Bitter {j} {meas} profile")

                        axindex1 += 1
                    if args.save:
                        directory = f"{args.dir}/"
                        if args.boxplotsdir:
                            directory = f"{args.boxplotsdir}/"
                        os.makedirs(f"{directory}res_profile", exist_ok=True)

                        fig.suptitle(
                            f"{cooling} {heatcorrelation} {friction} {meas} profiles",
                            fontsize=16,
                        )
                        fig.tight_layout()
                        fig.savefig(
                            f"{directory}res_profile/{cooling}_{heatcorrelation}_{friction}_{meas}_profiles.png"
                        )

    if args.show and "T" in args.measures:
        ax2.grid()
        ax2.legend(legend2)
    else:
        plt.close(fig2)

    return 0


def main():

    parser = options("", "")
    args = parser.parse_args()

    if args.debug:
        print(f"args={args}")

    if args.command == "plot_measures":
        plotm_show = args.show

        Table, Be = plot_measures(args, plotm_show)

        directory = ""
        if args.dir:
            directory = f"{args.dir}/"

        for meas in args.measures:  ## go through each measure
            ## save measure table as txt
            if args.save:
                name = (
                    f"{args.coolings}_{args.heatcorrelations}_{args.frictions}_{meas}"
                )
                name = name.replace("[", "").replace("]", "").replace("'", "")
                os.makedirs(f"{directory}{args.filter}res_plot", exist_ok=True)
                with open(
                    f"{directory}{args.filter}res_plot/{args.filter}{name}.txt",
                    "w",
                ) as f:
                    f.write(
                        tabulate(
                            Table[meas],
                            headers="keys",
                            tablefmt="pretty",
                            showindex=False,
                        )
                    )
                Table[meas].T.to_csv(
                    f"{directory}{args.filter}res_plot/{args.filter}{name}.csv",
                    index=True,
                )
            else:
                print("\n")
                ## print all measure tables
                print(
                    tabulate(
                        Table[meas],
                        headers="keys",
                        tablefmt="fancy_grid",
                        showindex=False,
                    )
                )

    if args.command == "plot_values":
        plot_values(args)

    if args.command == "plot_profiles":
        if args.mtype == "bitter":
            plot_profiles(args)
        else:
            print("plot_profiles is only available for bitters")

    if args.show:
        plt.show()

    return 0


if __name__ == "__main__":
    sys.exit(main())
