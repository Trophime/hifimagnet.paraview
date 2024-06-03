import pandas as pd
import os
import matplotlib.pyplot as plt

from typing import List

from paraview.simple import (
    PlotOverLine,
    CreateWriter,
)

from ..method import convert_data


def plotOr(
    input,
    r: list[float],
    z: float,
    fieldunits: dict,
    basedir: str,
    show: bool = True,
    printed: bool = True,
):
    os.makedirs(f"{basedir}/plots", exist_ok=True)

    [r0, r1] = r

    plotOverLine = PlotOverLine(registrationName="Oz", Input=input, Source="Line")
    # get params list
    if not printed:
        for prop in plotOverLine.ListProperties():
            print(f"plotOverLine': {prop}={plotOverLine.GetPropertyValue(prop)}")

    # init the 'Line' selected for 'Source'
    plotOverLine.Source.Point1 = [r0, z, 0]
    plotOverLine.Source.Point2 = [r1, z, 0]

    filename = f"{basedir}/plots/r0={r0}m-r1={r1}m-z={z}m.csv"
    export = CreateWriter(filename, proxy=plotOverLine)
    if not printed:
        for prop in export.ListProperties():
            print(f"export: {prop}={export.GetPropertyValue(prop)}")
    export.UpdateVTKObjects()  # is it needed?
    export.UpdatePipeline()

    # plot with matplotlib
    def plotOrField(
        file,
        key: str,
        z: float,
        fieldunits: dict,
        basedir: str,
        ax=None,
        show: bool = True,
    ):

        print(f"plotOrField: file={file}, key={key}", flush=True)
        keyinfo = key.replace("_Magnitude", "").split(".")
        # print(f"keyinfo={keyinfo}", flush=True)
        if len(keyinfo) == 2:
            (physic, fieldname) = keyinfo
        elif len(keyinfo) == 3:
            (toolbox, physic, fieldname) = keyinfo
        else:
            raise RuntimeError(f"{key}: cannot get keyinfo as splitted char")
        # print(f"physic={physic}, fieldname={fieldname}", flush=True)
        # print(f'fieldunits[fieldname]={fieldunits[fieldname]}"', flush=True)
        symbol = fieldunits[fieldname]["Symbol"]
        [in_unit, out_unit] = fieldunits[fieldname]["Units"]
        # print(f"in_units={in_unit}, out_units={out_unit}", flush=True)

        if ax is None:
            ax = plt.gca()

        # see vonmises-vs-theta.py and/or vonmises-vs-theta-plot-savedata.py
        csv = pd.read_csv(file)
        keycsv = csv[["arc_length", key]]

        # rename columns
        keycsv.rename(columns={"arc_length": "r"}, inplace=True)
        print(f"new keys: {keycsv.columns.values.tolist()}")

        # rescale columns to plot
        units = {fieldname: fieldunits[fieldname]["Units"]}
        values = keycsv[key].to_list()
        out_values = convert_data(units, values, fieldname)
        keycsv = keycsv.assign(key=[f"{val:.2f}" for val in out_values])

        keycsv.plot(x="r", y=key, marker="o", grid=True, ax=ax)

        plt.xlabel("r [m]")
        plt.ylabel(rf"{symbol} [{out_unit:~P}]")

        # ax.yaxis.set_major_locator(MaxNLocator(10))
        plt.title(f"{key}: z={z} ")

        if show:
            plt.show()
        else:
            plt.savefig(f"{basedir}/plots/{key}-vs-r-z={z}.png", dpi=300)
        plt.close()

    # requirements: create PointData from CellData
    for field in input.PointData:
        plotOrField(
            filename,
            field,
            z,
            fieldunits,
            basedir,
            ax=None,
            show=show,
        )
    pass


def plotOz(
    input,
    r: float,
    z: list[float],
    fieldunits: dict,
    basedir: str,
    show: bool = True,
    printed: bool = True,
):
    os.makedirs(f"{basedir}/plots", exist_ok=True)

    [z0, z1] = z

    plotOverLine = PlotOverLine(registrationName="Oz", Input=input, Source="Line")

    # init the 'Line' selected for 'Source'
    plotOverLine.Source.Point1 = [r, z0, 0]
    plotOverLine.Source.Point2 = [r, z1, 0]

    filename = f"{basedir}/plots/r={r}m-z0={z0}m-z1={z1}m.csv"
    export = CreateWriter(filename, proxy=plotOverLine)
    if not printed:
        for prop in export.ListProperties():
            print(f"export: {prop}={export.GetPropertyValue(prop)}")
    export.UpdateVTKObjects()  # is it needed?
    export.UpdatePipeline()

    # plot with matplotlib
    def plotOzField(
        file,
        key: str,
        z: float,
        fieldunits: dict,
        basedir: str,
        ax=None,
        show: bool = True,
    ):

        print(f"plotOrField: file={file}, key={key}", flush=True)
        keyinfo = key.replace("_Magnitude", "").split(".")
        # print(f"keyinfo={keyinfo}", flush=True)
        if len(keyinfo) == 2:
            (physic, fieldname) = keyinfo
        elif len(keyinfo) == 3:
            (toolbox, physic, fieldname) = keyinfo
        else:
            raise RuntimeError(f"{key}: cannot get keyinfo as splitted char")
        symbol = fieldunits[fieldname]["Symbol"]
        [in_unit, out_unit] = fieldunits[fieldname]["Units"]

        if ax is None:
            ax = plt.gca()

        # see vonmises-vs-theta.py and/or vonmises-vs-theta-plot-savedata.py
        csv = pd.read_csv(file)
        keycsv = csv[["arc_length", key]]

        # rename columns
        keycsv.rename(columns={"arc_length": "z"}, inplace=True)
        print(f"new keys: {keycsv.columns.values.tolist()}")

        # rescale columns to plot
        units = {fieldname: fieldunits[fieldname]["Units"]}
        values = keycsv[key].to_list()
        out_values = convert_data(units, values, fieldname)
        keycsv = keycsv.assign(key=[f"{val:.2f}" for val in out_values])

        keycsv.plot(x="z", y=key, marker="o", grid=True, ax=ax)

        plt.xlabel("z [m]")
        plt.ylabel(rf"{symbol} [{out_unit:~P}]")

        # ax.yaxis.set_major_locator(MaxNLocator(10))
        plt.title(f"{key}: r={r} ")

        if show:
            plt.show()
        else:
            plt.savefig(f"{basedir}/plots/{key}-vs-z-r={r}.png", dpi=300)
        plt.close()

    # requirements: create PointData from CellData
    for field in input.PointData:
        plotOzField(
            filename,
            field,
            r,
            fieldunits,
            basedir,
            ax=None,
            show=show,
        )
    pass
