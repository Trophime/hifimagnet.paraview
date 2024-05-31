import pandas as pd
import os
import gc

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from numpy import pi, arctan2
from math import pi, cos, sin


from typing import List

from paraview.simple import (
    PlotOnIntersectionCurves,
    CellDatatoPointData,
    PlotOverLine,
    Delete,
    CreateWriter,
    SetActiveSource,
)

from ..method import convert_data, resultinfo
from ..view import makeclip, makecylinderslice


def plotOr(
    input,
    r: list[float],
    theta: float,
    fieldunits: dict,
    basedir: str,
    show: bool = True,
    printed: bool = True,
):

    os.makedirs(f"{basedir}/plots", exist_ok=True)

    [r0, r1] = r
    radian = theta * pi / 180.0

    plotOverLine = PlotOverLine(registrationName="Oz", Input=input, Source="Line")

    # init the 'Line' selected for 'Source'
    plotOverLine.Source.Point1 = [r0 * cos(radian), r0 * sin(radian), 0]
    plotOverLine.Source.Point2 = [r1 * cos(radian), r1 * sin(radian), 0]

    filename = f"{basedir}/plots/r0={r0}m-r1={r1}m-theta={theta}deg.csv"
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
        theta: float,
        fieldunits: dict,
        basedir: str,
        ax=None,
        show: bool = True,
    ):

        print(f"plotOrField: file={file}, key={key}", flush=True)
        keyinfo = key.split(".")
        print(f"keyinfo={keyinfo}", flush=True)
        if len(keyinfo) == 2:
            (physic, fieldname) = keyinfo
        elif len(keyinfo) == 3:
            (toolbox, physic, fieldname) = keyinfo
        else:
            raise RuntimeError(f"{key}: cannot get keyinfo as splitted char")
        # print(f"physic={physic}, fieldname={fieldname}", flush=True)
        # print(f'fieldunits[fieldname]={fieldunits[fieldname]}"', flush=True)
        symbol = fieldunits[fieldname]["Symbol"]
        msymbol = symbol
        if "mSymbol" in fieldunits[fieldname]:
            msymbol = fieldunits[fieldname]["mSymbol"]
        [in_unit, out_unit] = fieldunits[fieldname]["Units"]
        # print(f"in_units={in_unit}, out_units={out_unit}", flush=True)

        if ax is None:
            ax = plt.gca()

        # see vonmises-vs-theta.py and/or vonmises-vs-theta-plot-savedata.py
        csv = pd.read_csv(file)
        keycsv = csv[["arc_length", key]]

        # rename columns
        keycsv.rename(columns={"arc_length": "r"}, inplace=True)
        # print(f"new keys: {keycsv.columns.values.tolist()}")

        # rescale columns to plot
        units = {fieldname: fieldunits[fieldname]["Units"]}
        values = keycsv[key].to_list()
        out_values = convert_data(units, values, fieldname)
        ndf = {key: [val for val in out_values]}

        keycsv[key] = ndf[key]
        keycsv.plot(x="r", y=key, marker="o", grid=True, ax=ax)

        plt.xlabel("r [m]")
        plt.ylabel(rf"{msymbol} [{out_unit:~P}]")

        # ax.yaxis.set_major_locator(MaxNLocator(10))
        plt.title(f"{key}: theta={theta} deg")

        if show:
            plt.show()
        else:
            plt.savefig(f"{basedir}/plots/{key}-vs-r-theta={theta}deg.png", dpi=300)
        plt.close()
        keycsv.to_csv(f"{basedir}/plots/{key}-vs-r-theta={theta}deg.csv")
        pass

    # requirements: create PointData from CellData
    for field in input.PointData:
        plotOrField(
            filename,
            field,
            theta,
            fieldunits,
            basedir,
            ax=None,
            show=show,
        )

    # remove temporary csv files
    os.remove(filename)

    Delete(plotOverLine)
    del plotOverLine
    pass


def plotTheta(
    input,
    r: float,
    fieldunits: dict,
    ignored_keys: List[str],
    basedir: str,
    show: bool = True,
    printed: bool = True,
    verbose: bool = False,
):
    """
    for theta, need to apply CellDataToPointData filter
    """
    os.makedirs(f"{basedir}/plots", exist_ok=True)

    print(f"plotTheta: r={r}", flush=True)
    cellDatatoPointData1 = CellDatatoPointData(
        registrationName="CellDatatoPointData", Input=input
    )

    # create clip with plane (howto give a color for each clip)
    print("cellDatatoPointDatalip up and down", flush=True)
    clip_down = makeclip(cellDatatoPointData1, "clip_down", invert=False)
    clip_up = makeclip(cellDatatoPointData1, "clip_up", invert=True)

    files = []
    for i, clip in enumerate([clip_down, clip_up]):
        slice = makecylinderslice(clip, f"slice{i}", r)
        SetActiveSource(slice)

        plotOnIntersectionCurve = PlotOnIntersectionCurves(
            registrationName=f"PlotOnIntersectionCurves{i}", Input=slice
        )
        # get params list
        if not printed:
            for prop in plotOnIntersectionCurve.ListProperties():
                print(
                    f"plotOnIntersectionCurve: {prop}={plotOnIntersectionCurve.GetPropertyValue(prop)}",
                    flush=True,
                )

        plotOnIntersectionCurve.SliceType = "Plane"
        # Properties modified on plotOnIntersectionCurves1.SliceType
        plotOnIntersectionCurve.SliceType.Origin = [0.0, 0.0, 0]
        plotOnIntersectionCurve.SliceType.Normal = [0.0, 0.0, 1.0]

        # plotOnIntersectionCurves.append(plotOnIntersectionCurve)
        export = CreateWriter(
            f"{basedir}/plots/r={r}m-{i}.csv", proxy=plotOnIntersectionCurve
        )
        if not printed:
            for prop in export.ListProperties():
                print(f"export: {prop}={export.GetPropertyValue(prop)}", flush=True)
        export.UpdateVTKObjects()  # is it needed?
        export.UpdatePipeline()
        files.append(f"{basedir}/plots/r={r}m-{i}.csv")

    # plot with matplotlib
    def plotThetaField(
        files,
        key: str,
        r: float,
        fieldunits: dict,
        basedir: str,
        ax=None,
        show: bool = True,
    ):

        print(f"plotThetaField: files={files}, key={key}, show={show}", flush=True)
        keyinfo = key.split(".")
        print(f"keyinfo={keyinfo}", flush=True)
        if len(keyinfo) == 2:
            (physic, fieldname) = keyinfo
        elif len(keyinfo) == 3:
            (toolbox, physic, fieldname) = keyinfo
        else:
            raise RuntimeError(f"{key}: cannot get keyinfo as splitted char")
        symbol = fieldunits[fieldname]["Symbol"]
        msymbol = symbol
        if "mSymbol" in fieldunits[fieldname]:
            msymbol = fieldunits[fieldname]["mSymbol"]
        [in_unit, out_unit] = fieldunits[fieldname]["Units"]

        if ax is None:
            ax = plt.gca()

        # see vonmises-vs-theta.py and/or vonmises-vs-theta-plot-savedata.py
        keycsv_dfs = []
        for file in files:
            csv = pd.read_csv(file)
            keycsv = csv[["Points:0", "Points:1", "arc_length", key]]

            # rename columns
            keycsv.rename(columns={"Points:0": "x", "Points:1": "y"}, inplace=True)
            # print(f"new keys: {keycsv.columns.values.tolist()}", flush=True)

            # rescale columns to plot
            theta = arctan2(keycsv["y"].to_numpy(), keycsv["x"].to_numpy()) * 180 / pi
            # print(f'theta (type={type(theta)}, nelem={theta.shape}): {theta}', flush=True)
            keycsv["theta"] = keycsv.apply(
                lambda row: arctan2(row.y, row.x) * 180 / pi,
                axis=1,
            )

            # Sort the rows of dataframe by 'Name' column
            keycsv = keycsv.sort_values(by="theta")
            # # drop last line
            # print(f"drop first line from {file}", flush=True)
            # keycsv.drop(index=keycsv.index[-1], axis=0, inplace=True)

            units = {fieldname: fieldunits[fieldname]["Units"]}
            values = keycsv[key].to_list()
            out_values = convert_data(units, values, fieldname)
            ndf = {key: [val for val in out_values]}

            keycsv[key] = ndf[key]
            keycsv_dfs.append(keycsv)

            del ndf

        df = pd.concat(keycsv_dfs)
        df = df.sort_values(by="theta")
        r_units = {"coord": fieldunits["coord"]["Units"]}
        mm = f'{fieldunits["coord"]["Units"][1]:~P}'
        r_mm = convert_data(r_units, r, "coord")
        df.to_csv(f"{basedir}/plots/{key}-vs-theta-r={r_mm}{mm}.csv")
        # print(f"df keys: {df.columns.values.tolist()}", flush=True)
        # assert key in df.columns.values.tolist(), f"{key} not in df_keys"
        # print(f"df={df}", flush=True)
        df.plot(x="theta", y=key, marker="o", grid=True, ax=ax)

        plt.xlabel(r"$\theta$ [deg]")
        plt.ylabel(rf"{msymbol} [{out_unit:~P}]")

        # x range
        x0 = -180
        x1 = 180
        plt.xlim(x0, x1)

        ax.yaxis.set_major_locator(MaxNLocator(10))
        plt.title(f"{key}: r={r_mm} {mm}")

        if show:
            plt.show()
        else:
            plt.savefig(f"{basedir}/plots/{key}-vs-theta-r={r_mm}{mm}.png", dpi=300)
        plt.close()

        print(f"{key} stats on r={r_mm}{mm}", flush=True)
        print(f"{df[key].describe()}", flush=True)
        pass

    # requirements: create PointData from CellData
    datadict = resultinfo(cellDatatoPointData1, ignored_keys)
    for field in datadict["PointData"]["Arrays"]:
        if not field in ignored_keys:
            kdata = datadict["PointData"]["Arrays"][field]
            Components = kdata["Components"]
            print(f"plotThetaField for {field} - components={Components}", flush=True)
            if Components > 1:
                for i in range(Components):
                    print(f"plotThetaField for {field}:{i} skipped", flush=True)
                # for i in range(Components):
                #     plotThetaField(files, f"{field}:{i}", r, z, ax=None, show=show)
            else:
                plotThetaField(files, field, r, fieldunits, basedir, ax=None, show=show)

    # remove temporary csv files
    for file in files:
        os.remove(file)

    Delete(cellDatatoPointData1)
    del cellDatatoPointData1

    # Force a garbage collection
    collected = gc.collect()
    if verbose:
        print(f"Garbage collector: collected {collected} objects.", flush=True)
