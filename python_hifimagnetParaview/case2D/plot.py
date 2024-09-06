import pandas as pd
import os
import gc

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from numpy import pi, arctan2
from math import pi, cos, sin

from paraview.simple import (
    CellDatatoPointData,
    PlotOverLine,
    Delete,
    CreateWriter,
    SetActiveSource,
)

from ..method import convert_data, resultinfo, showplot, plot_greySpace, keyinfo
from ..view import makeclip, makecylinderslice


def plotOr(
    input,
    r: list[float],
    theta: float,
    fieldunits: dict,
    ignored_keys: list[str],
    basedir: str,
    printed: bool = True,
    marker: str = None,
    axs: dict = None,  # dict of fig ax for each field
    greyspace: bool = False,
    argsfield: str = None,
) -> dict:
    """plot vs r for a given theta

    Args:
        input: paraview reader
        r (list[float]): [r_start, r_end]
        theta (float): angle of the plot in degrees
        fieldunits (dict): dict of field units
        ignored_keys (list[str]): list of ignored fields
        basedir (str): result directory
        printed (bool, optional): Defaults to True.
        marker (str, optional): plot on specific marker. Defaults to None.
        axs (dict, optional): dict containing fig,ax,legend for each exported fields. Defaults to None.
        greyspace (bool, optional): plot grey vertical bars to fill holes in plots (slits/channels). Defaults to False.
        argsfield (str, optional): selected field to display. Defaults to None.

    Returns:
        dict: contains fig,ax,legend for each exported fields
    """

    [r0, r1] = r
    radian = theta * pi / 180.0

    cellDatatoPointData1 = CellDatatoPointData(
        registrationName="CellDatatoPointData", Input=input
    )

    plotOverLine = PlotOverLine(registrationName="Oz", Input=cellDatatoPointData1)

    # init the 'Line' selected for 'Source'
    plotOverLine.Point1 = [r0 * cos(radian), r0 * sin(radian), 0]
    plotOverLine.Point2 = [r1 * cos(radian), r1 * sin(radian), 0]

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
        axs=None,
        marker: str = None,
        greyspace: bool = False,
    ):
        [fig, ax, legend] = axs
        print(f"plotOrField: file={file}, key={key}", flush=True)
        (toolbox, physic, fieldname) = keyinfo(key)
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
        keycsv["r"] = keycsv["r"] + r0
        # print(f"new keys: {keycsv.columns.values.tolist()}")

        # rescale columns to plot
        units = {fieldname: fieldunits[fieldname]["Units"]}
        values = keycsv[key].to_list()
        out_values = convert_data(units, values, fieldname)
        ndf = {key: [val for val in out_values]}

        keycsv[key] = ndf[key]
        keycsv.plot(x="r", y=key, marker=marker, grid=True, ax=ax)
        legend.append(f"theta={theta:.2f}deg")
        if greyspace:
            legend = plot_greySpace(keycsv, "r", key, ax, legend)

        ax.set_xlabel("r [m]", fontsize=18)
        ax.set_ylabel(rf"{msymbol} [{out_unit:~P}]", fontsize=18)

        # ax.yaxis.set_major_locator(MaxNLocator(10))

        keycsv.to_csv(f"{basedir}/plots/{key}-vs-r-theta={theta}deg.csv")
        return legend

    # requirements: create PointData from CellData
    # for field in input.PointData.keys():
    datadict = resultinfo(cellDatatoPointData1, ignored_keys)
    for field in datadict["PointData"]["Arrays"]:
        if not field in ignored_keys and (not argsfield or field.startswith(argsfield)):
            kdata = datadict["PointData"]["Arrays"][field]
            Components = kdata["Components"]
            print(f"plotOrField for {field} - components={Components}", flush=True)
            if Components > 1:
                for i in range(Components):
                    print(f"plotOrField for {field}:{i} skipped", flush=True)
            else:
                if field not in axs:  # if field not in dict -> create fig, ax,legend
                    fig, ax = plt.subplots(figsize=(12, 8))
                    axs[field] = [fig, ax, []]
                axs[field][2] = plotOrField(
                    filename,
                    field,
                    theta,
                    fieldunits,
                    basedir,
                    axs=axs[field],
                    marker=marker,
                    greyspace=greyspace,
                )

    # remove temporary csv files
    os.remove(filename)

    Delete(plotOverLine)
    del plotOverLine

    Delete(cellDatatoPointData1)
    del cellDatatoPointData1
    return axs


def plotTheta(
    input,
    r: float,
    fieldunits: dict,
    ignored_keys: list[str],
    basedir: str,
    printed: bool = True,
    verbose: bool = False,
    marker: str = None,
    axs: dict = None,  # dict of fig,ax for each field
    argsfield: str = None,
) -> dict:
    """plot along theta for a given r

    for theta, need to apply CellDataToPointData filter

    Args:
        input: paraview reader
        r (float): r coordinates in m
        fieldunits (dict): dict of field units
        ignored_keys (list[str]): list of ignored fields
        basedir (str): result directory
        printed (bool, optional): Defaults to True.
        verbose (bool, optional): _description_. Defaults to False.
        marker (str, optional): plot on specific marker. Defaults to None.
        axs (dict, optional): dict containing fig,ax,legend for each exported fields. Defaults to None.
        argsfield (str, optional): selected field to display. Defaults to None.

    Returns:
        dict: contains fig,ax,legend for each exported fields
    """

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
        if clip.PointData.keys():
            slice = makecylinderslice(clip, f"slice{i}", r)
            SetActiveSource(slice)

            export = CreateWriter(
                f"{basedir}/plots/r={r}m-{i}.csv",
                proxy=slice,
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
        axs=None,
        marker: str = None,
    ):
        [fig, ax, legend] = axs
        print(f"plotThetaField: files={files}, key={key}", flush=True)
        (toolbox, physic, fieldname) = keyinfo(key)
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
            keycsv = csv[["Points:0", "Points:1", key]]

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
        df.plot(x="theta", y=key, marker=marker, grid=True, ax=ax)
        legend.append(f"r={r_mm:.0f}{mm}")

        ax.set_xlabel(r"$\theta$ [deg]", fontsize=18)
        ax.set_ylabel(rf"{msymbol} [{out_unit:~P}]", fontsize=18)

        # x range
        # x0 = -180
        # x1 = 180
        # ax.set_xlim(x0, x1)

        ax.yaxis.set_major_locator(MaxNLocator(10))

        print(f"{key} stats on r={r_mm}{mm}", flush=True)
        print(f"{df[key].describe()}", flush=True)
        return legend

    # requirements: create PointData from CellData
    datadict = resultinfo(cellDatatoPointData1, ignored_keys)
    for field in datadict["PointData"]["Arrays"]:
        if not field in ignored_keys and (not argsfield or field.startswith(argsfield)):
            kdata = datadict["PointData"]["Arrays"][field]
            Components = kdata["Components"]
            print(f"plotThetaField for {field} - components={Components}", flush=True)
            if Components > 1:
                for i in range(Components):
                    print(f"plotThetaField for {field}:{i} skipped", flush=True)
                # for i in range(Components):
                #     plotThetaField(files, f"{field}:{i}", r, z, ax=None, show=show)
            else:
                if field not in axs:  # if field not in dict -> create fig, ax
                    fig, ax = plt.subplots(figsize=(12, 8))
                    axs[field] = [fig, ax, []]
                axs[field][2] = plotThetaField(
                    files,
                    field,
                    r,
                    fieldunits,
                    basedir,
                    axs=axs[field],
                    marker=marker,
                )

    # remove temporary csv files
    for file in files:
        os.remove(file)

    Delete(cellDatatoPointData1)
    del cellDatatoPointData1

    # Force a garbage collection
    collected = gc.collect()
    if verbose:
        print(f"Garbage collector: collected {collected} objects.", flush=True)

    return axs


def makeplot(args, cellsize, fieldunits: dict, ignored_keys: list[str], basedir: str):
    """different plot situations for 2D

    * if 2 args.r and args.theta: plot Or from r0 to r1 at theta=args.theta
    * if not args.theta and args.r: plot Otheta ar r=args.r

    Args:
        args: options
        cellsize: paraview reader
        fieldunits (dict): dict of field units
        ignored_keys (list[str]): list of ignored fields
        basedir (str): result directory
    """
    if args.r:
        title = ""
        if fieldunits["Current"]["Val"]:
            title = title + f"\nI={fieldunits['Current']['Val']}"
        if fieldunits["B0"]["Val"]:
            title = title + f"\nB0={fieldunits['B0']['Val']}T"
        if fieldunits["Bbg"]["Val"]:
            title = title + f"\nBackground field: {fieldunits['Bbg']['Val']}"
        if len(args.r) != 2 or not args.theta:
            figaxs = {}
            for r in args.r:
                figaxs = plotTheta(
                    cellsize,
                    r,
                    fieldunits,
                    ignored_keys,
                    basedir,
                    marker=args.plotsMarker,
                    axs=figaxs,
                    argsfield=args.field,
                )

            showplot(figaxs, f"-vs-theta", basedir, title=title, show=args.show)
            plt.close()

        if args.theta and len(args.r) == 2:
            figaxs = {}
            for theta in args.theta:
                figaxs = plotOr(
                    cellsize,
                    args.r,
                    theta,
                    fieldunits,
                    ignored_keys,
                    basedir,
                    marker=args.plotsMarker,
                    axs=figaxs,
                    greyspace=args.greyspace,
                    argsfield=args.field,
                )  # with r=[r1, r2]

            showplot(
                figaxs,
                f"-vs-r",
                basedir,
                title=title,
                show=args.show,
            )
            plt.close()
