import pandas as pd
import os
import gc
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from numpy import pi, arctan2
from math import pi, cos, sin

from paraview.simple import (
    PlotOnIntersectionCurves,
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
    z: float,
    fieldunits: dict,
    ignored_keys: list[str],
    basedir: str,
    printed: bool = True,
    marker: str = None,
    axs: dict = None,  # dict of fig,ax for each field
    greyspace: bool = False,
    argsfield: str = None,
) -> dict:
    """plot vs r for a given theta and a given z

    Args:
        input: paraview reader
        r (list[float]): [r_start, r_end]
        theta (float): angle of the plot in degree
        z (float): z coordinate of the plot in m
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
    # get params list
    if not printed:
        for prop in plotOverLine.ListProperties():
            print(
                f"plotOverLine': {prop}={plotOverLine.GetPropertyValue(prop)}",
                flush=True,
            )

    # init the 'Line' selected for 'Source'
    plotOverLine.Point1 = [r0 * cos(radian), r0 * sin(radian), z]
    plotOverLine.Point2 = [r1 * cos(radian), r1 * sin(radian), z]

    filename = f"{basedir}/plots/r0={r0}m-r1={r1}m-theta={theta}deg-z={z}m.csv"
    export = CreateWriter(filename, proxy=plotOverLine)
    if not printed:
        for prop in export.ListProperties():
            print(f"export: {prop}={export.GetPropertyValue(prop)}", flush=True)
    export.UpdateVTKObjects()  # is it needed?
    export.UpdatePipeline()

    # plot with matplotlib
    def plotOrField(
        file,
        key: str,
        theta: float,
        z: float,
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
        # print(f"new keys: {keycsv.columns.values.tolist()}", flush=True)

        # rescale columns to plot
        units = {fieldname: fieldunits[fieldname]["Units"]}
        values = keycsv[key].to_list()
        out_values = convert_data(units, values, fieldname)
        ndf = {key: [val for val in out_values]}

        keycsv[key] = ndf[key]
        keycsv.plot(x="r", y=key, marker=marker, grid=True, ax=ax)

        ax.set_xlabel("r [m]", fontsize=18)
        ax.set_ylabel(rf"{msymbol} [{out_unit:~P}]", fontsize=18)

        # ax.yaxis.set_major_locator(MaxNLocator(10))
        r_units = {"coord": fieldunits["coord"]["Units"]}
        mm = f'{fieldunits["coord"]["Units"][1]:~P}'
        z_mm = convert_data(r_units, z, "coord")
        legend.append(f"z={z_mm:.0f}{mm}")
        if greyspace:
            legend = plot_greySpace(keycsv, "r", key, ax, legend)

        keycsv.to_csv(f"{basedir}/plots/{key}-vs-r-theta={theta}deg-z={z_mm}{mm}.csv")
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
                if field not in axs:  # if field not in dict -> create fig, ax, legend
                    fig, ax = plt.subplots(figsize=(12, 8))
                    axs[field] = [fig, ax, []]
                axs[field][2] = plotOrField(
                    filename,
                    field,
                    theta,
                    z,
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


def plotOz(
    input,
    r: float,
    theta: float,
    z: list[float],
    fieldunits: dict,
    ignored_keys: list[str],
    basedir: str,
    printed: bool = True,
    marker: str = None,
    axs: dict = None,  # dict of fig,ax for each field
    argsfield: str = None,
) -> dict:
    """plot along z for a given r and for a given theta

    Args:
        input: paraview reader
        r (float): r coordinates in m
        theta (float): angle of the plot in degrees
        z (list[float]): [z_start, z_end]
        fieldunits (dict): dict of field units
        ignored_keys (list[str]): list of ignored fields
        basedir (str): result directory
        printed (bool, optional): Defaults to True.
        marker (str, optional): plot on specific marker. Defaults to None.
        axs (dict, optional): dict containing fig,ax,legend for each exported fields. Defaults to None.
        argsfield (str, optional): selected field to display. Defaults to None.

    Returns:
        dict: contains fig,ax,legend for each exported fields
    """
    [z0, z1] = z
    radian = theta * pi / 180.0

    cellDatatoPointData1 = CellDatatoPointData(
        registrationName="CellDatatoPointData", Input=input
    )

    plotOverLine = PlotOverLine(registrationName="Oz", Input=cellDatatoPointData1)

    # init the 'Line' selected for 'Source'
    plotOverLine.Point1 = [r * cos(radian), r * sin(radian), z0]
    plotOverLine.Point2 = [r * cos(radian), r * sin(radian), z1]

    filename = f"{basedir}/plots/r={r}m-theta={theta}deg-z0={z0}m-z1={z1}m.csv"
    export = CreateWriter(filename, proxy=plotOverLine)
    if not printed:
        for prop in export.ListProperties():
            print(f"export: {prop}={export.GetPropertyValue(prop)}", flush=True)
    export.UpdateVTKObjects()  # is it needed?
    export.UpdatePipeline()

    # plot with matplotlib
    def plotOzField(
        file,
        key: str,
        theta: float,
        z: list[float],
        fieldunits: dict,
        basedir: str,
        axs=None,
        marker: str = None,
    ):
        [fig, ax, legend] = axs
        print(f"plotOrField: file={file}, key={key}", flush=True)
        (toolbox, physic, fieldname) = keyinfo(key)
        symbol = fieldunits[fieldname]["Symbol"]
        msymbol = symbol
        if "mSymbol" in fieldunits[fieldname]:
            msymbol = fieldunits[fieldname]["mSymbol"]
        [in_unit, out_unit] = fieldunits[fieldname]["Units"]

        if ax is None:
            ax = plt.gca()

        # see vonmises-vs-theta.py and/or vonmises-vs-theta-plot-savedata.py
        csv = pd.read_csv(file)
        keycsv = csv[["arc_length", key]]

        # rename columns
        keycsv.rename(columns={"arc_length": "z"}, inplace=True)
        keycsv["z"] = keycsv["z"] + z0
        # print(f"new keys: {keycsv.columns.values.tolist()}", flush=True)

        # rescale columns to plot
        units = {fieldname: fieldunits[fieldname]["Units"]}
        values = keycsv[key].to_list()
        out_values = convert_data(units, values, fieldname)
        ndf = {key: [val for val in out_values]}

        keycsv[key] = ndf[key]
        keycsv.plot(x="z", y=key, marker=marker, grid=True, ax=ax)
        legend.append(f"theta={theta}deg")

        ax.set_xlabel("z [m]", fontsize=18)
        ax.set_ylabel(rf"{msymbol} [{out_unit:~P}]", fontsize=18)

        # ax.yaxis.set_major_locator(MaxNLocator(10))
        r_units = {"coord": fieldunits["coord"]["Units"]}
        mm = f'{fieldunits["coord"]["Units"][1]:~P}'
        r_mm = convert_data(r_units, r, "coord")

        keycsv.to_csv(f"{basedir}/plots/{key}-vs-z-theta={theta}deg-r={r_mm}{mm}.csv")
        return legend

    # requirements: create PointData from CellData
    # for field in input.PointData.keys():
    datadict = resultinfo(cellDatatoPointData1, ignored_keys)
    for field in datadict["PointData"]["Arrays"]:
        if not field in ignored_keys and (not argsfield or field.startswith(argsfield)):
            kdata = datadict["PointData"]["Arrays"][field]
            Components = kdata["Components"]
            print(f"plotOzField for {field} - components={Components}", flush=True)
            if Components > 1:
                for i in range(Components):
                    print(f"plotOzField for {field}:{i} skipped", flush=True)
            else:
                if field not in axs:  # if field not in dict -> create fig, ax
                    fig, ax = plt.subplots(figsize=(12, 8))
                    axs[field] = [fig, ax, []]
                axs[field][2] = plotOzField(
                    filename,
                    field,
                    theta,
                    z,
                    fieldunits,
                    basedir,
                    axs=axs[field],
                    marker=marker,
                )

    # remove temporary csv files
    os.remove(filename)

    Delete(plotOverLine)
    del plotOverLine

    Delete(cellDatatoPointData1)
    del cellDatatoPointData1
    return axs


# @profile
def plotTheta(
    input,
    r: float,
    z: float,
    fieldunits: dict,
    ignored_keys: list[str],
    basedir: str,
    printed: bool = True,
    verbose: bool = False,
    marker: str = None,
    axs: dict = None,  # dict of fig,ax for each field
    argsfield: str = None,
) -> dict:
    """plot along theta for a given r and for a given z

    Args:
        input: paraview reader
        r (float): r coordinates in m
        z (float): z coordinates in m
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

    print(f"plotTheta: r={r}, z={z}", flush=True)
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
        plotOnIntersectionCurve.SliceType.Origin = [0.0, 0.0, z]
        plotOnIntersectionCurve.SliceType.Normal = [0.0, 0.0, 1.0]

        # plotOnIntersectionCurves.append(plotOnIntersectionCurve)
        export = CreateWriter(
            f"{basedir}/plots/r={r}m-z={z}m-{i}.csv", proxy=plotOnIntersectionCurve
        )
        if not printed:
            for prop in export.ListProperties():
                print(f"export: {prop}={export.GetPropertyValue(prop)}", flush=True)
        export.UpdateVTKObjects()  # is it needed?
        export.UpdatePipeline()
        files.append(f"{basedir}/plots/r={r}m-z={z}m-{i}.csv")

    # plot with matplotlib
    def plotThetaField(
        files,
        key: str,
        r: float,
        z: float,
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
        z_mm = convert_data(r_units, z, "coord")
        df.to_csv(f"{basedir}/plots/{key}-vs-theta-r={r_mm}{mm}-z={z_mm}{mm}.csv")
        # print(f"df keys: {df.columns.values.tolist()}", flush=True)
        # assert key in df.columns.values.tolist(), f"{key} not in df_keys"
        # print(f"df={df}", flush=True)
        df.plot(x="theta", y=key, marker=marker, grid=True, ax=ax)
        legend.append(f"z={z_mm:.0f}{mm}")

        ax.set_xlabel(r"$\theta$ [deg]", fontsize=18)
        ax.set_ylabel(rf"{msymbol} [{out_unit:~P}]", fontsize=18)

        # x range
        x0 = -180
        x1 = 180
        ax.set_xlim(x0, x1)

        ax.yaxis.set_major_locator(MaxNLocator(10))

        print(f"{key} stats on r={r_mm}{mm}-z={z_mm}{mm}", flush=True)
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
                #     plotThetaField(files, f"{field}:{i}", r, z, ax=None)
            else:
                if field not in axs:  # if field not in dict -> create fig, ax
                    fig, ax = plt.subplots(figsize=(12, 8))
                    axs[field] = [fig, ax, []]
                axs[field][2] = plotThetaField(
                    files,
                    field,
                    r,
                    z,
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
    """different plot situations for 3D

    * if args.r and args.z: plot Otheta
    * if 2 args.r and args.theta and args.z: plot Or from r0 to r1 at theta=args.theta, z=args.z
    * if 2 args.z and args.theta and args.r: plot Oz from z0 to z1 at theta=args.theta, r=args.r

    Args:
        args: options
        cellsize: paraview reader
        fieldunits (dict): dict of field units
        ignored_keys (list[str]): list of ignored fields
        basedir (str): result directory
    """
    if args.r and args.z:
        title = ""
        if fieldunits["Current"]["Val"]:
            title = title + f"\nI={fieldunits['Current']['Val']}"
        if fieldunits["B0"]["Val"]:
            title = title + f"\nB0={fieldunits['B0']['Val']}T"
        if fieldunits["Bbg"]["Val"]:
            title = title + f"\nBackground field: {fieldunits['Bbg']['Val']}"
        for r in args.r:
            figaxs = {}
            for z in args.z:
                figaxs = plotTheta(
                    cellsize,
                    r,
                    z,
                    fieldunits,
                    ignored_keys,
                    basedir,
                    marker=args.plotsMarker,
                    axs=figaxs,
                    argsfield=args.field,
                )

            showplot(
                figaxs,
                f"-vs-theta-r={r}m",
                basedir,
                title=title,
                show=args.show,
            )
            plt.close()

            if args.theta and len(args.z) == 2:
                figaxs = {}
                for theta in args.theta:
                    figaxs = plotOz(
                        cellsize,
                        r,
                        theta,
                        args.z,
                        fieldunits,
                        ignored_keys,
                        basedir,
                        marker=args.plotsMarker,
                        axs=figaxs,
                        argsfield=args.field,
                    )  # with r: float, z=[z1,z2]

                showplot(
                    figaxs,
                    f"-vs-z-r={r}m",
                    basedir,
                    title=title,
                    show=args.show,
                )
                plt.close()
        if args.theta and len(args.r) == 2:
            for theta in args.theta:
                figaxs = {}
                for z in args.z:
                    figaxs = plotOr(
                        cellsize,
                        args.r,
                        theta,
                        z,
                        fieldunits,
                        ignored_keys,
                        basedir,
                        marker=args.plotsMarker,
                        axs=figaxs,
                        greyspace=args.greyspace,
                        argsfield=args.field,
                    )  # with r=[r1, r2], z: float

                showplot(
                    figaxs,
                    f"-vs-r-theta={theta}deg",
                    basedir,
                    title=title,
                    show=args.show,
                )
                plt.close()
