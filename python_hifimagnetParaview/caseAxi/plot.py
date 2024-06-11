import pandas as pd
import os
import matplotlib.pyplot as plt

from typing import List

from paraview.simple import CellDatatoPointData, PlotOverLine, CreateWriter, Delete

from ..method import convert_data, resultinfo, plot_greySpace


def plotOr(
    input,
    r: list[float],
    theta: float,
    z: float,
    fieldunits: dict,
    ignored_keys: List[str],
    basedir: str,
    printed: bool = True,
    marker: str = None,
    axs: dict = None,  # dict of fig,ax for each field
    greyspace: bool = False,
) -> dict:

    [r0, r1] = r

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
    plotOverLine.Point1 = [r0, z, 0]
    plotOverLine.Point2 = [r1, z, 0]

    filename = f"{basedir}/plots/r0={r0}m-r1={r1}m-z={z}m.csv"
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
        z: float,
        fieldunits: dict,
        axs=None,
        marker: str = None,
        greyspace: bool = False,
    ):
        [fig, ax, legend] = axs
        print(f"plotOrField: file={file}, key={key}", flush=True)
        keyinfo = key.replace("_Magnitude", "").split(".")
        # print(f"keyinfo={keyinfo}", flush=True)
        if len(keyinfo) == 1:
            fieldname = field
        elif len(keyinfo) == 2:
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

        z_units = {"coord": fieldunits["coord"]["Units"]}
        mm = f'{fieldunits["coord"]["Units"][1]:~P}'
        z_mm = convert_data(z_units, z, "coord")

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
        legend.append(f"z={z_mm:.0f}{mm}")
        if greyspace:
            legend = plot_greySpace(keycsv, "r", key, ax, legend)

        ax.set_xlabel("r [m]", fontsize=18)
        ax.set_ylabel(rf"{symbol} [{out_unit:~P}]", fontsize=18)

        # ax.yaxis.set_major_locator(MaxNLocator(10))
        return legend

    # requirements: create PointData from CellData
    # for field in input.PointData.keys():
    datadict = resultinfo(cellDatatoPointData1, ignored_keys)
    for field in datadict["PointData"]["Arrays"]:
        if not field in ignored_keys:
            kdata = datadict["PointData"]["Arrays"][field]
            Components = kdata["Components"]
            print(f"plotOrField for {field} - components={Components}", flush=True)
            if Components > 1:
                for i in range(Components):
                    print(f"plotOrField for {field}:{i} skipped", flush=True)
            else:
                if field not in axs:  # if field not in dict -> create fig, ax
                    fig, ax = plt.subplots(figsize=(12, 8))
                    axs[field] = [fig, ax, []]
                axs[field][2] = plotOrField(
                    filename,
                    field,
                    z,
                    fieldunits,
                    axs=axs[field],
                    marker=marker,
                    greyspace=greyspace,
                )

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
    ignored_keys: List[str],
    basedir: str,
    printed: bool = True,
    marker: str = None,
    axs: dict = None,  # dict of fig,ax for each field
) -> dict:

    [z0, z1] = z

    cellDatatoPointData1 = CellDatatoPointData(
        registrationName="CellDatatoPointData", Input=input
    )

    plotOverLine = PlotOverLine(registrationName="Oz", Input=cellDatatoPointData1)

    # init the 'Line' selected for 'Source'
    plotOverLine.Point1 = [r, z0, 0]
    plotOverLine.Point2 = [r, z1, 0]

    filename = f"{basedir}/plots/r={r}m-z0={z0}m-z1={z1}m.csv"
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
        r: float,
        fieldunits: dict,
        axs=None,
        marker: str = None,
    ):
        [fig, ax, legend] = axs
        print(f"plotOrField: file={file}, key={key}", flush=True)
        keyinfo = key.replace("_Magnitude", "").split(".")
        # print(f"keyinfo={keyinfo}", flush=True)
        if len(keyinfo) == 1:
            fieldname = field
        elif len(keyinfo) == 2:
            (physic, fieldname) = keyinfo
        elif len(keyinfo) == 3:
            (toolbox, physic, fieldname) = keyinfo
        else:
            raise RuntimeError(f"{key}: cannot get keyinfo as splitted char")
        symbol = fieldunits[fieldname]["Symbol"]
        [in_unit, out_unit] = fieldunits[fieldname]["Units"]

        if ax is None:
            ax = plt.gca()

        r_units = {"coord": fieldunits["coord"]["Units"]}
        mm = f'{fieldunits["coord"]["Units"][1]:~P}'
        r_mm = convert_data(r_units, r, "coord")

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
        legend.append(f"r={r_mm:.0f}{mm}")

        ax.set_xlabel("z [m]", fontsize=18)
        ax.set_ylabel(rf"{symbol} [{out_unit:~P}]", fontsize=18)

        # ax.yaxis.set_major_locator(MaxNLocator(10))
        return legend

    # requirements: create PointData from CellData
    # for field in input.PointData.keys():
    datadict = resultinfo(cellDatatoPointData1, ignored_keys)
    for field in datadict["PointData"]["Arrays"]:
        if not field in ignored_keys:
            kdata = datadict["PointData"]["Arrays"][field]
            Components = kdata["Components"]
            print(f"plotOzField for {field} - components={Components}", flush=True)
            if Components > 1:
                for i in range(Components):
                    print(f"plotOzField for {field}:{i} skipped", flush=True)
            else:
                if field not in axs:  # if field not in dict -> create fig, ax,legend
                    fig, ax = plt.subplots(figsize=(12, 8))
                    axs[field] = [fig, ax, []]
                axs[field][2] = plotOzField(
                    filename,
                    field,
                    r,
                    fieldunits,
                    axs=axs[field],
                    marker=marker,
                )
    return axs
