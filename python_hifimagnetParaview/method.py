import gc
import os
import re
import pandas as pd
from typing import List

from paraview.simple import (
    OpenDataFile,
    UpdatePipeline,
    Calculator,
    IntegrateVariables,
    CreateView,
    Show,
    ExportView,
    Delete,
)

from pint import Quantity

# Ignore warning for pint
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Quantity([])


def convert_data(
    units: dict, quantity: float | list[float], qtype: str, debug: bool = False
):
    """
    Returns quantity unit consistant with length unit
    """

    data = None
    if isinstance(quantity, float):
        data = Quantity(quantity, units[qtype][0]).to(units[qtype][1]).magnitude
        if debug:
            print(qtype, quantity, "data=", data, flush=True)
    elif isinstance(quantity, list):
        data = (
            Quantity(quantity, units[qtype][0]).to(units[qtype][1]).magnitude.tolist()
        )
    else:
        raise Exception(
            f"convert_data/quantity: unsupported type {type(quantity)} for {qtype}"
        )

    return data


def invert_convert_data(
    units: dict, quantity: float | list[float], qtype: str, debug: bool = False
):
    """
    Returns quantity unit consistant with length unit
    """

    data = None
    if isinstance(quantity, float):
        data = Quantity(quantity, units[qtype][1]).to(units[qtype][0]).magnitude
        if debug:
            print(qtype, quantity, "data=", data, flush=True)
    elif isinstance(quantity, list):
        data = (
            Quantity(quantity, units[qtype][1]).to(units[qtype][0]).magnitude.tolist()
        )
    else:
        raise Exception(
            f"convert_data/quantity: unsupported type {type(quantity)} for {qtype}"
        )

    return data


def selectBlocks(blockdata: list, excludes: list):
    """
    select Block name not in exludes

    blockdata: list of block names
    excludes: list of names to be excluded

    """
    isolants = []
    if "Isolant" in excludes:
        isolants = [
            r"/Root/\w+_B0",
            r"/Root/\w+_B4",
            r"/Root/\w*H\d+_Cu0",
            r"/Root/\w*H\d+_Cu21",
            r"/Root/\w*R\d+",
        ]
    selectedblocks = []
    for key in blockdata:
        found = False
        for excluded in excludes:
            if excluded in key:
                found = True
                break
        if isolants:
            for excluded in isolants:
                regex = re.compile(excluded)
                if regex.match(key):
                    found = True
                    break
        if not found:
            selectedblocks.append(key)

    return selectedblocks


def info(input):
    """
    returns info about input dataset
    """
    print("info:", flush=True)

    dataInfo = input.GetDataInformation()
    print(f"Memory: {dataInfo.DataInformation.GetMemorySize()}", flush=True)
    print(f"Bounds: {dataInfo.DataInformation.GetBounds()}", flush=True)
    print(f"Nodes: {dataInfo.GetNumberOfPoints()}", flush=True)
    print(f"Cells: {dataInfo.GetNumberOfCells()}", flush=True)
    print(f"pointData: {input.PointData[:]}", flush=True)
    print(f"cellData: {input.CellData[:]}", flush=True)
    print(f"fieldData: {input.FieldData[:]}", flush=True)
    return dataInfo


def resultinfo(input, ignored_keys: List[str], verbose: bool = False) -> dict:
    """
    returns a dict gathering info on PointData, CellData and FieldData
    """
    if verbose:
        print(f"resultinfo: input={input}", flush=True)

    datadict = {
        "PointData": {"TypeMode": "POINT", "AttributeMode": "Point Data", "Arrays": {}},
        "CellData": {"TypeMode": "CELL", "AttributeMode": "Cell Data", "Arrays": {}},
        "FieldData": {"TypeMode": "FIELD", "AttributeMode": "Field Data", "Arrays": {}},
    }

    if verbose:
        print("PointData:", flush=True)
    for key in input.PointData:
        datadict["PointData"]["Arrays"][key.Name] = getresultInfo(
            key, ignored_keys, verbose
        )
    if verbose:
        print("CellData:", flush=True)
    for key in input.CellData:
        datadict["CellData"]["Arrays"][key.Name] = getresultInfo(
            key, ignored_keys, verbose
        )
    if verbose:
        print("FieldData:", flush=True)
    for key in input.FieldData:
        datadict["FieldData"]["Arrays"][key.Name] = getresultInfo(
            key, ignored_keys, verbose
        )

    return datadict


def getresultInfo(
    key, ignored_keys: List[str], verbose: bool = False, printed: bool = True
) -> dict:
    """
    print info on key

    """

    name = key.Name
    components = key.GetNumberOfComponents()
    bounds = []
    if key.GetNumberOfComponents() > 1:
        # print(f"{key}: range={key.GetRange(-1)}", flush=True)
        bounds.append(key.GetRange(-1))
        for i in range(key.GetNumberOfComponents()):
            bounds.append(key.GetRange(i))
            # print(f"{key}: component={i}, range={bounds[-1]}", flush=True)
    else:
        bounds.append(key.GetRange())

    datadict = {
        "Components": components,
        "Bounds": bounds,
    }
    if name not in ignored_keys and verbose:
        print(f"getresultInfo {name}: datadict={datadict}", flush=True)

    return datadict


def getbounds(input):
    """
    returns bounds of input geometry
    """

    dataInfo = input.GetDataInformation()
    bounds = dataInfo.DataInformation.GetBounds()
    return bounds


def load(file: str, printed: bool = True):
    """
    create dataset from file
    """

    print(f"Load Ensight case: {file}", flush=True)
    input = OpenDataFile(file)
    UpdatePipeline()

    if not printed:
        print(f"properties: {input.ListProperties()}", flush=True)
        for prop in input.ListProperties():
            print(f"reader: {prop}={input.GetPropertyValue(prop)}", flush=True)

    return input


def momentN(input, key: str, nkey: str, order: int, AttributeType: str):
    """
    compute moment of order N
    """

    calculator1 = Calculator(registrationName=f"moment{order}", Input=input)
    calculator1.AttributeType = AttributeType  # 'Cell Data'
    calculator1.ResultArrayName = f"{nkey}_moment{order}"
    # print(f"momentN{order}(key={key}, order={order}): {calculator1.ResultArrayName}")

    # check Points_0 for r
    if order == 1:
        calculator1.Function = f"coordsX*{key}"
    else:
        calculator1.Function = f"coordsX*{key}^{order}"

    calculator1.UpdatePipeline()
    return calculator1


def integrateKeys(
    input, name: str, basedir: str, printed: bool = True, verbose: bool = False
):
    """
    compute integral of Data over area/volume

    to get values use a spreadsheet
    """
    os.makedirs(f"{basedir}/stats", exist_ok=True)
    integratedvalues = IntegrateVariables(Input=input)
    if not printed:
        for prop in integratedvalues.ListProperties():
            print(
                f"integratedvalues: {prop}={integratedvalues.GetPropertyValue(prop)}",
                flush=True,
            )

    spreadSheetView = CreateView("SpreadSheetView")
    spreadSheetView.FieldAssociation = "Cell Data"

    descriptiveDisplay = Show(
        integratedvalues, spreadSheetView, "SpreadSheetRepresentation"
    )

    filename = f"{basedir}/stats/{name}-integrals.csv"
    spreadSheetView.Update()
    export = ExportView(
        filename,
        view=spreadSheetView,
        RealNumberNotation="Scientific",
    )

    Delete(spreadSheetView)
    del spreadSheetView
    Delete(descriptiveDisplay)

    # Force a garbage collection
    collected = gc.collect()
    if verbose:
        print(
            f"integrateKeys: Garbage collector: collected {collected} objects.",
            flush=True,
        )

    return filename


def showplot(figaxs: dict, title: str, basedir: str, show: bool = True):
    for f in figaxs:
        axs = figaxs[f]
        axs[1].legend(axs[2], fontsize=18, loc="best")
        axs[1].set_title(f"{f}{title} ", fontsize=20)
        axs[1].grid(True, linestyle="--")
        axs[1].tick_params(axis="x", which="major", labelsize=15)
        axs[1].tick_params(axis="y", which="major", labelsize=15)
        if show:
            print(f"show {f}{title}", flush=True)
            axs[0].show()
        else:
            print(f"save {f}{title}.png", flush=True)
            axs[0].tight_layout()
            axs[0].savefig(f"{basedir}/plots/{f}{title}.png", dpi=300)


def plot_greySpace(df: pd.DataFrame, cx: str, cy: str, ax, legend):
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
