import gc
import os
import re
import pandas as pd

from paraview.simple import (
    OpenDataFile,
    UpdatePipeline,
    Calculator,
    IntegrateVariables,
    CreateView,
    Show,
    ExportView,
    Delete,
    ProbeLocation,
    SaveData,
)

from pint import Quantity

# Ignore warning for pint
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Quantity([])


def convert_data(
    units: dict, quantity: float | list[float], qtype: str, debug: bool = False
) -> float | list[float]:
    """Returns quantity unit consistant with length unit

    Args:
        units (dict): dict of the form {field: units}
        quantity (float | list[float]): quantity to convert
        qtype (str): name of quantity/field
        debug (bool, optional): print debug. Defaults to False.

    Raises:
        Exception: convert_data/quantity: unsupported type

    Returns:
        float | list[float]: converted quantity
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
) -> float | list[float]:
    """Returns quantity unit to its original unit consistant with length unit

    Args:
        units (dict): dict of the form {field: units}
        quantity (float | list[float]): quantity to convert
        qtype (str): name of quantity/field
        debug (bool, optional): print debug. Defaults to False.

     Raises:
        Exception: convert_data/quantity: unsupported type

    Returns:
        float | list[float]: converted quantity
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


def selectBlocks(blockdata: list, excludes: list[str]) -> list[str]:
    """select Block name not in exludes

    Args:
        blockdata (list): list of block names
        excludes (list[str]): list of block names to be excluded

    Returns:
        list[str]: list of included block names
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
    """returns info about input dataset

    Args:
        input: paraview reader

    Returns:
        dataInfo
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


def resultinfo(input, ignored_keys: list[str], verbose: bool = False) -> dict:
    """returns a dict gathering info on PointData, CellData and FieldData

    Args:
        input: paraview reader
        ignored_keys (list[str]): list of ignored key
        verbose (bool, optional): print verbose. Defaults to False.

    Returns:
        dict: info dictionnary
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
    key, ignored_keys: list[str], verbose: bool = False, printed: bool = True
) -> dict:
    """print info on key

    Args:
        key: key object
        ignored_keys (list[str]): list of ignored keys
        verbose (bool, optional): print verbose. Defaults to False.
        printed (bool, optional): Defaults to True.

    Returns:
        dict: info dictionnary
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


def keyinfo(key: str) -> tuple:
    keyinfo = key.split(".")
    if len(keyinfo) == 1:
        fieldname = key
        toolbox = None
        physic = None
    elif len(keyinfo) == 2:
        toolbox = None
        (physic, fieldname) = keyinfo
    elif len(keyinfo) == 3:
        (toolbox, physic, fieldname) = keyinfo
    else:
        raise RuntimeError(f"{key}: cannot get keyinfo as splitted char")

    return (toolbox, physic, fieldname)


def getbounds(input):
    """returns bounds of input geometry

    Args:
        input: paraview reader

    Returns:
        bounds
    """

    dataInfo = input.GetDataInformation()
    bounds = dataInfo.DataInformation.GetBounds()
    return bounds


def load(file: str, printed: bool = True):
    """create dataset from file

    Args:
        file (str): file name
        printed (bool, optional): Defaults to True.

    Returns:
        paraview reader
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
    """compute moment of order N

    Args:
        input: paraview reader
        key (str): field name
        nkey (str): field surname
        order (int): order of moment
        AttributeType (str): Point Data or Cell Data

    Returns:
        calculator
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
) -> str:
    """compute integral of Data over area/volume

    to get values use a spreadsheet

    Args:
        input: paraview reader
        name (str): _description_
        basedir (str): result directory
        printed (bool, optional): Defaults to True.
        verbose (bool, optional): print verbose. Defaults to False.

    Returns:
        str: file name
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


def showplot(
    figaxs: dict, suffix: str, basedir: str, title: str = "", show: bool = True
):
    """show or save plot for each field in figaxs dictionnary

    Args:
        figaxs (dict): dict containing fig,ax,legend for each exported fields
        suffix (str): file and title suffix
        basedir (str): result directory
        title (str): title suffix
        show (bool, optional): show plot (False=Save plot). Defaults to True.
    """

    for f in figaxs:
        axs = figaxs[f]
        axs[1].legend(axs[2], fontsize=18, loc="best")
        axs[1].set_title(f"{f}{suffix}{title}", fontsize=20)
        axs[1].grid(True, linestyle="--")
        axs[1].tick_params(axis="x", which="major", labelsize=15)
        axs[1].tick_params(axis="y", which="major", labelsize=15)
        if show:
            print(f"show {f}{suffix}", flush=True)
            axs[0].show()
        else:
            print(f"save {f}{suffix}.png", flush=True)
            axs[0].tight_layout()
            axs[0].savefig(f"{basedir}/plots/{f}{suffix}.png", dpi=300)


def plot_greySpace(df: pd.DataFrame, cx: str, cy: str, ax, legend: list[str]):
    """plot grey vertical bars to fill holes in plots (slits/channels) (works in 2D vs r)


    Args:
        df (pd.DataFrame): DataFrame of the plot
        cx (str): column name of xaxis in plot
        cy (str): column name of yaxis in plot
        ax: ax of the plot
        legend (list[str]):  legend (to hide grey bars)

    Returns:
        list[str]: updated legend
    """

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


def getcurrent(current: str, marker: str = None):
    if current.endswith(".csv"):
        df = pd.read_csv(current)
        for col in df.columns():
            if "intensity" in col or "Intensity" in col:
                current = f"{int(df[col].iloc[-1])}A"
                break
    else:
        current += "A"

    return current


def getB0(reader, fieldtype: dict, basedir: str, dim: int, axis: bool = False) -> float:
    B0 = None
    for f in fieldtype:
        print(f)
        print(fieldtype[f]["Type"])
        if fieldtype[f]["Type"] == "MagneticField":
            B0 = f
            break

    if not B0:
        return None
    # create a new 'Probe Location'
    probeLocation = ProbeLocation(
        registrationName="ProbeLocation",
        Input=reader,
        ProbeType="Fixed Radius Point Source",
    )

    # init the 'Fixed Radius Point Source' selected for 'ProbeType'
    probeLocation.ProbeType.Center = [0.0, 0.0, 0.0]

    for key in list(probeLocation.PointData.keys()) + list(
        probeLocation.CellData.keys()
    ):
        (toolbox, physic, fieldname) = keyinfo(key)
        if B0 == fieldname and key in list(probeLocation.CellData.keys()):
            SaveData(
                f"{basedir}insert-B0.csv",
                proxy=probeLocation,
                ChooseArraysToWrite=1,
                CellDataArrays=[key],
            )
            savedkey = key
        if B0 == fieldname and key in list(probeLocation.PointData.keys()):
            SaveData(
                f"{basedir}/insert-B0.csv",
                proxy=probeLocation,
                ChooseArraysToWrite=1,
                PointDataArrays=[key],
            )
            savedkey = key

    try:
        df = pd.read_csv(f"{basedir}/insert-B0.csv")

        if axis:
            B0 = abs(df[f"{savedkey}:1"].iloc[-1])
        elif dim == 2:
            B0 = abs(df[f"{savedkey}:1"].iloc[-1])
        elif dim == 3:
            B0 = abs(df[f"{savedkey}:2"].iloc[-1])
        return round(B0, 1)
    except:
        return None
