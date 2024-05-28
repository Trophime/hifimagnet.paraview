import numpy as np
import vtk
import gc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import hist
from math import pi, sqrt

from tabulate import tabulate

from paraview.simple import *
from paraview.simple import (
    CellSize,
    CellCenters,
    DescriptiveStatistics,
    ExtractBlock,
    ExtractSurface,
    BoundingRuler,
    PointDatatoCellData,
    CellDatatoPointData,
    Calculator,
    WarpByVector,
    MergeBlocks,
    PlotOverLine,
    Text,
    IntegrateVariables,
)

from paraview import servermanager as sm
from paraview.vtk.numpy_interface import dataset_adapter as dsa
from paraview.vtk.numpy_interface import algorithms as algs

import argparse
import os

from pint import UnitRegistry, Unit, Quantity

# Ignore warning for pint
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Quantity([])

# Pint configuration
ureg = UnitRegistry()
ureg.define("percent = 0.01 = %")
ureg.define("ppm = 1e-6")
ureg.default_system = "SI"
ureg.autoconvert_offset_to_baseunit = True

# set default output unit to millimeter
distance_unit = "millimeter"  # or "meter"

# fieldunits: dict( Quantity: symbol: str, units: [ in_unit, out_unit ]
# TODO
# * use different symbol for Paraview, MatPlotlib not necessary the same
#   unless latex is used for matplotlib
#  !! watchout matplotlib setup for latex support and packages required !!!
# * add an exclude list: list of blockname to be excluded
# ex: for Temperature exclude=['Air'], J exclude=['Air', '*Isolant']
# use regexp to select block to be be excluded

# build fieldunits when creating setup and store as a json
# load json

# use r"$\theta$" for displaying units in mathplotlib
fieldunits = {
    "VolumicMass": {
        "Symbol": "rho",
        "mSymbol": r"$\rho$",
        "Units": [
            ureg.kilogram / ureg.meter**3,
            ureg.kilogram / ureg.Unit(distance_unit) ** 3,
        ],
        "Exclude": ["Air"],
    },
    "ThermalConductivity": {
        "Symbol": "k",
        "Units": [
            ureg.watt / ureg.meter / ureg.kelvin,
            ureg.watt / ureg.Unit(distance_unit) / ureg.kelvin,
        ],
        "Exclude": ["Air"],
    },
    "ElectricalConductivity": {
        "Symbol": "Sigma",
        "mSymobol": r"$\sigma$",
        "Units": [
            ureg.siemens / ureg.meter,
            ureg.siemens / ureg.Unit(distance_unit),
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "YoungModulus": {
        "Symbol": "E",
        "Units": [
            ureg.pascal,
            ureg.megapascal,
        ],
        "Exclude": ["Air"],
    },
    "Length": {"Symbol": "", "Units": [ureg.meter, ureg.Unit(distance_unit)]},
    "Area": {
        "Symbol": "S",
        "Units": [
            ureg.meter**2,
            ureg.Unit(distance_unit) ** 2,
        ],
    },
    "Volume": {
        "Symbol": "V",
        "Units": [
            ureg.meter**3,
            ureg.Unit(distance_unit) ** 3,
        ],
    },
    "mu0": {
        "Symbol": "mu",
        "Units": [ureg.henry / ureg.meter, ureg.henry / ureg.Unit(distance_unit)],
    },
    "h": {
        "Symbol": "hw",
        "Units": [
            ureg.watt / ureg.meter**2 / ureg.kelvin,
            ureg.watt / ureg.Unit(distance_unit) ** 2 / ureg.kelvin,
        ],
    },
    "Flow": {
        "Symbol": "Q",
        "Units": [
            ureg.liter / ureg.second,
            ureg.Unit(distance_unit)
            * ureg.Unit(distance_unit)
            * ureg.Unit(distance_unit)
            / ureg.second,
        ],
    },
    "Current": {"Symbol": "I", "Units": [ureg.ampere, ureg.ampere]},
    "Power": {"Symbol": "W", "Units": [ureg.watt, ureg.megawatt]},
    "Qth": {
        "Symbol": "Qth",
        "Units": [
            ureg.watt / ureg.meter**3,
            ureg.megawatt / ureg.Unit(distance_unit) ** 3,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "temperature": {
        "Symbol": "T",
        "Units": [ureg.degK, ureg.degC],
        "Exclude": ["Air"],
    },
    "U": {
        "Symbol": "V",
        "Units": [ureg.volt, ureg.volt],
        "Exclude": ["Air", "Isolant"],
    },
    "Jth": {
        "Symbol": "J",
        "Units": [
            ureg.ampere / ureg.meter**2,
            ureg.ampere / ureg.Unit(distance_unit) ** 2,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "Fext": {
        "Symbol": "F_ext",
        "mSymbol": r"$F_{ext}$",
        "Units": [
            ureg.newton / ureg.meter**3,
            ureg.newton / ureg.Unit(distance_unit) ** 3,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "Flaplace": {
        "Symbol": "F",
        "Units": [
            ureg.newton / ureg.meter**3,
            ureg.newton / ureg.Unit(distance_unit) ** 3,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "atheta": {
        "Symbol": "A",
        "Units": [
            ureg.ampere / ureg.meter,
            ureg.ampere / ureg.Unit(distance_unit),
        ],
        "Exclude": [],
    },
    "Bg": {
        "Symbol": "B_Bg",
        "mSymbol": r"$B_{bg}$",
        "Units": [ureg.tesla, ureg.tesla],
        "Exclude": [],
    },
    "B": {
        "Symbol": "B",
        "Units": [ureg.tesla, ureg.tesla],
        "Exclude": [],
    },
    "displacement": {
        "Symbol": "u",
        "Units": [
            ureg.meter / ureg.second,
            ureg.Unit(distance_unit) / ureg.second,
        ],
        "Exclude": ["Air"],
    },
    "displacement_r": {
        "Symbol": "ur",
        "Units": [
            ureg.meter / ureg.second,
            ureg.Unit(distance_unit) / ureg.second,
        ],
        "Exclude": ["Air"],
    },
    "displacement_z": {
        "Symbol": "uz",
        "Units": [
            ureg.meter / ureg.second,
            ureg.Unit(distance_unit) / ureg.second,
        ],
        "Exclude": ["Air"],
    },
    "displacementnorm": {
        "Symbol": "u",
        "mSymbol": r"$\| u \|$",
        "Units": [
            ureg.meter / ureg.second,
            ureg.Unit(distance_unit) / ureg.second,
        ],
        "Exclude": ["Air"],
    },
    "Lame1": {
        "Symbol": "Lame1",
        "Units": [
            ureg.pascal / ureg.meter,
            ureg.megapascal / ureg.Unit(distance_unit),
        ],
        "Exclude": ["Air"],
    },
    "Lame2": {
        "Symbol": "Lame2",
        "Units": [
            ureg.pascal / ureg.meter,
            ureg.megapascal / ureg.Unit(distance_unit),
        ],
        "Exclude": ["Air"],
    },
    # "PoissonCoefficient": {}, # no units
    "strain_00": {
        "Symbol": "strain_00",
        "mSymbol": r"$\bar{\bar{\epsilon}}_{00}$",
        "Units": [
            ureg.dimensionless,
            ureg.dimensionless,
        ],
        "Exclude": ["Air"],
    },
    "strain_01": {
        "Symbol": "strain_01",
        "mSymbol": r"$\bar{\bar{\epsilon}}_{01}$",
        "Units": [
            ureg.dimensionless,
            ureg.dimensionless,
        ],
        "Exclude": ["Air"],
    },
    "strain_10": {
        "Symbol": "strain_10",
        "mSymbol": r"$\bar{\bar{\epsilon}}_{10}$",
        "Units": [
            ureg.dimensionless,
            ureg.dimensionless,
        ],
        "Exclude": ["Air"],
    },
    "strain_11": {
        "Symbol": "strain_11",
        "mSymbol": r"$\bar{\bar{\epsilon}}_{11}$",
        "Units": [
            ureg.dimensionless,
            ureg.dimensionless,
        ],
        "Exclude": ["Air"],
    },
    "stress_00": {
        "Symbol": "stress_00",
        "mSymbol": r"$\bar{\bar{\sigma}}_{00}$",
        "Units": [
            ureg.pascal,
            ureg.megapascal,
        ],
        "Exclude": ["Air"],
    },
    "stress_01": {
        "Symbol": "stress_01",
        "mSymbol": r"$\bar{\bar{\sigma}}_{01}$",
        "Units": [
            ureg.pascal,
            ureg.megapascal,
        ],
        "Exclude": ["Air"],
    },
    "stress_10": {
        "Symbol": "stress_10",
        "mSymbol": r"$\bar{\bar{\sigma}}_{10}$",
        "Units": [
            ureg.pascal,
            ureg.megapascal,
        ],
        "Exclude": ["Air"],
    },
    "stress_11": {
        "Symbol": "stress_11",
        "mSymbol": r"$\bar{\bar{\sigma}}_{11}$",
        "Units": [
            ureg.pascal,
            ureg.megapascal,
        ],
        "Exclude": ["Air"],
    },
    "HoopStrain": {
        "Symbol": "HoopStrain",
        "mSymbol": r"$\bar{\bar{\epsilon}}_{Hoop}$",
        "Units": [
            ureg.dimensionless,
            ureg.dimensionless,
        ],
        "Exclude": ["Air"],
    },
    "HoopStress": {
        "Symbol": "HoopStress",
        "mSymbol": r"$\bar{\bar{\sigma}}_{Hoop}$",
        "Units": [
            ureg.pascal,
            ureg.megapascal,
        ],
        "Exclude": ["Air"],
    },
    "Vonmises": {
        "Symbol": "VonMises",
        "mSymbol": r"$\bar{\bar{\sigma}}_{VonMises}$",
        "Units": [
            ureg.pascal,
            ureg.megapascal,
        ],
        "Exclude": ["Air"],
    },
}

ignored_keys = [
    "elasticity.Lame1",
    "elasticity.Lame2",
    "elasticity.PoissonCoefficient",
    "elasticity.YoungModulus",
    "Area",
    "AxiVolume",
    "r",
    "Cos",
    "Sin",
]


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
            print(qtype, quantity, "data=", data)
    elif isinstance(quantity, list):
        data = (
            Quantity(quantity, units[qtype][0]).to(units[qtype][1]).magnitude.tolist()
        )
    else:
        raise Exception(
            f"convert_data/quantity: unsupported type {type(quantity)} for {qtype}"
        )

    return data


def createStatsTable(stats: list, name: str, verbose: bool = False):
    import pandas as pd

    # TODO add a column with the Block Name
    # the column Block Name in cvs is not what I think
    # it is either "Primary Statistics" or "Derived Statistics"
    # remove "Derived Statistics" rows
    # why Standard Deviation and Variance are NaN??
    print(
        f"createStatsTable: aggregate stats by datatype and field ({len(stats)}) items",
        flush=True,
    )
    _dataset = {
        "PointData": {},
        "CellData": {},
        "FieldData": {},
    }

    for statdict in stats:
        for datatype in statdict:
            # print(f"datatype={datatype}", flush=True)
            for key, kdata in statdict[datatype]["Arrays"].items():
                # print(f"key={key} kdata={kdata.keys()}", flush=True)
                if "Stats" in kdata:
                    """
                    print(
                        tabulate(kdata["Stats"], headers="keys", tablefmt="psql"),
                        flush=True,
                    )
                    """

                    if not key in _dataset[datatype]:
                        # print(f"create _dataset[{datatype}][{key}]")
                        _dataset[datatype][key] = []

                    # print(f"append kdata[Stats] to _dataset[{datatype}][{key}]")
                    _dataset[datatype][key].append(kdata["Stats"])
                    """
                    print(
                        f"_dataset[{datatype}][{key}]: len() = {len(_dataset[datatype][key])}"
                    )
                    """

    # print("_dataset:", flush=True)
    dfs = []
    for datatype in _dataset:
        # print(f"datatype={datatype}")
        for key in _dataset[datatype]:
            # print(f"DescriptiveStats for datatype={datatype}, key={key}:", flush=True)
            # print(f"dataset: {len(_dataset[datatype][key])}")

            keyinfo = key.split(".")
            # print(f"keyinfo={keyinfo}", flush=True)
            if len(keyinfo) == 2:
                (physic, fieldname) = keyinfo
            elif len(keyinfo) == 3:
                (toolbox, physic, fieldname) = keyinfo
            else:
                raise RuntimeError(f"{key}: cannot get keyinfo as splitted char")

            # print(f"toolbox={toolbox}", flush=True)
            # print(f"physic={physic}", flush=True)
            # print(f"fieldname={fieldname}", flush=True)

            # Exclude row with Name in fieldunits[fieldname]['Exclude']
            excludeblocks = fieldunits[fieldname]["Exclude"]
            found = False
            for block in excludeblocks:
                if block in name:
                    found = True
                    break

            """
            for dset in _dataset[datatype][key]:
                print(tabulate(dset, headers="keys", tablefmt="psql"))
            """
            if not found:
                df = pd.concat(_dataset[datatype][key])
                # print(
                #     f"Aggregated DescriptiveStats for datatype={datatype}, key={key}:"
                # )

                # Reorder columns
                df = df[
                    [
                        "Variable",
                        "Name",
                        "Minimum",
                        "Mean",
                        "Maximum",
                        "Standard Deviation",
                        "M2",
                        "M3",
                        "M4",
                    ]
                ]

                if verbose:
                    print(
                        tabulate(df, headers="keys", tablefmt="psql", showindex=False)
                    )
                df.to_csv(f"{key}-descriptivestats-create.csv")
                dfs.append(df)

    total_df = pd.concat(dfs)
    print(tabulate(total_df, headers="keys", tablefmt="psql", showindex=False))
    total_df.to_csv(f"{name}-descriptiveAxistats-create.csv")

    pass


# plot with matplotlib
def plotHistoAxi(
    filename: str,
    name: str,
    key: str,
    BinCount: int,
    show: bool = True,
    verbose: bool = False,
):
    print(f"plotHistAxi: name={name}, key={key}, bin={BinCount}", flush=True)

    ax = plt.gca()
    csv = pd.read_csv(filename)
    keys = csv.columns.values.tolist()
    # print(f"plotHistAxi: keys={keys}", flush=True)
    # print("histo before scaling")
    # print(tabulate(csv, headers="keys", tablefmt="psql"))

    # get key unit
    keyinfo = key.replace("_Magnitude", "").split(".")
    # print(f"keyinfo={keyinfo}", flush=True)
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
    # print(f"in_units={in_unit}, out_units={out_unit}", flush=True)

    units = {fieldname: fieldunits[fieldname]["Units"]}
    values = csv[key].to_list()
    out_values = convert_data(units, values, fieldname)
    csv[key] = out_values  # [f"{val:.2f}" for val in out_values]

    counts, extend_bins, patches = hist(csv[key], bins=BinCount, weights=csv["AxiVol"])
    print(f"counts={counts}")
    print(f"extend_bins={extend_bins}")
    print(f"patches={patches}")

    total_key = "Fraction of total Volume [%]"
    plt.xlabel(rf"{msymbol}[{out_unit:~P}]")
    plt.ylabel(total_key)
    plt.xticks(rotation=45, ha="right")
    plt.title(f"{name}: {key}")
    plt.grid(True)
    # plt.legend(False)

    # if legend is mandatory, set legend to True above and comment out the following line
    # ax.legend([rf"{symbol}[{out_unit:~P}]"])
    ax.yaxis.set_major_formatter(lambda x, pos: f"{x:.1f}")
    ax.xaxis.set_major_formatter(lambda x, pos: f"{x:.2f}")
    show = False
    if show:
        plt.show()
    else:
        plt.tight_layout()
        plt.savefig(
            f'{name}-{key.replace("_Magnitude", "")}-histogram-matplotlib.png', dpi=300
        )
    plt.close()

    pass


def selectBlocks(blockdata: list, excludes: list):
    """
    select Block name not in exludes

    blockdata: list of block names
    excludes: list of names to be excluded

    """

    selectedblocks = []
    for key in blockdata:
        found = False
        for excluded in excludes:
            if excluded in key:
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
    print(f"Memory: {dataInfo.DataInformation.GetMemorySize()}")
    print(f"Bounds: {dataInfo.DataInformation.GetBounds()}")
    print(f"Nodes: {dataInfo.GetNumberOfPoints()}")
    print(f"Cells: {dataInfo.GetNumberOfCells()}")
    print(f"pointData: {input.PointData[:]}", flush=True)
    print(f"cellData: {input.CellData[:]}")
    print(f"fieldData: {input.FieldData[:]}")
    return dataInfo


def resultinfo(input, verbose: bool = False):
    """
    returns a dict gathering info on PointData, CellData and FieldData
    """
    print("resultinfo:", flush=True)

    datadict = {
        "PointData": {"TypeMode": "POINT", "AttributeMode": "Point Data", "Arrays": {}},
        "CellData": {"TypeMode": "CELL", "AttributeMode": "Cell Data", "Arrays": {}},
        "FieldData": {"TypeMode": "FIELD", "AttributeMode": "Field Data", "Arrays": {}},
    }

    if verbose:
        print("PointData:", flush=True)
    for key in input.PointData:
        datadict["PointData"]["Arrays"][key.Name] = getresultInfo(key, verbose)
    if verbose:
        print("CellData:", flush=True)
    for key in input.CellData:
        datadict["CellData"]["Arrays"][key.Name] = getresultInfo(key, verbose)
    if verbose:
        print("FieldData:", flush=True)
    for key in input.FieldData:
        datadict["FieldData"]["Arrays"][key.Name] = getresultInfo(key, verbose)

    return datadict


def getresultInfo(key, verbose: bool = False, printed: bool = True):
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
    if verbose:
        print(f"getresultInfo {name}: datadict={datadict}", flush=True)

    return datadict


def resultStats(input, name: str, Area: float, verbose: bool = False):
    """
    compute stats for PointData, CellData and FieldData

    returns a dict
    """

    datadict = resultinfo(input, verbose)
    print(
        f'resultStats[{name}]: datadict={datadict["PointData"]["Arrays"].keys()}',
        flush=True,
    )
    PointData_keys = list(datadict["PointData"]["Arrays"].keys())

    spreadSheetView = CreateView("SpreadSheetView")
    cellCenters1Display = Show(input, spreadSheetView, "SpreadSheetRepresentation")
    spreadSheetView.Update()

    filename = f"{name}-Axi-cellcenters-all.csv"
    ExportView(
        filename,
        view=spreadSheetView,
        RealNumberNotation="Scientific",
    )

    csv = pd.read_csv(filename)
    keys = csv.columns.values.tolist()
    print(f"read_csv: csv={filename}, keys={keys}", flush=True)

    PointData_keys = []
    for item in input.PointData[:]:
        suffix = ""
        if item.GetNumberOfComponents() > 1:
            suffix = "_Magnitude"

            for component in range(item.GetNumberOfComponents()):
                ignore = f"{item.Name}_{component}"
                if ignore not in ignored_keys:
                    ignored_keys.append(f"{item.Name}_{component}")

        PointData_keys.append(f"{item.Name}{suffix}")

    # print(f"PointData_keys={PointData_keys}")
    # print(f"csv keys={keys}")
    missing_keys = [key for key in PointData_keys if key not in keys]
    if missing_keys:
        print(f"missing_keys={missing_keys}")
        exit(1)

    for datatype in datadict:
        if datatype == "PointData":
            AttributeMode = datadict[datatype]["AttributeMode"]
            TypeMode = datadict[datatype]["TypeMode"]
            for key, kdata in datadict[datatype]["Arrays"].items():
                if not key in ignored_keys:
                    Components = kdata["Components"]
                    bounds = kdata["Bounds"]

                    keyinfo = key.split(".")
                    # print(f"keyinfo={keyinfo}", flush=True)
                    if len(keyinfo) == 2:
                        (physic, fieldname) = keyinfo
                    elif len(keyinfo) == 3:
                        (toolbox, physic, fieldname) = keyinfo
                    else:
                        raise RuntimeError(
                            f"{key}: cannot get keyinfo as splitted char"
                        )

                    symbol = fieldunits[fieldname]["Symbol"]
                    msymbol = symbol
                    if "mSymbol" in fieldunits[fieldname]:
                        msymbol = fieldunits[fieldname]["mSymbol"]
                    units = {fieldname: fieldunits[fieldname]["Units"]}
                    [in_unit, out_unit] = fieldunits[fieldname]["Units"]

                    if bounds[0][0] != bounds[0][1]:
                        stats = {
                            "Variable": [rf"{symbol}[{out_unit:~P}]"],
                            "Name": [name],
                            "Minimum": [convert_data(units, bounds[0][0], fieldname)],
                            "M1": [0],
                            "Maximum": [convert_data(units, bounds[0][1], fieldname)],
                            "Standard Deviation": [0],
                            "M2": [0],
                            "M3": [0],
                            "M4": [0],
                        }

                        # print(f"AxiStats for {key}", flush=True)
                        if Components > 1:
                            key = f"{key}_Magnitude"
                        for order in range(1, 5):
                            # print(f"Create moment{order} for {key}", flush=True)
                            units = {f"M{order}": fieldunits[fieldname]["Units"]}
                            if in_unit != ureg.kelvin and order > 1:
                                units[f"M{order}"] = [in_unit**order, out_unit**order]
                            else:
                                units[f"M{order}"] = [in_unit**order, in_unit**order]

                            res = (
                                2
                                * pi
                                * csv["Area"]
                                * csv["Points_0"]
                                * csv[key] ** order
                            )
                            value = res.sum() / Area
                            # print(
                            #     f"{key} M{order}: {value}, integ={res.sum()}, Area={Area}, name={name}")
                            #    type={type(value)}, type={type(value.item())}
                            out_res = convert_data(units, value.item(), f"M{order}")
                            # print(f"M{order}: {out_res}")
                            stats[f"M{order}"] = [out_res]

                        stats["Standard Deviation"] = [
                            sqrt(abs(stats["M1"][0] ** 2 - stats["M2"][0]))
                        ]
                        # rename M1 to Mean
                        stats["Mean"] = stats["M1"]
                        del stats["M1"]
                        # print(f"{key}: stats={stats}", flush=True)

                        kdata["Stats"] = pd.DataFrame.from_dict(stats)

                        if verbose:
                            print(f"resultStats: key={key}")
                            print(
                                tabulate(
                                    kdata["Stats"],
                                    headers="keys",
                                    tablefmt="psql",
                                )
                            )

    # remove: f"{name}-Axi-cellcenters-all.csv"
    os.remove(filename)

    # display stats
    return datadict


def resultHistos(
    input,
    name: str,
    Area: float,
    BinCount: int = 10,
    printed: bool = True,
    show: bool = False,
    verbose: bool = False,
):
    """
    histogram
    """
    print(f"resultHistos: name={name}, Area={Area}, BinCount={BinCount}", flush=True)

    spreadSheetView = CreateView("SpreadSheetView")
    cellCenters1Display = Show(input, spreadSheetView, "SpreadSheetRepresentation")
    spreadSheetView.Update()

    filename = f"{name}-Axi-cellcenters-all.csv"
    ExportView(
        filename,
        view=spreadSheetView,
        RealNumberNotation="Scientific",
    )

    csv = pd.read_csv(filename)
    keys = csv.columns.values.tolist()
    print(f"read_csv: csv={filename}, keys={keys}", flush=True)

    sum = csv["AxiVolume"].sum()
    csv["AxiVol"] = csv["AxiVolume"] / sum * 100
    print(f"Area={Area}, sum={sum}")

    # check that sum is roughtly equal to 1
    # print(f'check Sum(Fraction): {csv["Fraction of total Area [%]"].sum()}')
    eps = 1.0e-4
    error = abs(1 - csv["AxiVol"].sum() / 100.0)
    assert error <= eps, f"Check Sum(Fraction) failed (error={error} > eps={eps})"

    csv.to_csv(filename)
    print(f'Sum(AxiVol)={csv["AxiVol"].sum()}')

    datadict = resultinfo(input)
    for datatype in datadict:
        if datatype == "CellData":
            AttributeMode = datadict[datatype]["AttributeMode"]
            TypeMode = datadict[datatype]["TypeMode"]
            for key, kdata in datadict[datatype]["Arrays"].items():
                if not key in ignored_keys:
                    Components = kdata["Components"]
                    bounds = kdata["Bounds"]
                    if bounds[0][0] != bounds[0][1]:
                        keyname = key
                        if Components > 1:
                            keyname = f"{key}_Magnitude"
                        plotHistoAxi(
                            filename,
                            name,
                            keyname,
                            BinCount,
                            show=show,
                            verbose=verbose,
                        )

    Delete(spreadSheetView)
    del spreadSheetView

    # Force a garbage collection
    collected = gc.collect()
    print(f"resultsHistos: Garbage collector: collected {collected} objects.")

    # remove: f"{name}-Axi-cellcenters-all.csv"
    os.remove(filename)


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

    print(f"Load Ensight case: {file}")
    input = OpenDataFile(file)
    UpdatePipeline()

    if not printed:
        print(f"properties: {input.ListProperties()}")
        for prop in input.ListProperties():
            print(f"reader: {prop}={input.GetPropertyValue(prop)}")

    return input


def momentN(input, key: str, nkey: str, order: int, AttributeType: str):
    """
    compute moment of order N
    """

    calculator1 = Calculator(registrationName=f"moment{order}", Input=input)
    calculator1.AttributeType = AttributeType  # 'Cell Data'
    calculator1.ResultArrayName = f"{nkey}_moment{order}"

    # check Points_0 for r
    if order == 1:
        calculator1.Function = f"coordsX*{key}"
    else:
        calculator1.Function = f"coordsX*{key}^{order}"

    calculator1.UpdatePipeline()
    return calculator1


def integrateKeys(input, name: str, printed: bool = True):
    """
    compute integral of Data over area/volume

    to get values use a spreadsheet
    """

    integratedvalues = IntegrateVariables(Input=input)
    if not printed:
        for prop in integratedvalues.ListProperties():
            print(f"integratedvalues: {prop}={integratedvalues.GetPropertyValue(prop)}")

    spreadSheetView = CreateView("SpreadSheetView")
    descriptiveDisplay = Show(
        integratedvalues, spreadSheetView, "SpreadSheetRepresentation"
    )
    spreadSheetView.Update()
    export = ExportView(
        f"{name}-integrals.csv",
        view=spreadSheetView,
        RealNumberNotation="Scientific",
    )

    Delete(spreadSheetView)
    del spreadSheetView
    Delete(descriptiveDisplay)

    # Force a garbage collection
    collected = gc.collect()
    print(f"integrateKeys: Garbage collector: collected {collected} objects.")

    return integratedvalues


def meshinfo(
    input,
    ComputeStats: bool = True,
    ComputeHisto: bool = False,
    BinCount: int = 10,
    show: bool = False,
    verbose: bool = False,
    printed: bool = True,
):
    """
    display geometric info from input dataset
    """

    def scaleField(input, key: str, nkey: str, AttributeType: str, factor: float):
        """
        scale a Field by 1/factor

        AttributeType: 'Cell Data'
        """
        calculator1 = Calculator(registrationName="Calculator1", Input=input)
        calculator1.AttributeType = AttributeType  # 'Cell Data'
        calculator1.ResultArrayName = f"{nkey}"
        calculator1.Function = f'{key}/{factor}")'

        calculator1.UpdatePipeline()
        return calculator1

    def addField(input, key: str, nkey: str, AttributeType: str, factor: float):
        """
        add factor to Field

        AttributeType: 'Cell Data'
        """
        calculator1 = Calculator(registrationName="Calculator1", Input=input)
        calculator1.AttributeType = AttributeType  # 'Cell Data'
        calculator1.ResultArrayName = f"{nkey}"

        # if float is negatif
        if float < 0:
            calculator1.Function = f'{key}-{abs(factor)}")'
        else:
            calculator1.Function = f'{key}+{factor}")'

        calculator1.UpdatePipeline()
        return calculator1

    def createVectorNorm(
        input, key: str, nkey: str, AttributeType: str, printed: bool = True
    ):
        """
        create Norm of a vector

        AttributeType: 'Cell Data'
        """

        calculator1 = Calculator(registrationName="Calculator1", Input=input)
        calculator1.AttributeType = AttributeType  # 'Cell Data'
        calculator1.ResultArrayName = f"{key}norm"
        calculator1.Function = f'mag("{key}")'
        if not printed:
            for prop in calculator1.ListProperties():
                print(
                    f"Calculator: {prop}={calculator1.GetPropertyValue(prop)}",
                    flush=True,
                )

        calculator1.UpdatePipeline()
        return calculator1

    """
    # does not apply for CellData vector (no coordX)
    for field in input.PointData:
        if field.GetNumberOfComponents() == 3:
            print(
                f"create {field.Name}norm for {field.Name} PointData vector",
                flush=True,
            )
            tmp = createVectorNorm(
                tmp, field.Name, field.Name, "Point Data"
            )
    for field in input.CellData:
        if field.GetNumberOfComponents() == 3:
            print(f"create norm for {field.Name} CellData vector", flush=True)
            tmp = createVectorNorm(
                tmp, field.Name, field.Name, "Cell Data"
            )
    """

    def part_integrate(
        input, name, selected_blocks: list, merge: bool = True, verbose: bool = False
    ):
        """
        compute integral over input
        """
        print(f"part_integrate: name={name}")

        # try to add Moment here
        # convert CellData to PointData
        cellDatatoPointData = CellDatatoPointData(
            registrationName="CellDatatoPointData", Input=input
        )
        tmp = cellDatatoPointData

        if selected_blocks:
            extractBlock1 = ExtractBlock(registrationName="insert", Input=tmp)
            extractBlock1.Selectors = selected_blocks
            extractBlock1.UpdatePipeline()
            tmp = extractBlock1

        if merge:
            mergeBlocks1 = MergeBlocks(registrationName="MergeBlocks1", Input=tmp)
            mergeBlocks1.UpdatePipeline()
            tmp = mergeBlocks1

        drop_keys = []
        for field in input.PointData:
            fname = f'"{field.Name}"'
            if field.GetNumberOfComponents() > 1:
                fname = f'mag("{field.Name}")'
                for i in range(field.GetNumberOfComponents()):
                    drop_keys.append(f"{field.Name}_{i}")
            for order in range(1, 5):
                tmp = momentN(tmp, fname, field.Name, order, "Point Data")

        integratedvalues = integrateKeys(tmp, name, printed=False)
        print(f"IntegratedValues: PointData={integratedvalues.PointData[:]}")
        print(f"IntegratedValues: FieldData={integratedvalues.FieldData[:]}")
        drop_keys += ["Point ID", "Points_1", "Points_2", "Points_Magnitude"]
        csv = pd.read_csv(f"{name}-integrals.csv")
        csv.drop(columns=drop_keys, inplace=True)

        if verbose:
            print(
                tabulate(csv.transpose(), headers="keys", tablefmt="psql"),
                flush=True,
            )

        Delete(cellDatatoPointData)
        del cellDatatoPointData
        return csv

    # PointData to CellData
    pointDatatoCellData = PointDatatoCellData(
        registrationName="PointDatatoCellData", Input=input
    )

    print("Get mesh size")
    cellsize = CellSize(pointDatatoCellData)  # input

    # set some params
    cellsize.ComputeLength = 0
    cellsize.ComputeArea = 1
    cellsize.ComputeVolume = 0
    cellsize.ComputeVertexCount = 0
    cellsize.ComputeSum = 1
    # get params list
    if not printed:
        for prop in cellsize.ListProperties():
            print(f"cellsize: {prop}={cellsize.GetPropertyValue(prop)}")

    # apply
    cellsize.UpdatePipeline()

    cellcenters = CellCenters(registrationName="CellCenters", Input=cellsize)
    # Properties modified on cellCenters1
    cellcenters.VertexCells = 1
    cellcenters.UpdatePipeline()

    # Add Volume
    calculator1 = Calculator(registrationName=f"AxiVolume", Input=cellcenters)
    calculator1.AttributeType = "Point Data"
    calculator1.ResultArrayName = f"AxiVolume"
    calculator1.Function = f"2*{pi}*coordsX*Area"
    calculator1.UpdatePipeline()

    dataInfo = info(calculator1)
    dataset = sm.Fetch(calculator1)
    np_dataset = dsa.WrapDataObject(dataset)

    blockdata = {}
    # check dataset type
    print(f"type(dataset)={type(dataset)}", flush=True)
    if dataset.IsA("vtkUnstructuredGrid"):
        print("UnstructuredGrid")
        block = dataset
    elif dataset.IsA("vtkMultiBlockDataSet"):
        print("MultiBlockDataSet")
        block = dataset.GetBlock(0)

    Area_data = block.GetPointData().GetArray("AxiVolume")
    # print(f'block cellData[Area_data]: {Area_data}, type={type(Area_data)}')
    Area = block.GetPointData().GetArray("AxiVolume")
    # print(f'block fieldData[Area]: {Area}, type={type(Area)}')
    np_Area = np_dataset.PointData["AxiVolume"]
    # print(f'block fieldData[np_Area]: {np_Area}, type={type(np_Area)}, length={algs.shape(np_Area)}')
    tvol = algs.sum(np_Area)
    vunits = fieldunits["Volume"]["Units"]
    mm3 = f"{vunits[1]:~P}"
    tvol_mm3 = convert_data(
        {"Volume": vunits},
        tvol,
        "Volume",
    )
    print(
        f"block fieldData[np_Area]: total={tvol_mm3} {mm3}, parts={algs.shape(np_Area)}"
    )

    if dataset.IsA("vtkMultiBlockDataSet"):
        hierarchy = dataInfo.GetHierarchy()
        rootnode = hierarchy.GetRootNode()
        rootSelector = f"/{hierarchy.GetRootNodeName()}"
        blocks = hierarchy.GetNumberOfChildren(rootnode)
        print(f"Load blocks: {blocks}")

        sum_vol = 0
        for i in range(blocks):
            child = hierarchy.GetChild(rootnode, i)
            name = hierarchy.GetNodeName(child)
            rootChild = f"{rootSelector}/{name}"
            child_info = reader.GetSubsetDataInformation(0, child)  # rootChild)
            # print(f'block[{i}]: {name}, child={child}, childInfo: {child_info}')

            nodes = child_info.GetNumberOfPoints()
            cells = child_info.GetNumberOfCells()
            bounds = child_info.GetBounds()

            vol = np_Area.Arrays[i][0]
            vol_mm3 = convert_data(
                {"Volume": vunits},
                vol,
                "Volume",
            )

            # do not print volumes since they are wrong for unknown reason
            print(
                f"block[{i}]: {name}, nodes={nodes}, cells={cells}, bounds={bounds}"
            )  # vol={vol_mm3} {mm3},

            blockdata[rootChild] = {
                "name": name,
                "bounds": bounds,
                "nodes": nodes,
                "cells": cells,
                "Area": vol,
            }

            sum_vol += vol

        # # check tvol == Sum(vol)
        # if abs(1 - sum_vol / tvol) > 1.0e-3:
        #     print(f"Total volume != Sum(vol), tvol={tvol}, sum_vol={sum_vol}, error={abs(1-sum_vol/tvol)*100} %")

        # Compute Stats
        stats = []

        print("Data ranges:", flush=True)
        resultinfo(calculator1, verbose)
        Delete(pointDatatoCellData)
        del pointDatatoCellData

        # Force a garbage collection
        collected = gc.collect()
        print(f"meshinfo: Garbage collector: collected {collected} objects.")

        print("Data ranges without Air:", flush=True)
        selected_blocks = [block for block in blockdata.keys() if not "Air" in block]

        # TODO: need to restart from input and apply ... Calculator1
        extractBlock1 = ExtractBlock(registrationName="insert", Input=input)
        extractBlock1.Selectors = selected_blocks
        extractBlock1.UpdatePipeline()
        mergeBlocks1 = MergeBlocks(registrationName="MergeBlocks1", Input=extractBlock1)
        mergeBlocks1.UpdatePipeline()

        def part(
            pinput,
            name: str,
            computeHisto: bool,
            BinCount: int = 20,
            show: bool = False,
            verbose: bool = False,
        ):
            pointDatatoCellData = PointDatatoCellData(pinput)
            cellsize = CellSize(pointDatatoCellData)
            # set some params
            cellsize.ComputeLength = 0
            cellsize.ComputeArea = 1
            cellsize.ComputeVolume = 0
            cellsize.ComputeVertexCount = 0
            cellsize.ComputeSum = 1
            cellsize.UpdatePipeline()

            cellcenters = CellCenters(registrationName="CellCenters", Input=cellsize)
            # Properties modified on cellCenters1
            cellcenters.VertexCells = 1
            cellcenters.UpdatePipeline()

            # Add Volume
            calculator1 = Calculator(registrationName=f"AxiVolume", Input=cellcenters)
            calculator1.AttributeType = "Point Data"
            calculator1.ResultArrayName = f"AxiVolume"
            calculator1.Function = f"2*{pi}*coordsX*Area"
            calculator1.UpdatePipeline()

            dataset = sm.Fetch(calculator1)
            np_dataset = dsa.WrapDataObject(dataset)
            print(f"type(dataset)={type(dataset)}", flush=True)
            if dataset.IsA("vtkUnstructuredGrid"):
                print("UnstructuredGrid", flush=True)
                block = dataset  # .GetBlock(0)
            elif dataset.IsA("vtkPolyData"):
                print("vtkPolyData", flush=True)
                block = dataset  # .GetBlock(0)
            elif dataset.IsA("vtkMultiBlockDataSet"):
                print("MultiBlockDataSet", flush=True)
                block = dataset.GetBlock(0)

            Area_data = block.GetPointData().GetArray("AxiVolume")
            # print(f'block cellData[Area_data]: {Area_data}, type={type(Area_data)}')
            Area = block.GetPointData().GetArray("AxiVolume")
            # print(f'block fieldData[Area]: {Area}, type={type(Area)}')
            np_Area = np_dataset.PointData["AxiVolume"]
            # print(f'block fieldData[np_Area]: {np_Area}, type={type(np_Area)}, length={algs.shape(np_Area)}')
            vol = algs.sum(np_Area)
            vunits = fieldunits["Volume"]["Units"]
            mm3 = f"{vunits[1]:~P}"
            vol_mm3 = convert_data(
                {"Volume": vunits},
                vol,
                "Volume",
            )
            print(
                f"{name}: block fieldData[np_Area]: vol={vol_mm3} {mm3}, parts={algs.shape(np_Area)}"
            )

            # # check tvol == Sum(vol)
            # if abs(1 - vol / tvol) > 1.0e-3:
            #     print(f"insert Total volume != vol(insert), tvol={tvol}, vol={vol}, error={abs(1-vol/tvol)*100} %")

            statsdict = resultStats(calculator1, name, tvol, verbose)
            # print(f"insert statsdict: {statsdict}", flush=True)
            if ComputeHisto:
                resultHistos(
                    calculator1,
                    name,
                    vol,
                    BinCount=BinCount,
                    show=show,
                    verbose=verbose,
                )

            stats.append(statsdict)

            Delete(pointDatatoCellData)
            del pointDatatoCellData

            # Force a garbage collection
            collected = gc.collect()
            print(f"part: Garbage collector: collected {collected} objects.")

            return vol, statsdict

        vol, statsdict = part(mergeBlocks1, "insert", ComputeHisto, BinCount)
        stats.append(statsdict)

        icsv = part_integrate(
            input, "insert", selected_blocks, merge=True, verbose=True
        )
        print(f'insert: vol={vol}, ivol={icsv["Points_0"].to_list()[0] * 2 * pi}')
        # aggregate stats data
        createStatsTable([statsdict], "insert")

        if not ComputeStats:
            return cellsize, blockdata, statsdict

        print("Data ranges per block:", flush=True)
        sum_vol = 0
        for i, block in enumerate(blockdata.keys()):
            name = blockdata[block]["name"]
            print(f"block[{i}]: extract {block}, name={name}", flush=True)
            # TODO: need to restart from input and apply ... Calculator1
            extractBlock1 = ExtractBlock(registrationName=name, Input=input)
            extractBlock1.Selectors = [block]
            extractBlock1.UpdatePipeline()

            vol, statsdict = part(extractBlock1, name, ComputeHisto, BinCount)
            sum_vol += vol

            icsv = part_integrate(input, name, [block], merge=False, verbose=True)
            print(
                f'name={name}: vol={vol}, ivol={icsv["Points_0"].to_list()[0] * 2 * pi}'
            )

            stats.append(statsdict)
            Delete(extractBlock1)
            del extractBlock1

            # Force a garbage collection
            collected = gc.collect()
            print(f"loopblock: Garbage collector: collected {collected} objects.")

            # aggregate stats data
            createStatsTable([statsdict], name)

        # check tvol == Sum(vol)
        print(
            f"reworked Total volume: tvol={tvol}, sum_vol={sum_vol}, error={abs(1-sum_vol/tvol)*100} %"
        )
        if abs(1 - sum_vol / tvol) > 1.0e-3:
            print(
                f"reworked Total volume != Sum(vol), tvol={tvol}, sum_vol={sum_vol}, error={abs(1-sum_vol/tvol)*100} %"
            )

        # aggregate stats data
        createStatsTable(stats, "total")

    return input, blockdata, dict()


def deformed(input, factor: float = 1, printed: bool = True):
    """
    create deformed geometry
    """
    print(f"deformed view: factor={factor}", flush=True)

    warpByVector1 = WarpByVector(registrationName="WarpByVector1", Input=input)
    warpByVector1.Vectors = ["POINTS", "elasticity.displacement"]
    warpByVector1.ScaleFactor = factor
    # get params list
    if not printed:
        for prop in warpByVector1.ListProperties():
            print(f"warpByVector1: {prop}={warpByVector1.GetPropertyValue(prop)}")

    UpdatePipeline()
    return warpByVector1


################################################################
def displayField(
    input,
    selectedblocks: list,
    field: str,
    color,
    addruler: bool = True,
    renderView=None,
    filename: str = None,
    comment: str = None,
    grid: bool = False,
    printed: bool = True,
):
    """
    display field in renderview

    TODO: eventually add an annotation
    """

    print(f"displayField: field={field}, renderView={renderView}", flush=True)
    if renderView is None:
        renderView = CreateView("RenderView")
        print("createrenderView", flush=True)

    if comment is not None:
        print(f"add comment: {comment}", flush=True)
        text = Text(registrationName="Text1")
        text.Text = rf"{comment}"
        textDisplay = Show(text, renderView)
        textDisplay.WindowLocation = "Upper Center"
        textDisplay.FontSize = 24
        textDisplay.Bold = 1
        textDisplay.Italic = 1
        textDisplay.Color = [0.0, 0.0, 0.0]

    display = Show(input, renderView)

    # TODO if args.field is not Magnetic something
    display.BlockSelectors = selectedblocks
    display.ColorArrayName = color
    if grid:
        display.DataAxesGrid.GridAxesVisibility = 1

    # for vector: ColorBy(display, ('CELLS', 'magnetic_field', 'Z'))
    ColorBy(display, tuple(color))

    renderView.ResetCamera()
    renderView.GetActiveCamera()

    # Add BoundingRuler filter to get an idea of the dimension
    if addruler:
        print("Add ruler to see dimensions", flush=True)
        ruler = BoundingRuler(registrationName="BoundingRuler1", Input=input)
        ruler.Axis = "Y Axis"
        if not printed:
            for prop in ruler.ListProperties():
                print(f"ruler: {prop}={ruler.GetPropertyValue(prop)}")
        Show(ruler, renderView)  # Reset Camera

    resolution = [1200, 1200]
    renderView.ViewSize = resolution
    if not printed:
        for prop in renderView.ListProperties():
            print(f"renderView: {prop}={renderView.GetPropertyValue(prop)}")
    renderView.OrientationAxesVisibility = 1
    display.SetScalarBarVisibility(renderView, True)
    display.RescaleTransferFunctionToDataRange(True, False)

    # get color transfer function/color map for 'thermo_electricheattemperature'
    field_name = field.replace(".", "")
    LUT = GetColorTransferFunction(field_name)
    LUT.ScalarRangeInitialized = 1.0

    # LUT.RescaleTransferFunction(293.6058044433594, 397.88848876953125)
    # valid range for temperature but where do the range come from?
    # see post https://stackoverflow.com/questions/63028755/paraview-rescaling-colour-scheme-to-visible-data-in-range-in-python

    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    # LUT.ApplyPreset("Turbo", True)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(LUT, renderView)

    # get color legend/bar for LUT in view renderView1
    LUTColorBar = GetScalarBar(LUT)
    # LUTColorBar.Position = [0.9118075801749271, 0.01059135039717564]

    keyinfo = field.split(".")
    # print(f"keyinfo={keyinfo}", flush=True)
    if len(keyinfo) == 2:
        (physic, fieldname) = keyinfo
    elif len(keyinfo) == 3:
        (toolbox, physic, fieldname) = keyinfo
    else:
        raise RuntimeError(f"{field}: cannot get keyinfo as splitted char")
    symbol = fieldunits[fieldname]["Symbol"]
    msymbol = symbol
    if "mSymbol" in fieldunits[fieldname]:
        msymbol = fieldunits[fieldname]["mSymbol"]
    [in_unit, out_unit] = fieldunits[fieldname]["Units"]
    LUTColorBar.Title = rf"{msymbol} [{in_unit:~P}]"
    print(f"LUTColorBar.ComponentTitle={LUTColorBar.ComponentTitle}")
    # LUTColorBar.ComponentTitle = ''

    # Properties modified on LUTColorBar
    LUTColorBar.TitleColor = [0.0, 0.0, 0.0]
    # LUTColorBar.TitleFontFamily = 'Courier'
    LUTColorBar.TitleBold = 1
    LUTColorBar.TitleItalic = 1
    # LUTColorBar.TitleShadow = 1
    LUTColorBar.TitleFontSize = 24
    LUTColorBar.HorizontalTitle = 1
    LUTColorBar.LabelColor = [0.0, 0.0, 0.0]
    # LUTColorBar.LabelFontFamily = 'Courier'
    # LUTColorBar.LabelBold = 1
    # LUTColorBar.LabelItalic = 1
    # LUTColorBar.LabelShadow = 1
    LUTColorBar.LabelFontSize = 20
    LUTColorBar.ScalarBarThickness = 32
    LUTColorBar.ScalarBarLength = 0.5
    LUTColorBar.AutomaticLabelFormat = 0
    # LUTColorBar.LabelFormat = '%-#6.2g'
    # LUTColorBar.RangeLabelFormat = '%-#6.1f'
    # LUTColorBar.DataRangeLabelFormat = '%-#5.2f'
    # LUTColorBar.DrawDataRange = 1
    # LUTColorBar.AddRangeAnnotations = 1
    # LUTColorBar.DrawAnnotations = 0
    # LUTColorBar.TextPosition = "Ticks left/bottom, annotations right/top"

    # get opacity transfer function/opacity map
    PWF = GetOpacityTransferFunction(field_name)
    PWF.ScalarRangeInitialized = 1
    # PWF = RescaleTransferFunction(293.6058044433594, 397.88848876953125)
    renderView.Update()

    # save screenshot
    # TransparentBackground=1, need to set title and label colors to black
    if filename:
        save = SaveScreenshot(
            filename,
            renderView,
            ImageResolution=resolution,
            TransparentBackground=1,
        )

    if comment is not None:
        Delete(textDisplay)
    Delete(display)
    return renderView


################################################################
# create a 3D view
def makeAxiview(
    input,
    blockdata,
    field: str,
    color,
    suffix: str = None,
    addruler: bool = False,
    printed: bool = True,
):
    """
    create a 2D view
    """
    print(f"makeAxiview: field={field}", end="")
    if suffix:
        print(f", suffix={suffix}", end="")
    print(flush=True)
    print(f"blockdata={blockdata}", flush=True)

    keyinfo = field.split(".")
    # print(f"keyinfo={keyinfo}", flush=True)
    if len(keyinfo) == 2:
        (physic, fieldname) = keyinfo
    elif len(keyinfo) == 3:
        (toolbox, physic, fieldname) = keyinfo
    else:
        raise RuntimeError(f"{field}: cannot get keyinfo as splitted char")
    print(f"Exclude blocks = {fieldunits[fieldname]['Exclude']}", flush=True)

    selectedblocks = selectBlocks(
        list(blockdata.keys()), fieldunits[fieldname]["Exclude"]
    )
    if selectedblocks:
        print(f"input.Selectors = {selectedblocks}", flush=True)

    filename = f"{field}.png"
    if suffix is not None:
        filename = f"{field}{suffix}.png"

    renderView = displayField(
        input, selectedblocks, field, color, addruler=addruler, filename=filename
    )

    Delete(renderView)
    del renderView


#################################################################


def plotOr(
    input,
    r: list[float],
    z: float,
    show: bool = True,
    printed: bool = True,
):
    from math import pi, cos, sin

    [r0, r1] = r

    plotOverLine = PlotOverLine(registrationName="Oz", Input=input, Source="Line")
    # get params list
    if not printed:
        for prop in plotOverLine.ListProperties():
            print(f"plotOverLine': {prop}={plotOverLine.GetPropertyValue(prop)}")

    # init the 'Line' selected for 'Source'
    plotOverLine.Source.Point1 = [r0, z, 0]
    plotOverLine.Source.Point2 = [r1, z, 0]

    filename = f"r0={r0}m-r1={r1}m-z={z}m.csv"
    export = CreateWriter(filename, proxy=plotOverLine)
    if not printed:
        for prop in export.ListProperties():
            print(f"export: {prop}={export.GetPropertyValue(prop)}")
    export.UpdateVTKObjects()  # is it needed?
    export.UpdatePipeline()

    # plot with matplotlib
    def plotOrField(file, key: str, z: float, ax=None, show: bool = True):
        import pandas as pd
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator

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
            plt.savefig(f"{key}-vs-r-z={z}.png", dpi=300)
        plt.close()

    # requirements: create PointData from CellData
    for field in input.PointData:
        plotOrField(
            filename,
            field,
            z,
            ax=None,
            show=show,
        )
    pass


def plotOz(
    input,
    r: float,
    z: list[float],
    show: bool = True,
    printed: bool = True,
):

    from math import pi, cos, sin

    [z0, z1] = z

    plotOverLine = PlotOverLine(registrationName="Oz", Input=input, Source="Line")

    # init the 'Line' selected for 'Source'
    plotOverLine.Source.Point1 = [r, z0, 0]
    plotOverLine.Source.Point2 = [r, z1, 0]

    filename = f"r={r}m-z0={z0}m-z1={z1}m.csv"
    export = CreateWriter(filename, proxy=plotOverLine)
    if not printed:
        for prop in export.ListProperties():
            print(f"export: {prop}={export.GetPropertyValue(prop)}")
    export.UpdateVTKObjects()  # is it needed?
    export.UpdatePipeline()

    # plot with matplotlib
    def plotOzField(file, key: str, z: float, ax=None, show: bool = True):
        import pandas as pd
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator

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
            plt.savefig(f"{key}-vs-z-r={r}.png", dpi=300)
        plt.close()

    # requirements: create PointData from CellData
    for field in input.PointData:
        plotOzField(
            filename,
            field,
            r,
            ax=None,
            show=show,
        )
    pass


parser = argparse.ArgumentParser()
parser.add_argument("file", type=str, help="input case file (ex. Export.case)")
parser.add_argument("--field", type=str, help="select field to display", default="")
parser.add_argument("--stats", help="activate stats calculations", action="store_true")
parser.add_argument(
    "--histos", help="activate histograms calculations", action="store_true"
)
parser.add_argument("--bins", type=int, help="set bins number (default 10)", default=10)
parser.add_argument("--plots", help="activate plots calculations", action="store_true")
parser.add_argument("--views", help="activate views calculations", action="store_true")
parser.add_argument(
    "--channels", help="activate views calculations", action="store_true"
)
parser.add_argument("--r", nargs="*", type=float, help="select r in m to display")
parser.add_argument("--z", nargs="*", type=float, help="select z in m to display")
parser.add_argument(
    "--save",
    help="save graphs",
    action="store_true",
)
parser.add_argument(
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
args = parser.parse_args()
print(f"args: {args}")

# check paraview version
version = GetParaViewVersion()
print(f"Paraview version: {version}")

cwd = os.getcwd()

# args.file = "../../HL-31/test/hybride-Bh27.7T-Bb9.15T-Bs9.05T_HPfixed_BPfree/bmap/np_32/thermo-electric.exports/Export.case"
reader = load(args.file)
bounds = getbounds(reader)
print(f"bounds={bounds}")  # , type={type(bounds)}")
info(reader)
color = []
CellToData = False
if args.field:
    if args.field in list(reader.CellData.keys()):
        field = reader.CellData[args.field]
        # if field.GetNumberOfComponents() == 1:
        color = ["CELLS", args.field]
        # if field.GetNumberOfComponents() == 3:
        #    color = ["CELLS", args.field, "Magnitude"]
    if args.field in list(reader.PointData.keys()):
        field = reader.PointData[args.field]
        # if field.GetNumberOfComponents() == 1:
        color = ["POINTS", args.field]
        # if field.GetNumberOfComponents() == 3:
        #    color = ["POINTS", args.field, "Magnitude"]


# get Block info
cellsize, blockdata, statsdict = meshinfo(
    reader,
    ComputeStats=args.stats,
    ComputeHisto=args.histos,
    BinCount=args.bins,
    show=(not args.save),
    verbose=args.verbose,
)

# create 3D view
print(f"reader: PointData={list(reader.PointData.keys())}")
print(f"reader: CellData={list(reader.CellData.keys())}")


# deformed view
suffix = ""
found = False
for field in list(reader.PointData.keys()):
    if field.endswith("displacement"):
        found = True
        break
print(f"displacement found={found} in {list(reader.PointData.keys())}")

if found:
    # make3Dview(cellsize, blockdata, key, color, addruler=True)
    if args.channels:
        # compute channel deformation
        # use MeshLib see test-meshlib example
        for i, block in enumerate(blockdata.keys()):
            name = blockdata[block]["name"]
            actual_name = name.replace("/root/", "")
            if actual_name.startswith("H") and not actual_name.endswith("Isolant"):
                extractBlock1 = ExtractBlock(registrationName=name, Input=reader)
                extractBlock1.Selectors = [block]
                extractBlock1.UpdatePipeline()
                extractSurface1 = ExtractSurface(
                    registrationName="ExtractSurface1", Input=extractBlock1
                )

                SaveData(f"{actual_name}.stl", proxy=extractSurface1)
                Delete(extractBlock1)
                del extractBlock1

    geometry = deformed(cellsize, factor=1, printed=False)
    info(geometry)

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

                SaveData(f"{cwd}/{actual_name}-deformed.stl", proxy=extractSurface1)
                Delete(extractBlock1)
                del extractBlock1
            else:
                print(" ignored", flush=True)
    suffix = "-deformed"
    cellsize = geometry

if args.plots:
    print(f"plots: r={args.r}, z={args.z}")
    # plotOr(reader, r, z, show=(not args.save))# with r=[r1, r2], z: float
    # plotOz(reader, r, z, show=(not args.save)) # with r: float, z=[z1,z2]

if args.views:
    if not args.field:
        for vkey in list(cellsize.PointData.keys()):
            if not vkey in ignored_keys:
                color = ["POINTS", vkey]
                makeAxiview(
                    cellsize,
                    blockdata,
                    vkey,
                    color,
                    suffix="deformed",
                    addruler=False,
                )
        for vkey in list(cellsize.CellData.keys()):
            if not vkey in ignored_keys:
                color = ["CELLS", vkey]
                makeAxiview(
                    cellsize, blockdata, vkey, color, suffix="deformed", addruler=False
                )
    else:
        if args.key in cellsize.PointData.keys():
            color = ["POINTS", args.key]
        elif args.key in cellsize.CellData.keys():
            color = ["CELLS", args.key]
            makeAxiview(
                cellsize, blockdata, args.key, color, suffix="deformed", addruler=False
            )

# for magnetfield:
#   - view contour for magnetic potential (see pv-contours.py)
#   - view glyph for MagneticField
#   - compute spherical expansion coefficients

# # save vtkjs

# # export view
# ExportView('/home/LNCMI-G/trophime/export-scene.vtkjs', view=renderView1, ParaViewGlanceHTML='')
