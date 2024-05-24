import numpy as np
import vtk
import gc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import hist
from math import pi

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
    "Power": {"Symbol": "W", "Units": [ureg.watt, ureg.watt]},
    "Temperature": {"Symbol": "T", "Units": [ureg.degK, ureg.degC]},
    "electric_potential": {"Symbol": "V", "Units": [ureg.volt, ureg.volt]},
    "current_density": {
        "Symbol": "J",
        "Units": [
            ureg.ampere / ureg.meter**2,
            ureg.ampere / ureg.Unit(distance_unit) ** 2,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "current_density_x": {
        "Symbol": "Jx",
        "Units": [
            ureg.ampere / ureg.meter**2,
            ureg.ampere / ureg.Unit(distance_unit) ** 2,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "current_density_y": {
        "Symbol": "Jy",
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
    "Fh": {
        "Symbol": "F",
        "Units": [
            ureg.newton / ureg.meter**3,
            ureg.newton / ureg.Unit(distance_unit) ** 3,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "magnetic_potential": {
        "Symbol": "A",
        "Units": [
            ureg.ampere / ureg.meter,
            ureg.ampere / ureg.Unit(distance_unit),
        ],
    },
    "Bg_MagField": {
        "Symbol": "B_Bg",
        "mSymbol": r"$B_{bg}$",
        "Units": [ureg.tesla, ureg.tesla],
    },
    "magnetic_field": {"Symbol": "B", "Units": [ureg.tesla, ureg.tesla]},
    "displacement": {
        "Symbol": "u",
        "Units": [
            ureg.meter / ureg.second,
            ureg.Unit(distance_unit) / ureg.second,
        ],
        "Exclude": ["Air"],
    },
    "displacement_x": {
        "Symbol": "ux",
        "Units": [
            ureg.meter / ureg.second,
            ureg.Unit(distance_unit) / ureg.second,
        ],
        "Exclude": ["Air"],
    },
    "displacement_y": {
        "Symbol": "uy",
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
    "princial_stress_0": {
        "Symbol": "principal_stress_0",
        "mSymbol": r"$\bar{\bar{\sigma}}_0$",
        "Units": [
            ureg.pascal,
            ureg.megapascal,
        ],
        "Exclude": ["Air"],
    },
    "princial_stress_1": {
        "Symbol": "principal_stress_1",
        "mSymbol": r"$\bar{\bar{\sigma}}_1$",
        "Units": [
            ureg.pascal,
            ureg.megapascal,
        ],
        "Exclude": ["Air"],
    },
    "von_mises_criterions": {
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
            # print(f"DescriptiveStats for datatype={datatype}, key={key}:")
            # print(f"dataset: {len(_dataset[datatype][key])}")

            (physic, fieldname) = key.split(".")
            units = {fieldname: fieldunits[fieldname]["Units"]}

            # Exclude row with Name in fieldunits[fieldname]['Exclude']
            excludeblocks = fieldunits[fieldname]["Exclude"]
            found = False
            for block in excludeblocks:
                if block in name:
                    found = True
                    break

            # for dset in _dataset[datatype][key]:
            #     print(tabulate(dset, headers="keys", tablefmt="psql"))
            if not found:
                df = pd.concat(_dataset[datatype][key])
                print(
                    f"Aggregated DescriptiveStats for datatype={datatype}, key={key}:"
                )

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
                        "Variance",
                    ]
                ]
                # how to: rewrite tab contents using symbol and units
                values = df["Variable"].to_list()
                out_values = [rf"{fieldunits[fieldname]['Symbol']}" for value in values]
                df = df.assign(Variable=out_values)
                for column in ["Minimum", "Mean", "Maximum", "Standard Deviation"]:
                    values = df[column].to_list()
                    out_values = convert_data(units, values, fieldname)
                    df = df.assign(column=[f"{val:.2f}" for val in out_values])

                # watch out:
                # variance is square of (std)
                # M2: moment of order 2 (square of mean)
                # M3: moment of order 3 (cube of mean)
                # M4: moment of order 4 (cube of mean)

                if verbose:
                    print(
                        tabulate(df, headers="keys", tablefmt="psql", showindex=False)
                    )
                df.to_csv(f"{key}-descriptivestats-create.csv")
                dfs.append(df)

    total_df = pd.concat(dfs)
    print(tabulate(total_df, headers="keys", tablefmt="psql", showindex=False))
    total_df.to_csv(f"{name}-descriptivestats-create.csv")

    pass


def createTable(file: str, key: str, name: str):
    import pandas as pd

    csv = pd.read_csv(file)
    keys = csv.columns.values.tolist()
    print(f"createTable: file={file}, key={key}", flush=True)
    # print(tabulate(csv, headers="keys", tablefmt="psql"))

    # drop following keys
    csv.rename(columns={"Block Name": "BlockName"}, inplace=True)
    dropped_keys = ["Row ID", "Cardinality", "Kurtosis", "Skewness", "Sum"]
    csv.drop(columns=dropped_keys, inplace=True)

    # print("createTable: post-process stats table", flush=True)
    # print(tabulate(csv, headers="keys", tablefmt="psql"))

    # for each row "Derived Stats" copy value to "Primary Statistics"
    primary = csv.query("BlockName == 'Primary Statistics'").dropna(axis="columns")
    # print(tabulate(primary, headers="keys", tablefmt="psql"))
    primary.drop(columns=["BlockName"], inplace=True)
    # print("primary:", tabulate(primary, headers="keys", tablefmt="psql"))
    derived = csv.query("BlockName == 'Derived Statistics'").dropna(axis="columns")
    derived.reset_index(drop=True, inplace=True)
    # print("derived:", tabulate(derived, headers="keys", tablefmt="psql"))
    derived.drop(columns=["BlockName"], inplace=True)
    stats_ = primary.join(derived)
    # print(tabulate(stats_, headers="keys", tablefmt="psql"))

    (nrows, ncols) = stats_.shape
    stats_["Name"] = [name for i in range(nrows)]
    # print("join:\n", tabulate(stats_, headers="keys", tablefmt="psql"))
    return stats_


# plot with matplotlib
def plotHistoAxi(file, name: str, key: str, BinCount: int, show: bool = True):

    ax = plt.gca()
    csv = pd.read_csv(file)
    keys = csv.columns.values.tolist()
    # print("histo before scaling")
    # print(tabulate(csv, headers="keys", tablefmt="psql"))

    # get key unit
    print(f"plotHisto: file={file}, key={key}", flush=True)
    (physic, fieldname) = key.split(".")
    symbol = fieldunits[fieldname]["Symbol"]
    msymbol = symbol
    if "mSymbol" in fieldunits[fieldname]:
        msymbol = fieldunits[fieldname]["mSymbol"]
    [in_unit, out_unit] = fieldunits[fieldname]["Units"]
    # print(f"in_units={in_unit}, out_units={out_unit}", flush=True)

    units = {fieldname: fieldunits[fieldname]["Units"]}
    values = csv[key].to_list()
    out_values = convert_data(units, values, fieldname)
    csv[field] = [f"{val:.2f}" for val in out_values]

    hist(csv[key], bins=BinCount, weights=csv["AxiVol"], ax=ax, rot=45)

    plt.xlabel(rf"{msymbol}[{out_unit:~P}]")
    plt.ylabel("Fraction of total Area [%]")
    plt.title(f"{name}: {key}")
    plt.grid(True)
    # plt.legend(False)

    # if legend is mandatory, set legend to True above and comment out the following line
    # ax.legend([rf"{symbol}[{out_unit:~P}]"])
    ax.yaxis.set_major_formatter(lambda x, pos: f"{x:.1f}")
    show = False
    if show:
        plt.show()
    else:
        plt.tight_layout()
        plt.savefig(f"{name}-{key}-histogram-matplotlib.png", dpi=300)
    plt.close()

    # rename columns for tabulate
    csv.rename(
        columns={
            "bin_extents": rf"{symbol} [{out_unit:~P}]",
            "Area_total": "Fraction of total Area [%]",
        },
        inplace=True,
    )
    print(tabulate(csv, headers="keys", tablefmt="psql", showindex=False))

    # check that sum is roughtly equal to 1
    # print(f'check Sum(Fraction): {csv["Fraction of total Area [%]"].sum()}')
    eps = 1.0e-4
    error = abs(1 - csv["Fraction of total Area [%]"].sum() / 100.0)
    print(f"error={error}, check: {(abs(error) <= eps)}, eps={eps}", flush=True)
    assert error <= eps, "Check Sum(Fraction) failed"

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
    print(f"getresultInfo {name}: datadict={datadict}", flush=True)

    return datadict


def resultStats(input, name: str, Area: float, verbose: bool = False):
    """
    compute stats for PointData, CellData and FieldData

    returns a dict
    """

    datadict = resultinfo(input, verbose)
    print(f"resultStats[{name}]: datadict={datadict}", flush=True)

    # PointData to CellData
    # Cellsize
    # cellcenters
    pointDatatoCellData = PointDatatoCellData(
        registrationName="CellDatatoPointData", Input=input
    )

    cellSize1 = CellSize(registrationName="CellSize1", Input=pointDatatoCellData)
    # Properties modified on cellSize1 for 3D
    cellSize1.ComputeVertexCount = 0
    cellSize1.ComputeLength = 0
    cellSize1.ComputeArea = 1  # for 2D
    cellSize1.ComputeSum = 0
    cellSize1.UpdatePipeline()
    # print(f"cellSize1 CellData: {cellSize1.CellData[:]}")

    cellCenters1 = CellCenters(registrationName="CellCenters1", Input=cellSize1)
    # Properties modified on cellCenters1
    cellCenters1.VertexCells = 1

    # add moments for each scalar CellData
    datadict = resultinfo(cellCenters1)
    for datatype in datadict:
        if datatype == "CellData":
            AttributeMode = datadict[datatype]["AttributeMode"]
            TypeMode = datadict[datatype]["TypeMode"]
            for key, kdata in datadict[datatype]["Arrays"].items():
                if not key in ignored_keys:
                    Components = kdata["Components"]
                    bounds = kdata["Bounds"]
                    if bounds[0][0] != bounds[0][1]:
                        if Components == 1:
                            tmp = cellCenters1
                            for order in range(1, 5):
                                cellCenters1 = momentN(
                                    tmp, key, key, order, AttributeMode
                                )

    # compute integrals
    integrateKeys(cellCenters1)

    # reload csv
    csv = pd.read_csv(f"{name}-integrals.csv")

    # post-process csv: divide by Area from CellSize1
    csv = csv / Area

    # drop columns
    # store res in same table as in resultstats (see getresultstats)

    for datatype in datadict:
        if datatype == "PointData":
            AttributeMode = datadict[datatype]["AttributeMode"]
            TypeMode = datadict[datatype]["TypeMode"]
            for key, kdata in datadict[datatype]["Arrays"].items():
                if not key in ignored_keys:
                    Components = kdata["Components"]
                    bounds = kdata["Bounds"]
                    if bounds[0][0] != bounds[0][1]:
                        if not "Stats" in kdata:
                            kdata["Stats"] = {}

                        kdata["Stats"] = # pandas query from csv
                        """
                        print("resultStats:")
                        print(
                            tabulate(
                                kdata["Stats"],
                                headers="keys",
                                tablefmt="psql",
                            )
                        )
                        """

    # display stats
    return datadict


def getresultStats(
    input, name: str, key: str, AttributeMode: str, printed: bool = True
):
    """
    compute stats for key

    in Axi, we cannot use directly DescriptiveStatistics as we need to account for the Jacobian
    of the transformation Cartesian to Cylindrical

    instead we compute integral
    """

    print(
        f"getresultStats: name={name}, key={key}, AttributeMode={AttributeMode}",
        flush=True,
    )

    # if field is a vector, create a field for its norm

    # statistics
    statistics = DescriptiveStatistics(input)
    # get params list
    if not printed:
        for prop in statistics.ListProperties():
            print(f"DescriptiveStatistics: {prop}={statistics.GetPropertyValue(prop)}")
    statistics.VariablesofInterest = [key]
    statistics.AttributeMode = AttributeMode
    # statistics.UpdatePipeline()
    # stats_info = statistics.GetDataInformation()
    # SetActiveSource(statistics)

    spreadSheetView = CreateView("SpreadSheetView")
    descriptiveStatisticsDisplay = Show(
        statistics, spreadSheetView, "SpreadSheetRepresentation"
    )
    spreadSheetView.Update()
    export = ExportView(
        f"{name}-{key}-descriptivestats.csv",
        view=spreadSheetView,
        RealNumberNotation="Scientific",
    )

    # if not printed:
    #     # get params list
    #     for prop in export.ListProperties():
    #         print(f'ExportView: {prop}={export.GetPropertyValue(prop)}')

    # # not working
    export = CreateWriter(f"{name}-{key}-descriptivestats-create.csv", proxy=statistics)
    export.FieldAssociation = "Row Data"
    # get params list
    if not printed:
        for prop in export.ListProperties():
            print(f"CreateWriter: {prop}={export.GetPropertyValue(prop)}")
    export.UpdateVTKObjects()  # is it needed?
    export.UpdatePipeline()

    Delete(descriptiveStatisticsDisplay)
    Delete(statistics)
    del statistics
    Delete(spreadSheetView)
    del spreadSheetView

    csv = createTable(f"{name}-{key}-descriptivestats.csv", key, name)

    return csv


def resultHistos(
    input,
    name: str,
    Area: float,
    BinCount: int = 10,
    printed: bool = True,
):
    """
    histogram
    """
    print(f"resultHistos: name={name}", flush=True)

    pointDatatoCellData = PointDatatoCellData(
        registrationName="CellDatatoPointData", Input=input
    )

    cellSize1 = CellSize(registrationName="CellSize1", Input=pointDatatoCellData)
    # Properties modified on cellSize1 for 3D
    cellSize1.ComputeVertexCount = 0
    cellSize1.ComputeLength = 0
    cellSize1.ComputeArea = 1  # for 2D
    cellSize1.ComputeSum = 0
    cellSize1.UpdatePipeline()
    # print(f"cellSize1 CellData: {cellSize1.CellData[:]}")

    cellCenters1 = CellCenters(registrationName="CellCenters1", Input=cellSize1)
    # Properties modified on cellCenters1
    cellCenters1.VertexCells = 1

    spreadSheetView = CreateView("SpreadSheetView")
    cellCenters1Display = Show(
        cellCenters1, spreadSheetView, "SpreadSheetRepresentation"
    )
    spreadSheetView.Update()
    ExportView(
        f"{name}-Axi-cellcenters-all.csv",
        view=spreadSheetView,
        RealNumberNotation="Scientific",
    )

    csv = pd.read_csv("Axi-B-cellcenters.csv")
    keys = csv.columns.values.tolist()
    print(f"read_csv: csv={name}-Axi-B-cellcenters.csv, keys={keys}", flush=True)

    csv["AxiSurf"] = csv["Points_0"] * csv["Area"] * (2 * pi) / 1.0e-9
    sum = csv["AxiSurf"].sum()
    csv["AxiVol"] = csv["AxiSurf"] / sum
    print(f"Sum(AxiSurf)={sum}")
    print(f'Sum(AxiVol)={csv["AxiVol"].sum()}')

    datadict = resultinfo(cellCenters1)
    for datatype in datadict:
        if datatype == "CellData":
            AttributeMode = datadict[datatype]["AttributeMode"]
            TypeMode = datadict[datatype]["TypeMode"]
            for key, kdata in datadict[datatype]["Arrays"].items():
                if not key in ignored_keys:
                    Components = kdata["Components"]
                    bounds = kdata["Bounds"]
                    if bounds[0][0] != bounds[0][1]:
                        if Components == 1:
                            plotHistoAxi(
                                f"{name}-Axi-B-cellcenters.csv", name, key, BinCount
                            )

    Delete(cellCenters1Display)
    Delete(spreadSheetView)
    Delete(cellCenters1)
    Delete(cellSize1)
    del spreadSheetView
    del cellCenters1
    del cellSize1

    Delete(pointDatatoCellData)
    del pointDatatoCellData

    # Force a garbage collection
    collected = gc.collect()
    print(f"Garbage collector: collected {collected} objects.")


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
    calculator1.Function = f'Points_0*{key}**{order}")'

    calculator1.UpdatePipeline()
    return calculator1


def integrateKeys(input):
    """
    compute integral of Data over area/volume

    to get values use a spreadsheet
    """

    integratevalues = IntegrateVariables(Input=input)

    spreadSheetView = CreateView("SpreadSheetView")
    descriptiveDisplay = Show(
        integratevalues, spreadSheetView, "SpreadSheetRepresentation"
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
    Delete(integratevalues)
    del integratevalues

    # Force a garbage collection
    collected = gc.collect()
    print(f"Garbage collector: collected {collected} objects.")


def meshinfo(
    input, ComputeStats: bool = True, verbose: bool = False, printed: bool = True
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

    # for vector
    print("Add Norm for vectors:")
    calculator = input

    for field in input.PointData:
        if field.GetNumberOfComponents() == 3:
            print(
                f"create {field.Name}norm for {field.Name} PointData vector",
                flush=True,
            )
            calculator = createVectorNorm(
                calculator, field.Name, field.Name, "Point Data"
            )

    for field in input.CellData:
        if field.GetNumberOfComponents() == 2:
            print(f"create norm for {field.Name} CellData vector", flush=True)
            calculator = createVectorNorm(
                calculator, field.Name, field.Name, "Cell Data"
            )
            # cannot create Ur and Ut for CellData vector, apply only to PointData vector

    print("Get mesh size")
    cellsize = CellSize(calculator)  # input

    # set some params
    cellsize.ComputeLength = 0
    cellsize.ComputeArea = 0
    cellsize.ComputeVertexCount = 0
    cellsize.ComputeSum = 1
    # get params list
    if not printed:
        for prop in cellsize.ListProperties():
            print(f"cellsize: {prop}={cellsize.GetPropertyValue(prop)}")

    # apply
    cellsize.UpdatePipeline()
    dataInfo = info(cellsize)

    dataset = sm.Fetch(cellsize)
    np_dataset = dsa.WrapDataObject(dataset)

    blockdata = {}
    # check dataset type
    if dataset.IsA("vtkUnstructuredGrid"):
        print("UnstructuredGrid")
    elif dataset.IsA("vtkMultiBlockDataSet"):
        print("MultiBlockDataSet")
        block = dataset.GetBlock(0)
        Area_data = block.GetCellData().GetArray("Area")
        # print(f'block cellData[Area_data]: {Area_data}, type={type(Area_data)}')
        Area = block.GetFieldData().GetArray("Area")
        # print(f'block fieldData[Area]: {Area}, type={type(Area)}')
        np_Area = np_dataset.FieldData["Area"]
        # print(f'block fieldData[np_Area]: {np_Area}, type={type(np_Area)}, length={algs.shape(np_Area)}')
        print(
            f"block fieldData[np_Area]: total={algs.sum(np_Area)}, parts={algs.shape(np_Area)}"
        )

        hierarchy = dataInfo.GetHierarchy()
        rootnode = hierarchy.GetRootNode()
        rootSelector = f"/{hierarchy.GetRootNodeName()}"
        blocks = hierarchy.GetNumberOfChildren(rootnode)
        print(f"Load blocks: {blocks}")

        for i in range(blocks):
            child = hierarchy.GetChild(rootnode, i)
            name = hierarchy.GetNodeName(child)
            rootChild = f"{rootSelector}/{name}"
            child_info = reader.GetSubsetDataInformation(0, child)  # rootChild)
            # print(f'block[{i}]: {name}, child={child}, childInfo: {child_info}')

            nodes = child_info.GetNumberOfPoints()
            cells = child_info.GetNumberOfCells()
            print(
                f"block[{i}]: {name}, nodes={nodes}, cells={cells}, vol={np_Area.Arrays[i][0]}"
            )

            blockdata[rootChild] = {
                "name": name,
                "nodes": nodes,
                "cells": cells,
                "Area": np_Area.Arrays[i][0],
            }

        # Compute Stats
        stats = []

        print("Data ranges:", flush=True)
        resultinfo(cellsize, verbose)

        print("Data ranges without Air:", flush=True)
        extractBlock1 = ExtractBlock(registrationName="insert", Input=cellsize)
        extractBlock1.Selectors = [
            block for block in blockdata.keys() if not "Air" in block
        ]
        extractBlock1.UpdatePipeline()
        Areas = [blockdata[block]["Area"] for block in extractBlock1.Selectors]
        mergeBlocks1 = MergeBlocks(registrationName="MergeBlocks1", Input=extractBlock1)
        mergeBlocks1.UpdatePipeline()

        statsdict = resultStats(mergeBlocks1, "insert", sum(Areas))
        resultHistos(mergeBlocks1, "insert", sum(Areas), BinCount=20)

        stats.append(statsdict)
        extractBlock1.UpdatePipeline()
        Delete(extractBlock1)
        del extractBlock1

        # aggregate stats data
        createStatsTable([statsdict], "insert")

        if not ComputeStats:
            return cellsize, blockdata, statsdict

        print("Data ranges per block:", flush=True)
        for i, block in enumerate(blockdata.keys()):
            name = blockdata[block]["name"]
            print(f"block[{i}]: extract {block}, name={name}", flush=True)
            extractBlock1 = ExtractBlock(registrationName=name, Input=cellsize)
            extractBlock1.Selectors = [block]
            extractBlock1.UpdatePipeline()
            statsdict = resultStats(extractBlock1, name, blockdata[block]["Area"])
            resultHistos(mergeBlocks1, name, blockdata[block]["Area"], BinCount=20)

            stats.append(statsdict)
            Delete(extractBlock1)
            del extractBlock1

            # for testing purpose only
            # if i == 1:
            #     break

        createStatsTable(stats, name)

    return cellsize, blockdata


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

    warpByVector1.UpdatePipeline()
    return warpByVector1


################################################################
def setCamera(
    renderView,
    Up: tuple = None,
    Angle: float = 30,
    pProjection: bool = True,
    roll: float = 0,
    elevation: float = 0,
    azimuth: float = 0,
):
    """
    adapt camera settings

    ref:
    https://docs.paraview.org/en/latest/ReferenceManual/customizingParaView.html#camera-settings
    https://docs.paraview.org/en/latest/Tutorials/ClassroomTutorials/pythonAndBatchParaViewAndPython.html#control-the-camera
    """

    renderView.ResetCamera()
    camera = renderView.GetActiveCamera()

    if Up is not None:
        # Set
        position = camera.GetPosition()
        focalPoint = camera.GetFocalPoint()
        viewDir = (
            position[0] - focalPoint[0],
            position[1] - focalPoint[1],
            position[2] - focalPoint[2],
        )
        print(
            f"Reset view: viewDir={viewDir}, viewUp={Up}, viewAngle={Angle}", flush=True
        )
        camera.SetViewUp(Up[0], Up[1], Up[2])
        camera.SetViewAngle(Angle)
        camera.SetParallelProjection(pProjection)
        if pProjection:
            camera.SetParallelScale(1)

        print(f"Adjust Camera: Roll={roll}, Elevation={elevation}", flush=True)
        camera.Roll(roll)  # rotate around ??
        # camera.Yaw(45)   # rotate around ??
        # camera.Pitch(45) # rotate around ??
        # camera.Azimuth(45) # rotate around viewUp?
        camera.Elevation(elevation)  # rotate around perpendicular viewUp x Oz?

    else:
        camera.OrthogonalizeViewUp()
        camera.Azimuth(azimuth)

    renderView.Update()
    camera = renderView.GetActiveCamera()
    position = camera.GetPosition()
    focalPoint = camera.GetFocalPoint()
    viewUp = camera.GetViewUp()
    viewAngle = camera.GetViewAngle()
    parallelProjection = camera.GetParallelProjection()
    print(f"position: {position}", flush=True)
    print(f"focalPoint: {focalPoint}", flush=True)
    print(f"viewUp: {viewUp}", flush=True)
    print(f"viewAngle: {viewAngle}", flush=True)
    print(f"parallelProjection: {parallelProjection}", flush=True)
    print(f"Roll: {camera.GetRoll()}", flush=True)

    # print(f"help={dir(camera)}")


def displayField(
    input,
    selectedblocks: list,
    field: str,
    color,
    addruler: bool = True,
    renderView=None,
    filename: str = None,
    viewUp: tuple = None,
    viewAngle: float = 30,
    parallelProjection: bool = False,
    roll: float = 0,
    elevation: float = 0,
    azimuth: float = 0,
    comment: str = None,
    grid: bool = False,
    polargrid: bool = False,
    printed: bool = True,
):
    """
    display field in renderview

    Azimuth - rotate around the vertical axis.
    Elevation - rotate around the horizontal axis in the plane of the screen.
    Roll - rotate around the axis coming out of the screen.
    View angle - basically a zoom in.
    Camera position - where the camera is.
    Focal point - where the camera is looking.
    View Up - I don't know what this is (default = (0, 1, 0) for 3D, view from +Oz)

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
    if polargrid:
        display.PolarAxes.Visibility = 1
        display.PolarAxes.MaximumAngle = 360.0

    # for vector: ColorBy(display, ('CELLS', 'magnetic_field', 'Z'))
    ColorBy(display, tuple(color))

    # Add BoundingRuler filter to get an idea of the dimension
    if addruler:
        print("Add ruler to see dimensions", flush=True)
        ruler = BoundingRuler(registrationName="BoundingRuler1", Input=input)
        ruler.Axis = "Y Axis"
        if not printed:
            for prop in ruler.ListProperties():
                print(f"ruler: {prop}={ruler.GetPropertyValue(prop)}")
        Show(ruler, renderView)  # Reset Camera

    setCamera(
        renderView, viewUp, viewAngle, parallelProjection, roll, elevation, azimuth
    )  # adjustCamera)

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

    (physic, fieldname) = key.split(".")
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
    print(f"make2Dview: field={field}", end="")
    if suffix:
        print(f", suffix={suffix}", end="")
    print(flush=True)
    print(f"blockdata={blockdata}", flush=True)

    (physic, fieldname) = key.split(".")
    print(f"Exclude blocks = {fieldunits[fieldname]['Exclude']}", flush=True)

    selectedblocks = selectBlocks(
        list(blockdata.keys()), fieldunits[fieldname]["Exclude"]
    )
    print(f"input.Selectors = {selectedblocks}", flush=True)

    filename = f"{field}.png"
    if suffix is not None:
        filename = f"{field}-{suffix}.png"

    renderView = displayField(
        input,
        selectedblocks,
        field,
        color,
        addruler=addruler,
        filename=filename,
        viewUp=(0, 1, 0),
        viewAngle=30,
        parallelProjection=False,
        roll=90,
        elevation=300,
    )

    Delete(boxclip)
    del boxclip
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
        (physic, fieldname) = key.split(".")
        print(f"physic={physic}, fieldname={fieldname}", flush=True)
        print(f'fieldunits[fieldname]={fieldunits[fieldname]}"', flush=True)
        symbol = fieldunits[fieldname]["Symbol"]
        [in_unit, out_unit] = fieldunits[fieldname]["Units"]
        print(f"in_units={in_unit}, out_units={out_unit}", flush=True)

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
        (physic, fieldname) = key.split(".")
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
    "--save",
    help="save graphs",
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
cellsize, blockdata, statsdict = meshinfo(reader, ComputeStats=args.stats)

# create 3D view
if not args.field:
    key = None
    for vkey in list(reader.PointData.keys()):
        if not vkey in ignored_keys:
            key = vkey
            break

    field = reader.PointData[key]
    color = ["POINTS", key]
    print(f"force field to {key}", flush=True)

# make3Dview(cellsize, blockdata, key, color, addruler=True)

# deformed view
geometry = deformed(cellsize, factor=1)

# compute channel deformation
# use MeshLib see test-meshlib example
for i, block in enumerate(blockdata.keys()):
    name = blockdata[block]["name"]
    actual_name = name.replace("/root/", "")
    if actual_name.startswith("H") and not actual_name.endswith("Isolant"):
        extractBlock1 = ExtractBlock(registrationName=name, Input=geometry)
        extractBlock1.Selectors = [block]
        extractBlock1.UpdatePipeline()
        extractSurface1 = ExtractSurface(
            registrationName="ExtractSurface1", Input=extractBlock1
        )

        SaveData(f"{actual_name}-deformed.stl", proxy=extractSurface1)
        Delete(extractBlock1)
        del extractBlock1


makeAxiview(geometry, blockdata, key, color, suffix="deformed", addruler=False)

# for magnetfield:
#   - view contour for magnetic potential (see pv-contours.py)
#   - view glyph for MagneticField
#   - compute spherical expansion coefficients

# # save vtkjs

# # export view
# ExportView('/home/LNCMI-G/trophime/export-scene.vtkjs', view=renderView1, ParaViewGlanceHTML='')
