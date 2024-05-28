import numpy as np
import vtk

from tabulate import tabulate

from paraview.simple import *
from paraview.simple import (
    CellSize,
    Clip,
    Slice,
    PlotOnIntersectionCurves,
    Histogram,
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
)

from paraview import servermanager as sm
from paraview.vtk.numpy_interface import dataset_adapter as dsa
from paraview.vtk.numpy_interface import algorithms as algs

import pandas as pd

pd.options.mode.copy_on_write = True

import argparse
import os
import gc

from pint import UnitRegistry, Quantity

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
# !!! watchout matplotlib setup for latex support and packages required !!!

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
    "sigma": {
        "Symbol": "Sigma",
        "mSymobol": r"$\sigma$",
        "Units": [
            ureg.siemens / ureg.meter,
            ureg.siemens / ureg.Unit(distance_unit),
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "nu": {
        "Symbol": "Poisson",
        "Units": [
            ureg.dimensionless,
            ureg.dimensionless,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "EE": {
        "Symbol": "YoungModulus",
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
    "E": {
        "Symbol": "E",
        "Units": [ureg.volt / ureg.meter, ureg.volt / ureg.Unit(distance_unit)],
        "Exclude": ["Air", "Isolant"],
    },
    "E_ur": {
        "Symbol": "Er",
        "Units": [
            ureg.volt / ureg.meter**2,
            ureg.volt / ureg.Unit(distance_unit) ** 2,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "E_ut": {
        "Symbol": "Et",
        "mSymbol": r"$E_{\theta}$",
        "Units": [
            ureg.volt / ureg.meter**2,
            ureg.volt / ureg.Unit(distance_unit) ** 2,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "Enorm": {
        "Symbol": "E",
        "mSymbol": r"$\| E \|$",
        "Units": [
            ureg.volt / ureg.meter**2,
            ureg.volt / ureg.Unit(distance_unit) ** 2,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "J": {
        "Symbol": "J",
        "Units": [
            ureg.ampere / ureg.meter**2,
            ureg.ampere / ureg.Unit(distance_unit) ** 2,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "J_ur": {
        "Symbol": "Jr",
        "Units": [
            ureg.ampere / ureg.meter**2,
            ureg.ampere / ureg.Unit(distance_unit) ** 2,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "J_ut": {
        "Symbol": "Jt",
        "mSymbol": r"$J_{\theta}$",
        "Units": [
            ureg.ampere / ureg.meter**2,
            ureg.ampere / ureg.Unit(distance_unit) ** 2,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "Jnorm": {
        "Symbol": "J",
        "mSymbol": r"$\| J \|$",
        "Units": [
            ureg.ampere / ureg.meter**2,
            ureg.ampere / ureg.Unit(distance_unit) ** 2,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "F_laplace": {
        "Symbol": "F",
        "Units": [
            ureg.newton / ureg.meter**3,
            ureg.newton / ureg.Unit(distance_unit) ** 3,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "F_laplace_ur": {
        "Symbol": "Fr",
        "Units": [
            ureg.newton / ureg.meter**2,
            ureg.newton / ureg.Unit(distance_unit) ** 2,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "F_laplace_ut": {
        "Symbol": "Ft",
        "mSymbol": r"$F_{\theta}$",
        "Units": [
            ureg.newton / ureg.meter**2,
            ureg.newton / ureg.Unit(distance_unit) ** 2,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "F_laplacenorm": {
        "Symbol": "F",
        "mSymbol": r"$\| F \|$",
        "Units": [
            ureg.newton / ureg.meter**2,
            ureg.newton / ureg.Unit(distance_unit) ** 2,
        ],
        "Exclude": ["Air", "Isolant"],
    },
    "Bz": {
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
    "displacement_ur": {
        "Symbol": "ur",
        "Units": [
            ureg.meter / ureg.second,
            ureg.Unit(distance_unit) / ureg.second,
        ],
        "Exclude": ["Air"],
    },
    "displacement_ut": {
        "Symbol": "ut",
        "mSymbol": r"$u_{\theta}$",
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
    "stress_T": {
        "Symbol": "stress_T",
        "mSymbol": r"$\bar{\bar{\sigma}}_{T}$",
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
    "cfpdes.expr.nu",
    "cfpdes.expr.EE",
    "Area",
    "Volume",
    "r",
    "Cos",
    "Sin",
    "cfpdes.pid",
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


def createStatsTable(stats: list, name: str, verbose: bool = False) -> pd.DataFrame:

    # TODO add a column with the Block Name
    # the column Block Name in cvs is not what I think
    # it is either "Primary Statistics" or "Derived Statistics"
    # remove "Derived Statistics" rows
    # why Standard Deviation and Variance are NaN??
    print(
        f"createStatsTable: name={name}, aggregate stats by datatype and field ({len(stats)}) items",
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
                    if not key in _dataset[datatype]:
                        # print(f"create _dataset[{datatype}][{key}]")
                        _dataset[datatype][key] = []

                    # print(f"append kdata[Stats] to _dataset[{datatype}][{key}]")
                    _dataset[datatype][key].append(kdata["Stats"])

    # print("_dataset:", flush=True)
    dfs = []
    for datatype in _dataset:
        # print(f"datatype={datatype}")
        for key in _dataset[datatype]:
            # print(f"DescriptiveStats for datatype={datatype}, key={key}:")
            # print(f"dataset: {len(_dataset[datatype][key])}")

            keyinfo = key.split(".")
            print(f"keyinfo={keyinfo}", flush=True)
            if len(keyinfo) == 2:
                (physic, fieldname) = keyinfo
            elif len(keyinfo) == 3:
                (toolbox, physic, fieldname) = keyinfo
            else:
                raise RuntimeError(f"{key}: cannot get keyinfo as splitted char")
            units = {fieldname: fieldunits[fieldname]["Units"]}
            # print(f"physic={physic}, fieldname={fieldname}", flush=True)
            # print(f'fieldunits[fieldname]={fieldunits[fieldname]}"', flush=True)
            symbol = fieldunits[fieldname]["Symbol"]
            [in_unit, out_unit] = fieldunits[fieldname]["Units"]
            # print(f"in_units={in_unit}, out_units={out_unit}", flush=True)

            """
            # Exclude row with Name in fieldunits[fieldname]['Exclude']
            excludeblocks = fieldunits[fieldname]["Exclude"]
            found = False
            for block in excludeblocks:
                if block in name:
                    found = True
                    print(f"ignore block: {block}", flush=True)
                    break
            """

            # for dset in _dataset[datatype][key]:
            #     print(tabulate(dset, headers="keys", tablefmt="psql"))
            df = pd.concat(_dataset[datatype][key])
            # print(f"Aggregated DescriptiveStats for datatype={datatype}, key={key}:")

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

            # how to: rewrite tab contents using symbol and units
            # see: https://github.com/pandas-dev/pandas/issues/29435
            values = df["Variable"].to_list()
            # print(f"values={values}", flush=True)
            out_values = [value.replace(key, rf"{symbol}") for value in values]
            out_values = [rf"{value} [{out_unit:~P}]" for value in out_values]
            df = df.assign(Variable=out_values)
            # print(f"df[Variable]={df['Variable'].to_list()}", flush=True)
            # print(f"out_values={out_values}", flush=True)
            ndf = {}
            for column in ["Minimum", "Mean", "Maximum", "Standard Deviation"]:
                values = df[column].to_list()
                # print(f"{column}:", flush=True)
                # print(f"values={values}", flush=True)
                out_values = convert_data(units, values, fieldname)
                # print(f"out_values={out_values}", flush=True)
                # print(
                #     f"format out_values={[f'{val:.2f}' for val in out_values]}",
                #     flush=True,
                # )
                ndf[column] = [f"{val:.3f}" for val in out_values]
            # print(f'ndf={ndf}')
            scaled_df = pd.DataFrame.from_dict(ndf)
            for column in ["Minimum", "Mean", "Maximum", "Standard Deviation"]:
                df[column] = scaled_df[column]
                # print(f'df[{column}]={df[column].to_list()}', flush=True)

            # watch out:
            # M2: moment of order 2 (square of mean)
            # M3: moment of order 3 (cube of mean)
            # M4: moment of order 4 (cube of mean)
            # print("change units for Moments")
            ndf = {}
            Munits = [in_unit, out_unit]
            if in_unit == ureg.kelvin:
                Munits = [in_unit, in_unit]
            # print(f"units={units}, type={type(units)}", flush=True)
            # print(f"Munits={Munits}, type={type(Munits)}", flush=True)
            for column in ["M2", "M3", "M4"]:
                # print(f"column={column}", flush=True)
                values = df[column].to_list()
                # print(f"values={values}", flush=True)

                MomentUnits = {column: []}
                for i, unit in enumerate(Munits):
                    # print(f"i={i}", flush=True)
                    unit_ = fieldunits[fieldname]["Units"][i]
                    if fieldunits[fieldname]["Units"][0] == ureg.kelvin:
                        unit_ = ureg.kelvin
                    # print(f"units[{i}]={unit_}", flush=True)
                    # print(f"Munits[{i}]={Munits[i]}", flush=True)
                    MomentUnits[column].append(Munits[i] * unit_)
                    # print(
                    #     f"MomentUnits[{column}]={MomentUnits[column][-1]}",
                    #     flush=True,
                    # )
                    Munits[i] = Munits[i] * unit_
                # print(f"MomentUnits[{column}]={MomentUnits[column]}", flush=True)

                out_values = convert_data(MomentUnits, values, column)
                ndf[column] = [f"{val:.3f}" for val in out_values]
                del MomentUnits[column]

            scaled_df = pd.DataFrame.from_dict(ndf)
            for column in ["M2", "M3", "M4"]:
                # print(f"df[{column}]={df[column].to_list()}", flush=True)
                df[column] = scaled_df[column]
                # print(f"scaled df[{column}]={df[column].to_list()}", flush=True)

            if verbose:
                print(
                    tabulate(df, headers="keys", tablefmt="psql", showindex=False),
                    flush=True,
                )
            df.to_csv(f"{key}-descriptivestats.csv")
            dfs.append(df)

    total_df = pd.DataFrame()
    if dfs:
        total_df = pd.concat(dfs)
        print(
            tabulate(total_df, headers="keys", tablefmt="psql", showindex=False),
            flush=True,
        )
        total_df.to_csv(f"{name}-descriptivestats.csv")

        # remove temporary csv files
        for datatype in _dataset:
            for key in _dataset[datatype]:
                os.remove(f"{key}-descriptivestats.csv")
                os.remove(f"{name}-{key}-descriptivestats.csv")

    return total_df


def createTable(file: str, key: str, name: str, verbose: bool = False):

    csv = pd.read_csv(file)
    keys = csv.columns.values.tolist()
    if verbose:
        print(f"createTable: file={file}, key={key}", flush=True)
    # print(tabulate(csv, headers="keys", tablefmt="psql"))

    # drop following keys
    csv.rename(columns={"Block Name": "BlockName"}, inplace=True)
    dropped_keys = ["Row ID", "Cardinality", "Kurtosis", "Skewness", "Sum", "Variance"]
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
def plotHisto(
    file, name: str, key: str, Area: float, show: bool = True, verbose: bool = False
):
    import matplotlib.pyplot as plt

    ax = plt.gca()
    csv = pd.read_csv(file)
    csv = csv[["bin_extents", "Area_total"]]
    keys = csv.columns.values.tolist()
    # print("histo before scaling")
    # print(tabulate(csv, headers="keys", tablefmt="psql"))

    # get key unit
    if verbose:
        print(f"plotHisto: file={file}, key={key}", flush=True)
    keyinfo = key.split(".")
    # print(f"keyinfo={keyinfo}", flush=True)
    if len(keyinfo) == 2:
        (physic, fieldname) = keyinfo
    elif len(keyinfo) == 3:
        (toolbox, physic, fieldname) = keyinfo
    symbol = fieldunits[fieldname]["Symbol"]
    msymbol = symbol
    if "mSymbol" in fieldunits[fieldname]:
        msymbol = fieldunits[fieldname]["mSymbol"]
    [in_unit, out_unit] = fieldunits[fieldname]["Units"]
    # print(f"in_units={in_unit}, out_units={out_unit}", flush=True)

    units = {fieldname: fieldunits[fieldname]["Units"]}
    values = csv["bin_extents"].to_list()
    out_values = convert_data(units, values, fieldname)
    csv = csv.assign(bin_extents=out_values)  # [f"{val:.2f}" for val in out_values])
    csv["Area_total"] = csv["Area_total"] / Area * 100

    csv.plot.bar(
        x="bin_extents",
        y="Area_total",
        xlabel=rf"{msymbol}[{out_unit:~P}]",
        ylabel="Fraction of total Area [%]",
        title=f"{name}: {key}",
        grid=True,
        legend=False,
        rot=45,
        ax=ax,
    )

    # if legend is mandatory, set legend to True above and comment out the following line
    # ax.legend([rf"{symbol}[{out_unit:~P}]"])
    ax.yaxis.set_major_formatter(lambda x, pos: f"{x:.1f}")
    ax.xaxis.set_major_formatter(lambda x, pos: f"{x:.3f}")
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
    if verbose:
        print(
            tabulate(csv, headers="keys", tablefmt="psql", showindex=False), flush=True
        )

    # check that sum is roughtly equal to 1
    # print(f'check Sum(Fraction): {csv["Fraction of total Area [%]"].sum()}')
    eps = 1.0e-4
    error = abs(1 - csv["Fraction of total Area [%]"].sum() / 100.0)
    if error > eps:
        print(
            f"histogram name={name}, key={key}: Check Sum(Fraction) failed : error={error}, eps={eps}"
        )
    # assert error <= eps, f"Check Sum(Fraction) failed : error={error}, eps={eps}"

    csv.to_csv(f"{name}-{key}-histogram-matplotlib.csv")
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


def resultinfo(input, verbose: bool = False) -> dict:
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


def getresultInfo(key, verbose: bool = False, printed: bool = True) -> dict:
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


def resultStats(
    input,
    name: str,
    Area: float,
    histo: bool = False,
    BinCount: int = 10,
    printed: bool = True,
    show: bool = False,
    verbose: bool = False,
) -> dict:
    """
    compute stats for PointData, CellData and FieldData

    returns a dict
    """

    datadict = resultinfo(input, verbose)
    if verbose:
        print(f"resultStats[{name}]: datadict={datadict}", flush=True)
    for datatype in datadict:
        if datatype != "FieldData":
            AttributeMode = datadict[datatype]["AttributeMode"]
            TypeMode = datadict[datatype]["TypeMode"]
            for key, kdata in datadict[datatype]["Arrays"].items():
                if not key in ignored_keys:

                    found = False
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

                    for excluded in fieldunits[fieldname]["Exclude"]:
                        if excluded in name:
                            found = True
                            if verbose:
                                print(f"ignore block: {name}", flush=True)
                            break

                    if not found:
                        Components = kdata["Components"]
                        bounds = kdata["Bounds"]
                        if abs(bounds[0][0] - bounds[0][1]) <= 1.0e-3:
                            if not "Stats" in kdata:
                                # print(f"\t{key}: create kdata[Stats]")
                                kdata["Stats"] = {}

                            kdata["Stats"] = getresultStats(
                                input, name, key, AttributeMode, verbose=verbose
                            )

                        if histo:
                            getresultHisto(
                                input,
                                name,
                                Area,
                                key,
                                TypeMode,
                                Components,
                                BinCount=BinCount,
                                show=show,
                                verbose=verbose,
                            )

    # display stats
    return datadict


def getresultStats(
    input,
    name: str,
    key: str,
    AttributeMode: str,
    printed: bool = True,
    verbose: bool = False,
):
    """
    compute stats for key
    """

    if verbose:
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

    filename = f"{name}-{key}-descriptivestats.csv"
    export = ExportView(
        filename,
        view=spreadSheetView,
        RealNumberNotation="Scientific",
    )

    # if not printed:
    #     # get params list
    #     for prop in export.ListProperties():
    #         print(f'ExportView: {prop}={export.GetPropertyValue(prop)}')

    Delete(descriptiveStatisticsDisplay)
    Delete(spreadSheetView)
    del spreadSheetView
    Delete(statistics)
    del statistics

    csv = createTable(filename, key, name)

    return csv


def getresultHisto(
    input,
    name: str,
    Area: float,
    key: str,
    TypeMode: str,
    Components: int = 1,
    BinCount: int = 10,
    printed: bool = True,
    verbose: bool = False,
):
    """
    histogram
    """
    print(f"getresultHisto: name={name}, key={key}, TypeMode={TypeMode}", flush=True)

    # convert pointdata to celldata
    if TypeMode == "POINT":
        pointDatatoCellData = PointDatatoCellData(
            registrationName="pointDatatoCellData", Input=input
        )

    elif TypeMode == "CELL":
        pointDatatoCellData = input

    else:
        print(
            'resultHisto: not applicable for {key["Name"]} - unsupported data type {key["Type"]}'
        )
        return

    cellSize1 = CellSize(registrationName="CellSize1", Input=pointDatatoCellData)
    # Properties modified on cellSize1 for 3D
    cellSize1.ComputeVertexCount = 0
    cellSize1.ComputeLength = 0
    cellSize1.ComputeArea = 1  # for 3D
    cellSize1.ComputeSum = 0
    cellSize1.UpdatePipeline()
    # print(f"cellSize1 CellData: {cellSize1.CellData[:]}")

    histogram1 = Histogram(registrationName="Histogram1", Input=cellSize1)
    histogram1.SelectInputArray = ["CELLS", key]
    # for scalar comment out Component
    if Components > 1:
        histogram1.Component = "Magnitude"
    histogram1.CalculateAverages = 1
    histogram1.CenterBinsAroundMinAndMax = 1
    histogram1.BinCount = BinCount
    if not printed:
        # get params list
        for prop in histogram1.ListProperties():
            print(f"Histogram: {prop}={histogram1.GetPropertyValue(prop)}")
    # TODO from key range
    # histogram1.CustomBinRanges = [min, max]

    # Properties modified on histogram1
    histogram1.CalculateAverages = 1

    filename = f"{name}-{key}-histogram.csv"
    export = CreateWriter(filename, proxy=histogram1)
    if not printed:
        for prop in export.ListProperties():
            print(f"export: {prop}={export.GetPropertyValue(prop)}")
    export.UpdateVTKObjects()  # is it needed?
    export.UpdatePipeline()
    Delete(histogram1)
    del histogram1

    plotHisto(filename, name, key, Area)

    # remove temporary csv files
    os.remove(filename)

    if TypeMode == "POINT":
        Delete(pointDatatoCellData)
        del pointDatatoCellData


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


def meshinfo(
    input,
    ComputeStats: bool = True,
    ComputeHisto: bool = False,
    BinCount: int = 10,
    show: bool = False,
    verbose: bool = False,
    printed: bool = True,
) -> tuple:
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

    def rectTocylField(input, key: str, nkey: str, AttributeType: str):
        """
        compute r and theta component of a vector Field

        vr = ux * cos + uy * sin
        vt = ux* -sin + uy * cos

        with cos = x/sqrt(x²+y²) and sin = y/sqrt(x²+y²)

        x: coordsX
        y : coordsY

        ux: "thermo_electric.electric.current_density_X"

        must be: AttributeType: 'Point Data'
        """

        if AttributeType != "Point Data":
            raise RuntimeError(
                f"rectTocylField: {key} - unsupported AttributeType: {AttributeType}"
            )
        inputDataPoint = [field.Name for field in input.PointData]
        print(
            f"rectTocylField: input PointData = {inputDataPoint}",
            flush=True,
        )

        # check if r exists already
        if "r" in inputDataPoint:
            # skip the next 3 steps: aka calculator1 to calculator3
            calculator3 = input

        else:
            calculator1 = Calculator(registrationName="Calculator1", Input=input)
            calculator1.AttributeType = AttributeType  # 'Cell Data'
            calculator1.ResultArrayName = "r"
            calculator1.Function = "sqrt(coordsX*coordsX+coordsY*coordsY)"
            calculator1.UpdatePipeline()

            calculator2 = Calculator(Input=calculator1)
            calculator2.AttributeType = AttributeType  # 'Cell Data'
            calculator2.ResultArrayName = "Cos"
            calculator2.Function = "coordsX/r"
            calculator2.UpdatePipeline()

            calculator3 = Calculator(Input=calculator2)
            calculator3.AttributeType = AttributeType  # 'Cell Data'
            calculator3.ResultArrayName = "Sin"
            calculator3.Function = "coordsY/r"
            calculator3.UpdatePipeline()

        calculator4 = Calculator(Input=calculator3)
        calculator4.AttributeType = AttributeType  # 'Cell Data'
        calculator4.ResultArrayName = f"{key}_ur"
        calculator4.Function = f'"{key}_X"*Cos+"{key}_Y"*Sin'
        calculator4.UpdatePipeline()

        calculator5 = Calculator(Input=calculator4)
        calculator5.AttributeType = AttributeType  # 'Cell Data'
        calculator5.ResultArrayName = f"{key}_ut"
        calculator5.Function = f'-"{key}_X"*Sin+"{key}_Y"*Cos'

        calculator5.UpdatePipeline()
        return calculator5

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

    # rectTocyl: need CellDataToPointData before
    # for temperature add, for forces and densities norm, rescale
    cellDatatoPointData1 = CellDatatoPointData(
        registrationName="CellDatatoPointData", Input=input
    )

    # for vector
    print("Add Norm for vectors and RectToCyl:")
    calculator = cellDatatoPointData1

    for field in cellDatatoPointData1.PointData:
        if field.GetNumberOfComponents() > 1:
            print(
                f"create {field.Name}norm for {field.Name} PointData vector",
                flush=True,
            )
            calculator = createVectorNorm(
                calculator, field.Name, field.Name, "Point Data"
            )
            print(
                f"create {field.Name}ur and {field.Name}ut for {field.Name} PointData vector",
                flush=True,
            )
            calculator = rectTocylField(
                calculator, field.Name, field.Name, "Point Data"
            )

    print("Get mesh size")
    cellsize = CellSize(calculator)  # input

    # set some params
    cellsize.ComputeLength = 0
    cellsize.ComputeArea = 1
    cellsize.ComputeVolume = 1
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
        tvol = algs.sum(np_Area)
        vunits = fieldunits["Area"]["Units"]
        mm2 = f"{vunits[1]:~P}"
        tvol_mm2 = convert_data(
            {"Area": vunits},
            tvol,
            "Area",
        )
        print(
            f"block fieldData[np_Area]: total={tvol_mm2} {mm2}, parts={algs.shape(np_Area)}"
        )

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
            vol = np_Area.Arrays[i][0]
            vol_mm2 = convert_data(
                {"Area": vunits},
                vol,
                "Area",
            )
            print(
                f"block[{i}]: {name}, nodes={nodes}, cells={cells}, vol={vol_mm2} {mm2}"
            )

            blockdata[rootChild] = {
                "name": name,
                "nodes": nodes,
                "cells": cells,
                "Area": vol,
            }

            sum_vol += vol

        # check tvol == Sum(vol)
        if abs(1 - sum_vol / tvol) > 1.0e-3:
            raise RuntimeError(f"Total area != Sum(area), error={abs(1-sum_vol/tvol)}")

        # Compute Stats
        stats = []

        print("Data ranges:", flush=True)
        datadict = resultinfo(cellsize)

        print("Data ranges without Air:", flush=True)
        extractBlock1 = ExtractBlock(registrationName="insert", Input=cellsize)
        extractBlock1.Selectors = [
            block for block in blockdata.keys() if not "Air" in block
        ]
        extractBlock1.UpdatePipeline()
        Areas = [blockdata[block]["Area"] for block in extractBlock1.Selectors]
        mergeBlocks1 = MergeBlocks(registrationName="MergeBlocks1", Input=extractBlock1)
        mergeBlocks1.UpdatePipeline()
        statsdict = resultStats(
            mergeBlocks1,
            "insert",
            sum(Areas),
            histo=ComputeHisto,
            BinCount=BinCount,
            show=show,
            verbose=verbose,
        )
        if verbose:
            print(f"insert statsdict={statsdict}")
        stats.append(statsdict)
        extractBlock1.UpdatePipeline()
        Delete(extractBlock1)
        del extractBlock1

        # Force a garbage collection
        collected = gc.collect()
        if verbose:
            print(f"Garbage collector: collected {collected} objects.")

        # aggregate stats data
        createStatsTable([statsdict], "insert")

        if not ComputeStats:
            return cellsize, blockdata, statsdict

        if len(blockdata.keys()) > 1:
            print("Data ranges per block:", flush=True)
            for i, block in enumerate(blockdata.keys()):
                name = blockdata[block]["name"]
                print(f"block[{i}]: extract {block}, name={name}", flush=True)
                extractBlock1 = ExtractBlock(registrationName=name, Input=cellsize)
                extractBlock1.Selectors = [block]
                extractBlock1.UpdatePipeline()
                statsdict = resultStats(
                    extractBlock1,
                    name,
                    blockdata[block]["Area"],
                    histo=ComputeHisto,
                    BinCount=BinCount,
                    show=show,
                    verbose=verbose,
                )
                stats.append(statsdict)
                Delete(extractBlock1)
                del extractBlock1

                # Force a garbage collection
                collected = gc.collect()
                if verbose:
                    print(f"Garbage collector: collected {collected} objects.")

                # aggregate stats data
                createStatsTable([statsdict], name)

            # aggregate stats data
            createStatsTable(stats, "total")

    return cellsize, blockdata, dict()


def deformed(input, factor: float = 1, printed: bool = True):
    """
    create deformed geometry
    """
    print(f"deformed view: factor={factor}", flush=True)

    print(f"input: PointData={list(input.PointData.keys())}")
    print(f"input: CellData={list(input.CellData.keys())}")

    warpByVector1 = WarpByVector(registrationName="WarpByVector1", Input=input)
    warpByVector1.Vectors = ["POINTS", "elasticity.displacement"]
    warpByVector1.ScaleFactor = factor
    # get params list
    if not printed:
        for prop in warpByVector1.ListProperties():
            print(f"warpByVector1: {prop}={warpByVector1.GetPropertyValue(prop)}")

    warpByVector1.UpdatePipeline()
    print(f"warpByVector1: PointData={list(warpByVector1.PointData.keys())}")
    print(f"warpByVector1: CellData={list(warpByVector1.CellData.keys())}")

    UpdatePipeline()
    return warpByVector1


def makeclip(input, name: str, invert: bool = True, printed: bool = True):
    """
    create a plane clip from input dataset
    """
    print(f"makeclip: name={name}, invert={invert}", flush=True)

    clip = Clip(registrationName=name, Input=input)
    clip.ClipType = "Plane"
    clip.HyperTreeGridClipper = "Plane"

    # init the 'Plane' selected for 'ClipType'
    clip.ClipType.Origin = [0.0, 0.0, 0.0]
    clip.ClipType.Normal = [0.0, 1.0, 0.0]
    clip.Invert = 1
    if invert:
        clip.Invert = 0

    # get params list
    if not printed:
        for prop in clip.ListProperties():
            print(f"clip: {prop}={clip.GetPropertyValue(prop)}")

    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip.HyperTreeGridClipper.Origin = [0.0, 0.0, 0.0]

    clip.UpdatePipeline()
    return clip


def makecylinderslice(input, name: str, r: float, printed: bool = True):
    """
    create a cylinder slice from input dataset
    """
    print(f"makecylinderslice: name={name}, r={r}", flush=True)

    # create a new 'Slice'
    slice1 = Slice(registrationName=name, Input=input)

    # Properties modified on slice1
    slice1.SliceType = "Cylinder"

    # Properties modified on slice1.SliceType
    slice1.SliceType.Center = [0.0, 0.0, 0.0]
    slice1.SliceType.Axis = [0.0, 0.0, 1.0]
    slice1.SliceType.Radius = r

    # get params list
    if not printed:
        for prop in slice1.ListProperties():
            print(f"slice: {prop}={slice1.GetPropertyValue(prop)}")

    slice1.UpdatePipeline()
    return slice1


def makeplaneOrOzslice(input, name: str, theta: float = 0.0):
    """
    create a OrOz plane slice from input dataset

    input:
    name:
    theta: angle of normal in degrees
    """

    from math import cos, sin, pi

    print(f"makeplaneOrOzslice: name={name}, theta={theta}", flush=True)
    radian = theta * pi / 180.0 + pi / 2.0

    # create a new 'Slice'
    slice1 = Slice(registrationName=name, Input=input)

    # Properties modified on slice1
    slice1.SliceType = "Plane"

    # Properties modified on slice1.SliceType
    slice1.SliceType.Origin = [0.0, 0.0, 0.0]
    slice1.SliceType.Normal = [cos(radian), sin(radian), 0.0]

    slice1.UpdatePipeline()
    return slice1


def makeplaneslice(input, name: str, z: float = 0.0):
    """
    create a plane slice from input dataset
    """
    print(f"makeplaneslice: name={name}, z={z}", flush=True)

    # create a new 'Slice'
    slice1 = Slice(registrationName=name, Input=input)

    # Properties modified on slice1
    slice1.SliceType = "Plane"

    # Properties modified on slice1.SliceType
    slice1.SliceType.Origin = [0.0, 0.0, z]
    slice1.SliceType.Normal = [0.0, 0.0, 1.0]

    slice1.UpdatePipeline()
    return slice1


def makesphereslice(
    input,
    name: str,
    radius: float,
    r: float = 0,
    theta: float = 0,
    z: float = 0,
):
    """
    create a plane slice from input dataset
    """
    from math import pi, cos, sin

    radian = theta * pi / 180.0

    print(
        f"makesphereslice: name={name}, radius={radius}, r={r}, theta={theta}, z={z}",
        flush=True,
    )

    # create a new 'Slice'
    slice1 = Slice(registrationName=name, Input=input)

    # Properties modified on slice1
    slice1.SliceType = "Sphere"

    # Properties modified on slice1.SliceType
    slice1.SliceType.Center = [r * cos(radian), r * sin(radian), z]
    slice1.SliceType.Radius = radius

    slice1.UpdatePipeline()
    return slice1


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
    polargrid: bool = False,
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
    if polargrid:
        display.PolarAxes.Visibility = 1
        # display.PolarAxes.MaximumAngle = 360.0

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
    print(f"keyinfo={keyinfo}", flush=True)
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
def make2Dview(
    input,
    blockdata,
    field: str,
    color,
    suffix: str = None,
    addruler: bool = False,
    printed: bool = True,
):
    """
    create a 3D view
    """
    print(f"make2Dview: field={field}", end="")
    if suffix:
        print(f", suffix={suffix}", end="")
    print(flush=True)
    print(f"blockdata={blockdata}", flush=True)

    keyinfo = field.split(".")
    print(f"keyinfo={keyinfo}", flush=True)
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
        input,
        selectedblocks,
        field,
        color,
        addruler=addruler,
        filename=filename,
    )

    Delete(renderView)
    del renderView


def plotOr(
    input,
    r: list[float],
    theta: float,
    show: bool = True,
    printed: bool = True,
):
    from math import pi, cos, sin

    [r0, r1] = r
    radian = theta * pi / 180.0

    plotOverLine = PlotOverLine(registrationName="Oz", Input=input, Source="Line")

    # init the 'Line' selected for 'Source'
    plotOverLine.Source.Point1 = [r0 * cos(radian), r0 * sin(radian), 0]
    plotOverLine.Source.Point2 = [r1 * cos(radian), r1 * sin(radian), 0]

    filename = f"r0={r0}m-r1={r1}m-theta={theta}deg.csv"
    export = CreateWriter(filename, proxy=plotOverLine)
    if not printed:
        for prop in export.ListProperties():
            print(f"export: {prop}={export.GetPropertyValue(prop)}")
    export.UpdateVTKObjects()  # is it needed?
    export.UpdatePipeline()

    # plot with matplotlib
    def plotOrField(file, key: str, theta: float, ax=None, show: bool = True):
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator

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
            plt.savefig(f"{key}-vs-r-theta={theta}deg.png", dpi=300)
        plt.close()
        keycsv.to_csv(f"{key}-vs-r-theta={theta}deg.csv")
        pass

    # requirements: create PointData from CellData
    for field in input.PointData:
        plotOrField(
            filename,
            field,
            theta,
            ax=None,
            show=show,
        )

    # remove temporary csv files
    os.remove(filename)

    Delete(plotOverLine)
    del plotOverLine
    pass


def plotTheta(
    input, r: float, show: bool = True, printed: bool = True, verbose: bool = False
):
    """
    for theta, need to apply CellDataToPointData filter
    """

    print(f"plotTheta: r={r}")
    cellDatatoPointData1 = CellDatatoPointData(
        registrationName="CellDatatoPointData", Input=input
    )

    # create clip with plane (howto give a color for each clip)
    print("cellDatatoPointDatalip up and down")
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
                    f"plotOnIntersectionCurve: {prop}={plotOnIntersectionCurve.GetPropertyValue(prop)}"
                )

        plotOnIntersectionCurve.SliceType = "Plane"
        # Properties modified on plotOnIntersectionCurves1.SliceType
        plotOnIntersectionCurve.SliceType.Origin = [0.0, 0.0, 0]
        plotOnIntersectionCurve.SliceType.Normal = [0.0, 0.0, 1.0]

        # plotOnIntersectionCurves.append(plotOnIntersectionCurve)
        export = CreateWriter(f"r={r}m-{i}.csv", proxy=plotOnIntersectionCurve)
        if not printed:
            for prop in export.ListProperties():
                print(f"export: {prop}={export.GetPropertyValue(prop)}")
        export.UpdateVTKObjects()  # is it needed?
        export.UpdatePipeline()
        files.append(f"r={r}m-{i}.csv")

    # plot with matplotlib
    def plotThetaField(files, key: str, r: float, ax=None, show: bool = True):
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
        from numpy import pi, arctan2

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
            # print(f"new keys: {keycsv.columns.values.tolist()}")

            # rescale columns to plot
            theta = arctan2(keycsv["y"].to_numpy(), keycsv["x"].to_numpy()) * 180 / pi
            # print(f'theta (type={type(theta)}, nelem={theta.shape}): {theta}')
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
        df.to_csv(f"{key}-vs-theta-r={r_mm}{mm}.csv")
        # print(f"df keys: {df.columns.values.tolist()}")
        # assert key in df.columns.values.tolist(), f"{key} not in df_keys"
        # print(f"df={df}")
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
            plt.savefig(f"{key}-vs-theta-r={r_mm}{mm}.png", dpi=300)
        plt.close()

        print(f"{key} stats on r={r_mm}{mm}")
        print(f"{df[key].describe()}")
        pass

    # requirements: create PointData from CellData
    datadict = resultinfo(cellDatatoPointData1)
    for field in datadict["PointData"]["Arrays"]:
        if not field in ignored_keys:
            kdata = datadict["PointData"]["Arrays"][field]
            Components = kdata["Components"]
            print(f"plotThetaField for {field} - components={Components}")
            if Components > 1:
                for i in range(Components):
                    print(f"plotThetaField for {field}:{i} skipped")
                # for i in range(Components):
                #     plotThetaField(files, f"{field}:{i}", r, z, ax=None, show=show)
            else:
                plotThetaField(files, field, r, ax=None, show=show)

    # remove temporary csv files
    for file in files:
        os.remove(file)

    Delete(cellDatatoPointData1)
    del cellDatatoPointData1

    # Force a garbage collection
    collected = gc.collect()
    if verbose:
        print(f"Garbage collector: collected {collected} objects.")


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
parser.add_argument("--r", nargs="*", type=float, help="select r in m to display")
parser.add_argument(
    "--theta", nargs="*", type=float, help="select theta in deg to display"
)
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

# get current working directory
cwd = os.getcwd()

# check paraview version
version = GetParaViewVersion()
print(f"Paraview version: {version}")

# args.file = "../../HL-31/test/hybride-Bh27.7T-Bb9.15T-Bs9.05T_HPfixed_BPfree/bmap/np_32/thermo-electric.exports/Export.case"
reader = load(args.file)
# print(f"help(reader) = {dir(reader)}")
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
        # if field.GetNumberOfComponents() == 2:
        #    color = ["CELLS", args.field, "Magnitude"]
    if args.field in list(reader.PointData.keys()):
        field = reader.PointData[args.field]
        # if field.GetNumberOfComponents() == 1:
        color = ["POINTS", args.field]
        # if field.GetNumberOfComponents() == 2:
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


# Plots
if args.plots:
    if args.r:
        for r in args.r:
            plotTheta(cellsize, r, show=(not args.save))
            for theta in args.theta:
                plotOr(cellsize, r, theta, show=(not args.save))
    # add plotOr
    # add plotOz

# When dealing with elasticity
suffix = ""
geometry = None
datadict = resultinfo(cellsize)
found = False
for field in list(reader.PointData.keys()):
    if field.endswith("displacement"):
        found = True
        break
print(f"displacement found={found} in {list(reader.PointData.keys())}")

if found:
    geometry = deformed(cellsize, factor=1)

    suffix = "-deformed"
    cellsize = geometry

# Views
if args.views:
    if not args.field:
        for vkey in list(cellsize.PointData.keys()):
            if not vkey in ignored_keys:
                color = ["POINTS", vkey]
                make2Dview(
                    cellsize, blockdata, vkey, color, suffix="deformed", addruler=False
                )
        for vkey in list(cellsize.CellData.keys()):
            if not vkey in ignored_keys:
                color = ["CELLS", vkey]
                make2Dview(
                    cellsize, blockdata, vkey, color, suffix="deformed", addruler=False
                )
    else:
        if args.key in cellsize.PointData.keys():
            color = ["POINTS", args.key]
        elif args.key in cellsize.CellData.keys():
            color = ["CELLS", args.key]
            make2Dview(
                cellsize,
                blockdata,
                args.key,
                color,
                suffix="deformed",
                addruler=False,
            )

# for magnetfield:
#   - view contour for magnetic potential (see pv-contours.py)
#   - view glyph for MagneticField
#   - compute spherical expansion coefficients

# # save vtkjs

# # export view
# ExportView('/home/LNCMI-G/trophime/export-scene.vtkjs', view=renderView1, ParaViewGlanceHTML='')
