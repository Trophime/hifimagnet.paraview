import pandas as pd
import os

from tabulate import tabulate
from typing import List
from math import pi, sqrt

from paraview.simple import (
    DescriptiveStatistics,
    Delete,
    ExportView,
    CreateView,
    Show,
)

from .method import convert_data, resultinfo
from .histo import getresultHisto


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


def createStatsTable(
    stats: list, name: str, fieldunits: dict, basedir: str, ureg, verbose: bool = False
) -> pd.DataFrame:

    os.makedirs(f"{basedir}/stats", exist_ok=True)
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
                        # print(f"create _dataset[{datatype}][{key}]", flush=True)
                        _dataset[datatype][key] = []

                    # print(f"append kdata[Stats] to _dataset[{datatype}][{key}]", flush=True)
                    _dataset[datatype][key].append(kdata["Stats"])

    # print("_dataset:", flush=True)
    dfs = []
    for datatype in _dataset:
        # print(f"datatype={datatype}", flush=True)
        for key in _dataset[datatype]:
            # print(f"DescriptiveStats for datatype={datatype}, key={key}:", flush=True)
            # print(f"dataset: {len(_dataset[datatype][key])}", flush=True)

            toolbox = ""
            keyinfo = key.split(".")
            print(f"keyinfo={keyinfo}", flush=True)
            if len(keyinfo) == 2:
                (physic, fieldname) = keyinfo
            elif len(keyinfo) == 3:
                (toolbox, physic, fieldname) = keyinfo
            else:
                raise RuntimeError(f"{key}: cannot get keyinfo as splitted char")
            units = {fieldname: fieldunits[fieldname]["Units"]}
            # print(f"toolxbox={toolbox}, physic={physic}, fieldname={fieldname}", flush=True)
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
            # print(f"Aggregated DescriptiveStats for datatype={datatype}, key={key}:", flush=True)

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
            # print(f'ndf={ndf}', flush=True)
            scaled_df = pd.DataFrame.from_dict(ndf)
            for column in ["Minimum", "Mean", "Maximum", "Standard Deviation"]:
                df[column] = scaled_df[column]
                # print(f'df[{column}]={df[column].to_list()}', flush=True)

            # watch out:
            # M2: moment of order 2 (square of mean)
            # M3: moment of order 3 (cube of mean)
            # M4: moment of order 4 (cube of mean)
            # print("change units for Moments", flush=True)
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
            df.to_csv(f"{basedir}/stats/{key}-descriptivestats.csv")
            dfs.append(df)

    total_df = pd.DataFrame()
    if dfs:
        total_df = pd.concat(dfs)
        print(
            tabulate(total_df, headers="keys", tablefmt="psql", showindex=False),
            flush=True,
        )
        total_df.to_csv(f"{basedir}/stats/{name}-descriptivestats.csv")

        # remove temporary csv files
        for datatype in _dataset:
            for key in _dataset[datatype]:
                try:
                    os.remove(f"{basedir}/stats/{key}-descriptivestats.csv")
                    os.remove(f"{basedir}/stats/{name}-{key}-descriptivestats.csv")
                except:
                    pass

    return total_df


def getresultStats(
    input,
    name: str,
    key: str,
    AttributeMode: str,
    basedir: str,
    printed: bool = True,
    verbose: bool = False,
):
    """
    compute stats for key
    """
    os.makedirs(f"{basedir}/stats", exist_ok=True)

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
            print(
                f"DescriptiveStatistics: {prop}={statistics.GetPropertyValue(prop)}",
                flush=True,
            )
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

    filename = f"{basedir}/stats/{name}-{key}-descriptivestats.csv"
    export = ExportView(
        filename,
        view=spreadSheetView,
        RealNumberNotation="Scientific",
    )

    # if not printed:
    #     # get params list
    #     for prop in export.ListProperties():
    #         print(f'ExportView: {prop}={export.GetPropertyValue(prop)}', flush=True)

    Delete(descriptiveStatisticsDisplay)
    Delete(spreadSheetView)
    del spreadSheetView
    Delete(statistics)
    del statistics

    csv = createTable(filename, key, name, verbose)

    return csv


def resultStats(
    input,
    name: str,
    dim: int,
    AreaorVolume: float,
    fieldunits: dict,
    ignored_keys: List[str],
    ureg,
    basedir: str,
    histo: bool = False,
    BinCount: int = 10,
    show: bool = False,
    verbose: bool = False,
    axis: bool = False,
) -> dict:
    """
    compute stats for PointData, CellData and FieldData

    returns a dict
    """
    datadict = resultinfo(input, ignored_keys, verbose)
    if verbose:
        print(f"resultStats[{name}]: datadict={datadict}", flush=True)

    if axis:
        PointData_keys = list(datadict["PointData"]["Arrays"].keys())

        spreadSheetView = CreateView("SpreadSheetView")
        cellCenters1Display = Show(input, spreadSheetView, "SpreadSheetRepresentation")
        spreadSheetView.Update()

        filename = f"{basedir}/{name}-Axi-cellcenters-all.csv"
        ExportView(
            filename,
            view=spreadSheetView,
            RealNumberNotation="Scientific",
        )

        csv = pd.read_csv(filename)
        keys = csv.columns.values.tolist()
        # print(f"read_csv: csv={filename}, keys={keys}", flush=True)

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
                    if axis:
                        Components = kdata["Components"]
                        bounds = kdata["Bounds"]
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
                                "Minimum": [
                                    convert_data(units, bounds[0][0], fieldname)
                                ],
                                "M1": [0],
                                "Maximum": [
                                    convert_data(units, bounds[0][1], fieldname)
                                ],
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
                                    units[f"M{order}"] = [
                                        in_unit**order,
                                        out_unit**order,
                                    ]
                                else:
                                    units[f"M{order}"] = [
                                        in_unit**order,
                                        in_unit**order,
                                    ]

                                res = (
                                    2
                                    * pi
                                    * csv["Area"]
                                    * csv["Points_0"]
                                    * csv[key] ** order
                                )
                                value = res.sum() / AreaorVolume
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

                    else:
                        for excluded in fieldunits[fieldname]["Exclude"]:
                            if excluded in name:
                                found = True
                                if verbose:
                                    print(f"ignore block: {name}", flush=True)
                                break

                        if not found:
                            Components = kdata["Components"]
                            bounds = kdata["Bounds"]
                            if bounds[0][0] != bounds[0][1]:
                                if not "Stats" in kdata:
                                    # print(f"\t{key}: create kdata[Stats]", flush=True)
                                    kdata["Stats"] = {}

                                kdata["Stats"] = getresultStats(
                                    input,
                                    name,
                                    key,
                                    AttributeMode,
                                    basedir,
                                    verbose=verbose,
                                )

                                if histo:
                                    getresultHisto(
                                        input,
                                        name,
                                        dim,
                                        AreaorVolume,
                                        fieldunits,
                                        key,
                                        TypeMode,
                                        basedir,
                                        Components,
                                        BinCount=BinCount,
                                        show=show,
                                        verbose=verbose,
                                    )

    if axis:
        # remove: f"{basedir}/{name}-Axi-cellcenters-all.csv"
        os.remove(filename)
    # display stats
    return datadict
