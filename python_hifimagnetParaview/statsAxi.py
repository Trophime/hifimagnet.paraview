import pandas as pd
import os

from tabulate import tabulate
from typing import List
from math import pi, sqrt

from paraview.simple import (
    ExportView,
    CreateView,
    Show,
)

from .method import convert_data, resultinfo


def createStatsTable(
    stats: list, name: str, fieldunits: dict, basedir: str, verbose: bool = False
):

    # TODO add a column with the Block Name
    # the column Block Name in cvs is not what I think
    # it is either "Primary Statistics" or "Derived Statistics"
    # remove "Derived Statistics" rows
    # why Standard Deviation and Variance are NaN??
    os.makedirs(f"{basedir}/stats", exist_ok=True)
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
            if len(keyinfo) == 1:
                fieldname = key
            elif len(keyinfo) == 2:
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
                df.to_csv(f"{basedir}/stats/{key}-descriptivestats-create.csv")
                dfs.append(df)

    total_df = pd.concat(dfs)
    print(tabulate(total_df, headers="keys", tablefmt="psql", showindex=False))
    total_df.to_csv(f"{basedir}/stats/{name}-descriptiveAxistats-create.csv")

    pass


def resultStats(
    input,
    name: str,
    dim: int,
    Area: float,
    fieldunits: dict,
    ignored_keys: List[str],
    ureg,
    basedir: str,
    histo: bool = False,
    BinCount: int = 10,
    show: bool = False,
    verbose: bool = False,
) -> dict:
    """
    compute stats for PointData, CellData and FieldData

    returns a dict
    """
    os.makedirs(f"{basedir}/stats", exist_ok=True)

    datadict = resultinfo(input, ignored_keys, verbose)
    if verbose:
        print(
            f'resultStats[{name}]: datadict={datadict["PointData"]["Arrays"].keys()}',
            flush=True,
        )
    PointData_keys = list(datadict["PointData"]["Arrays"].keys())

    spreadSheetView = CreateView("SpreadSheetView")
    cellCenters1Display = Show(input, spreadSheetView, "SpreadSheetRepresentation")
    spreadSheetView.Update()

    filename = f"{basedir}/stats/{name}-Axi-cellcenters-all.csv"
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
        if datatype == "PointData":
            AttributeMode = datadict[datatype]["AttributeMode"]
            TypeMode = datadict[datatype]["TypeMode"]
            for key, kdata in datadict[datatype]["Arrays"].items():
                if not key in ignored_keys:
                    Components = kdata["Components"]
                    bounds = kdata["Bounds"]

                    keyinfo = key.split(".")
                    # print(f"keyinfo={keyinfo}", flush=True)
                    if len(keyinfo) == 1:
                        fieldname = key
                    elif len(keyinfo) == 2:
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

    # remove: f"{basedir}/stats/{name}-Axi-cellcenters-all.csv"
    os.remove(filename)

    # display stats
    return datadict
