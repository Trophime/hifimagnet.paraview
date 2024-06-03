import pandas as pd
import os
import gc
import matplotlib.pyplot as plt
from matplotlib.pyplot import hist

from typing import List

from paraview.simple import (
    Delete,
    CreateView,
    Show,
    ExportView,
)

from .method import convert_data, resultinfo


# plot with matplotlib
def plotHistoAxi(
    filename: str,
    name: str,
    key: str,
    fieldunits: dict,
    basedir: str,
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
    # csv[key] = out_values
    csv[key] = [f"{val:.3f}" for val in out_values]

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
    # ax.xaxis.set_major_formatter(lambda x, pos: f"{x:.3f}")
    show = False
    if show:
        plt.show()
    else:
        plt.tight_layout()
        plt.savefig(
            f'{basedir}/histograms/{name}-{key.replace("_Magnitude", "")}-histogram-matplotlib.png',
            dpi=300,
        )
    plt.close()

    pass


def resultHistos(
    input,
    name: str,
    Area: float,
    fieldunits: dict,
    ignored_keys: List[str],
    basedir: str,
    BinCount: int = 10,
    printed: bool = True,
    show: bool = False,
    verbose: bool = False,
):
    """
    histogram
    """
    os.makedirs(f"{basedir}/histograms", exist_ok=True)
    print(f"resultHistos: name={name}, Area={Area}, BinCount={BinCount}", flush=True)

    spreadSheetView = CreateView("SpreadSheetView")
    cellCenters1Display = Show(input, spreadSheetView, "SpreadSheetRepresentation")
    spreadSheetView.Update()

    filename = f"{basedir}/histograms/{name}-Axi-cellcenters-all.csv"
    ExportView(
        filename,
        view=spreadSheetView,
        RealNumberNotation="Scientific",
    )

    csv = pd.read_csv(filename)
    keys = csv.columns.values.tolist()
    # print(f"read_csv: csv={filename}, keys={keys}", flush=True)

    sum = csv["AxiVolume"].sum()
    csv["AxiVol"] = csv["AxiVolume"] / sum * 100
    # print(f"Area={Area}, sum={sum}")

    # check that sum is roughtly equal to 1
    # print(f'check Sum(Fraction): {csv["Fraction of total Area [%]"].sum()}')
    eps = 1.0e-4
    error = abs(1 - csv["AxiVol"].sum() / 100.0)
    assert error <= eps, f"Check Sum(Fraction) failed (error={error} > eps={eps})"

    csv.to_csv(filename)
    # print(f'Sum(AxiVol)={csv["AxiVol"].sum()}')

    datadict = resultinfo(input, ignored_keys)
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
                            fieldunits,
                            basedir,
                            BinCount,
                            show=show,
                            verbose=verbose,
                        )

    Delete(spreadSheetView)
    del spreadSheetView

    # Force a garbage collection
    collected = gc.collect()
    if verbose:
        print(f"resultsHistos: Garbage collector: collected {collected} objects.")

    # remove: f"{basedir}/histograms/{name}-Axi-cellcenters-all.csv"
    # os.remove(filename)
