import pandas as pd
import numpy as np
import os
import gc
import matplotlib.pyplot as plt
from matplotlib.pyplot import hist

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
    """plot histogramms

    Args:
        filename (str): csv file containing datas for hist
        name (str): block name (aka `feelpp` marker) / insert
        key (str): field name
        fieldunits (dict): dict field units
        basedir (str): result directory
        BinCount (int): number of bins in histogram
        show (bool, optional): show histogramms. Defaults to True.
        verbose (bool, optional): print verbose. Defaults to False.
    """
    print(f"plotHistAxi: name={name}, key={key}, bin={BinCount}", flush=True)

    ax = plt.gca()
    csv = pd.read_csv(filename)
    keys = csv.columns.values.tolist()
    # print(f"plotHistAxi: keys={keys}", flush=True)
    # print("histo before scaling",flush=True)
    # print(tabulate(csv, headers="keys", tablefmt="psql"))

    # get key unit
    keyinfo = key.replace("_Magnitude", "").split(".")
    # print(f"keyinfo={keyinfo}", flush=True)
    if len(keyinfo) == 1:
        fieldname = key.replace("_Magnitude", "")
    elif len(keyinfo) == 2:
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
    csv[key] = [f"{val:.2E}" for val in out_values]

    counts, extend_bins, patches = hist(
        np.array(csv[key], float),
        bins=BinCount,
        weights=csv["AxiVol"],
        rwidth=0.5,
    )
    print(f"counts={counts}", flush=True)
    print(f"extend_bins={extend_bins}", flush=True)
    print(f"patches={patches}", flush=True)

    ticks = [(patch._x0 + patch._x0 + patch._width) / 2 for patch in patches]
    total_key = "Fraction of total Volume [%]"
    plt.xlabel(rf"{msymbol}[{out_unit:~P}]")
    plt.ylabel(total_key)
    plt.xticks(ticks, rotation=45, ha="right")
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

    df_histo_plt = pd.DataFrame()
    df_histo_plt[rf"{symbol} [{out_unit:~P}]"] = ticks
    df_histo_plt["Fraction of total Volume [%]"] = counts
    df_histo_plt.to_csv(
        f"{basedir}/histograms/{name}-{key.replace('_Magnitude', '')}-histogram-matplotlib.csv"
    )

    pass


def resultHistos(
    input,
    name: str,
    Area: float,
    fieldunits: dict,
    ignored_keys: list[str],
    basedir: str,
    BinCount: int = 10,
    printed: bool = True,
    show: bool = False,
    verbose: bool = False,
):
    """histogramms

    Args:
        input: paraview reader
        name (str): block name (aka `feelpp` marker) / insert
        Area (float): total area
        fieldunits (dict): dictionnary of field units
        ignored_keys (list[str]): list of ignored keys
        basedir (str): result directory
        BinCount (int, optional): number of bins in histograms. Defaults to 10.
        printed (bool, optional): Defaults to True.
        show (bool, optional): show histogramms. Defaults to False.
        verbose (bool, optional): print verbose. Defaults to False.
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
    # print(f"Area={Area}, sum={sum}",flush=True)

    # check that sum is roughtly equal to 1
    # print(f'check Sum(Fraction): {csv["Fraction of total Area [%]"].sum()}',flush=True)
    eps = 1.0e-4
    error = abs(1 - csv["AxiVol"].sum() / 100.0)
    assert error <= eps, f"Check Sum(Fraction) failed (error={error} > eps={eps})"

    csv.to_csv(filename)
    # print(f'Sum(AxiVol)={csv["AxiVol"].sum()}',flush=True)

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
        print(
            f"resultsHistos: Garbage collector: collected {collected} objects.",
            flush=True,
        )

    # remove: f"{basedir}/histograms/{name}-Axi-cellcenters-all.csv"
    # os.remove(filename)
