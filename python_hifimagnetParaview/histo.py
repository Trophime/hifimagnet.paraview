import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

from tabulate import tabulate

from paraview.simple import (
    CellSize,
    Histogram,
    PointDatatoCellData,
    Delete,
    CreateWriter,
)

from .method import convert_data


# plot with matplotlib
def plotHisto(
    file,
    name: str,
    key: str,
    fieldunits: dict,
    AreaorVolume: float,
    basedir: str,
    dim: int,
    show: bool = True,
    verbose: bool = False,
):

    if dim == 2:
        grandeur = "Area"
    elif dim == 3:
        grandeur = "Volume"

    ax = plt.gca()
    csv = pd.read_csv(file)

    csv = csv[["bin_extents", f"{grandeur}_total"]]
    keys = csv.columns.values.tolist()
    # print("histo before scaling", flush=True)
    # print(tabulate(csv, headers="keys", tablefmt="psql"))

    # get key unit
    if verbose:
        print(f"plotHisto: file={file}, key={key}", flush=True)
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
    values = csv["bin_extents"].to_list()
    out_values = convert_data(units, values, fieldname)
    # csv = csv.assign(bin_extents=out_values)
    csv = csv.assign(bin_extents=np.array([f"{val:.2E}" for val in out_values], float))
    csv[f"{grandeur}_total"] = csv[f"{grandeur}_total"] / AreaorVolume * 100

    csv.plot.bar(
        x="bin_extents",
        y=f"{grandeur}_total",
        xlabel=rf"{msymbol}[{out_unit:~P}]",
        ylabel=f"Fraction of total {grandeur} [%]",
        title=f"{name}: {key}",
        grid=True,
        legend=False,
        rot=45,
        ax=ax,
    )

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
            f"{basedir}/histograms/{name}-{key}-histogram-matplotlib.png", dpi=300
        )
    plt.close()

    # rename columns for tabulate
    csv.rename(
        columns={
            "bin_extents": rf"{symbol} [{out_unit:~P}]",
            f"{grandeur}_total": f"Fraction of total {grandeur} [%]",
        },
        inplace=True,
    )
    if verbose:
        print(
            tabulate(csv, headers="keys", tablefmt="psql", showindex=False), flush=True
        )

    # check that sum is roughtly equal to 1
    # print(f'check Sum(Fraction): {csv[f"Fraction of total {grandeur} [%]"].sum()}')
    eps = 1.0e-4
    error = abs(1 - csv[f"Fraction of total {grandeur} [%]"].sum() / 100.0)
    if error > eps:
        print(
            f"histogram name={name}, key={key}: Check Sum(Fraction) failed : error={error}, eps={eps}",
            flush=True,
        )
    # assert error <= eps, f"Check Sum(Fraction) failed : error={error}, eps={eps}"

    csv.to_csv(f"{basedir}/histograms/{name}-{key}-histogram-matplotlib.csv")
    pass


def getresultHisto(
    input,
    name: str,
    dim: int,
    AreaorVolume: float,
    fieldunits: dict,
    key: str,
    TypeMode: str,
    basedir: str,
    Components: int = 1,
    BinCount: int = 10,
    printed: bool = True,
    show: bool = False,
    verbose: bool = False,
):
    """
    histogram
    """
    os.makedirs(f"{basedir}/histograms", exist_ok=True)
    print(
        f"getresultHisto: name={name}, key={key}, TypeMode={TypeMode}, Components={Components}, BinCount={BinCount}, show={show}",
        flush=True,
    )

    # convert pointdata to celldata
    if TypeMode == "POINT":
        pointDatatoCellData = PointDatatoCellData(
            registrationName="pointDatatoCellData", Input=input
        )

    elif TypeMode == "CELL":
        pointDatatoCellData = input

    else:
        print(
            f'resultHisto: not applicable for {key["Name"]} - unsupported data type {key["Type"]}',
            flush=True,
        )
        return

    cellSize1 = CellSize(registrationName="CellSize1", Input=pointDatatoCellData)
    # Properties modified on cellSize1 for 3D
    cellSize1.ComputeVertexCount = 0
    cellSize1.ComputeLength = 0
    if dim == 2:
        cellSize1.ComputeArea = 1  # for 2D
    elif dim == 3:
        cellSize1.ComputeArea = 0
        cellSize1.ComputeVolume = 1  # for 3D
    cellSize1.ComputeSum = 0
    cellSize1.UpdatePipeline()
    # print(f"cellSize1 CellData: {cellSize1.CellData[:]}", flush=True)

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
            print(f"Histogram: {prop}={histogram1.GetPropertyValue(prop)}", flush=True)
    # TODO from key range
    # histogram1.CustomBinRanges = [min, max]

    # Properties modified on histogram1
    histogram1.CalculateAverages = 1

    filename = f"{basedir}/histograms/{name}-{key}-histogram.csv"
    export = CreateWriter(filename, proxy=histogram1)
    if not printed:
        for prop in export.ListProperties():
            print(f"export: {prop}={export.GetPropertyValue(prop)}", flush=True)
    export.UpdateVTKObjects()  # is it needed?
    export.UpdatePipeline()
    Delete(histogram1)
    del histogram1

    plotHisto(
        filename,
        name,
        key,
        fieldunits,
        AreaorVolume,
        basedir,
        dim,
        show=show,
        verbose=verbose,
    )

    # remove temporary csv files
    # os.remove(filename)

    if TypeMode == "POINT":
        Delete(pointDatatoCellData)
        del pointDatatoCellData
