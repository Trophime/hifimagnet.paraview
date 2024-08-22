import gc
import pandas as pd

from tabulate import tabulate
from math import pi, sqrt

from paraview.simple import (
    CellSize,
    ExtractBlock,
    CellDatatoPointData,
    Calculator,
    MergeBlocks,
    Delete,
    PointDatatoCellData,
    CellCenters,
)
from paraview import servermanager as sm
from paraview.vtk.numpy_interface import dataset_adapter as dsa
from paraview.vtk.numpy_interface import algorithms as algs

from .method import convert_data, info, resultinfo, momentN, integrateKeys, keyinfo
from .statsAxi import resultStats, createStatsTable
from .histoAxi import resultHistos
from .meshinfo import createVectorNorm


def part_integrate(
    input,
    name: str,
    selected_blocks: list[str],
    basedir: str,
    merge: bool = True,
    verbose: bool = False,
) -> pd.DataFrame:
    """compute integral over input

    Args:
        input: paraview reader
        name (str): block name
        selected_blocks (list[str]): block selection
        basedir (str): result directory
        merge (bool, optional): merge selected blocks. Defaults to True.
        verbose (bool, optional): print verbose. Defaults to False.

    Returns:
        pd.DataFrame: integral dataframe
    """
    print(f"part_integrate: name={name}", flush=True)

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
    for field in tmp.PointData:
        fname = f'"{field.Name}"'
        # print(f"fname={fname}, components={field.GetNumberOfComponents()}",flush=True)
        if field.GetNumberOfComponents() > 1:
            fname = f'mag("{field.Name}")'
            for i in range(field.GetNumberOfComponents()):
                drop_keys.append(f"{field.Name}_{i}")
            drop_keys.append(f"{field.Name}_Magnitude")
        else:
            drop_keys.append(field.Name)

        for order in range(1, 5):
            tmp = momentN(tmp, fname, field.Name, order, "Point Data")

    # Add field for AxiVol
    tmp = Calculator(registrationName=f"AxiVol", Input=tmp)
    tmp.AttributeType = "Point Data"
    tmp.ResultArrayName = f"AxiVol"
    tmp.Function = "coordsX"
    tmp.UpdatePipeline()

    # print(f"tmp: {tmp.PointData[:]}", flush=True)

    # convert CellData to PointData
    pointDatatoCellData = PointDatatoCellData(
        registrationName="PointDatatoCellData", Input=tmp
    )
    tmp = pointDatatoCellData
    # print(f"tmp pointDatatoCellData: {tmp.CellData[:]}", flush=True)

    # get integrated values
    csvfile = integrateKeys(tmp, name, basedir, verbose=verbose)
    csv = pd.read_csv(csvfile)
    # print(f"integrateKeys: csv={list(csv.keys())}",flush=True)
    drop_keys += ["Cell ID", "Area", "Cell Type"]
    # print(f"drop_keys={drop_keys}", flush=True)
    csv.drop(columns=drop_keys, inplace=True)

    # divide columns by AxiVol
    keys = list(csv.keys())  # .remove("AxiVol")
    for key in csv.keys():
        if key != "AxiVol":
            csv[key] = csv[key].div(csv["AxiVol"], axis=0)

    if verbose:
        print(
            tabulate(csv.transpose(), headers="keys", tablefmt="psql"),
            flush=True,
        )

    Delete(cellDatatoPointData)
    del cellDatatoPointData

    # Force a garbage collection
    collected = gc.collect()
    if verbose:
        print(
            f"resultsHistos: Garbage collector: collected {collected} objects.",
            flush=True,
        )

    return csv


def part(
    pinput,
    name: str,
    fieldunits: dict,
    ignored_keys: list[str],
    ureg,
    basedir: str,
    ComputeHisto: bool,
    BinCount: int = 20,
    show: bool = False,
    verbose: bool = False,
):
    """stats & histos for a part

    Args:
        pinput: paraview reader
        name (str): block name
        fieldunits (dict): dict of diel units
        ignored_keys (list[str]): list of ignored fields
        ureg: pint unit registry
        basedir (str): result directory
        ComputeHisto (bool): compute histograms
        BinCount (int, optional): number of bins in histogram. Defaults to 20.
        show (bool, optional): show histogramms. Defaults to False.
        verbose (bool, optional): print verbose. Defaults to False.

    Returns:
        vol, statsdict
    """
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
    # print(f"type(dataset)={type(dataset)}", flush=True)
    if dataset.IsA("vtkUnstructuredGrid"):
        # print("UnstructuredGrid", flush=True)
        block = dataset  # .GetBlock(0)
    elif dataset.IsA("vtkPolyData"):
        # print("vtkPolyData", flush=True)
        block = dataset  # .GetBlock(0)
    elif dataset.IsA("vtkMultiBlockDataSet"):
        # print("MultiBlockDataSet", flush=True)
        block = dataset.GetBlock(0)

    Area_data = block.GetPointData().GetArray("AxiVolume")
    # print(f'block cellData[Area_data]: {Area_data}, type={type(Area_data)}',flush=True)
    Area = block.GetPointData().GetArray("AxiVolume")
    # print(f'block fieldData[Area]: {Area}, type={type(Area)}',flush=True)
    np_Area = np_dataset.PointData["AxiVolume"]
    # print(f'block fieldData[np_Area]: {np_Area}, type={type(np_Area)}, length={algs.shape(np_Area)}',flush=True)
    vol = algs.sum(np_Area)
    vunits = fieldunits["Volume"]["Units"]
    mm3 = f"{vunits[1]:~P}"
    vol_mm3 = convert_data(
        {"Volume": vunits},
        vol,
        "Volume",
    )
    print(
        f"{name}: block fieldData[np_Area]: vol={vol_mm3} {mm3}, parts={algs.shape(np_Area)}",
        flush=True,
    )

    # # check tvol == Sum(vol)
    # if abs(1 - vol / tvol) > 1.0e-3:
    #     print(f"insert Total volume != vol(insert), tvol={tvol}, vol={vol}, error={abs(1-vol/tvol)*100} %",flush=True)

    statsdict = resultStats(
        calculator1,
        name,
        2,
        vol,
        fieldunits,
        ignored_keys,
        ureg,
        basedir,
        verbose=verbose,
    )
    # print(f"insert statsdict: {statsdict}", flush=True)
    if ComputeHisto:
        resultHistos(
            calculator1,
            name,
            vol,
            fieldunits,
            ignored_keys,
            basedir,
            BinCount=BinCount,
            show=show,
            verbose=verbose,
        )

    Delete(pointDatatoCellData)
    del pointDatatoCellData

    # Force a garbage collection
    collected = gc.collect()
    if verbose:
        print(f"part: Garbage collector: collected {collected} objects.", flush=True)

    return vol, statsdict


def cylField(input, key: str, nkey: str, AttributeType: str):
    """compute r and theta component of a vector Field

    vr = ux
    vz = uy


    ux: "thermo_electric.electric.current_density_X"

    Args:
        input: paraview reader
        key (str): field name
        nkey (str): _description_
        AttributeType (str): 'Point Data'

    Raises:
        RuntimeError: cylField: key - unsupported AttributeType

    Returns:
        cylindric paraview reader
    """

    if AttributeType != "Point Data":
        raise RuntimeError(
            f"cylField: {key} - unsupported AttributeType: {AttributeType}"
        )
    # inputDataPoint = [field.Name for field in input.PointData]
    # print(
    #     f"cylField: input PointData = {inputDataPoint}",
    #     flush=True,
    # )

    calculator2 = Calculator(Input=input)
    calculator2.AttributeType = AttributeType  # 'Cell Data'
    calculator2.ResultArrayName = f"{key}_r"
    calculator2.Function = f'"{key}_X"'
    calculator2.UpdatePipeline()

    calculator3 = Calculator(Input=calculator2)
    calculator3.AttributeType = AttributeType  # 'Cell Data'
    calculator3.ResultArrayName = f"{key}_z"
    calculator3.Function = f'"{key}_Y"'

    calculator3.UpdatePipeline()
    return calculator3


def meshinfo(
    input,
    dim: int,
    fieldunits: dict,
    ignored_keys: list[str],
    basedir: str,
    ureg,
    ComputeStats: bool = True,
    ComputeHisto: bool = False,
    BinCount: int = 10,
    show: bool = False,
    verbose: bool = False,
    printed: bool = True,
) -> tuple:
    """display geometric info from input dataset

    Args:
        input (_type_): paraview reader
        dim (int): geometry dimmension
        fieldunits (dict): dictionnary of field units
        ignored_keys (list[str]): list of ignored keys
        basedir (str): result directory
        ureg: pint unit registry
        ComputeStats (bool, optional): compute statistics. Defaults to True.
        ComputeHisto (bool, optional): compute histograms. Defaults to False.
        BinCount (int, optional): number of bins in histograms. Defaults to 10.
        show (bool, optional): show histograms. Defaults to False.
        verbose (bool, optional): print verbose. Defaults to False.
        printed (bool, optional): Defaults to True.

    Returns:
        cellsize: updated paraview reader
        blockdata (dict): dict of blocks data
        stats (dict): dict of statistics
    """

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
    cellDatatoPointData1 = CellDatatoPointData(
        registrationName="CellDatatoPointData", Input=input
    )

    # for vector
    print("Add Norm for vectors and CylFields:", flush=True)
    calculator = cellDatatoPointData1

    for field in cellDatatoPointData1.PointData:
        if (dim == 2 and field.GetNumberOfComponents() > 1) or (
            field.GetNumberOfComponents() == dim
        ):
            print(
                f"create {field.Name}norm for {field.Name} PointData vector",
                flush=True,
            )
            calculator = createVectorNorm(
                calculator, field.Name, field.Name, "Point Data"
            )

            print(
                f"create {field.Name}_r and {field.Name}_z for {field.Name} PointData vector",
                flush=True,
            )
            calculator = cylField(calculator, field.Name, field.Name, "Point Data")

    # PointData to CellData
    pointDatatoCellData = PointDatatoCellData(
        registrationName="PointDatatoCellData", Input=calculator
    )

    print("Get mesh size", flush=True)
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
            print(f"cellsize: {prop}={cellsize.GetPropertyValue(prop)}", flush=True)

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
    # print(f"type(dataset)={type(dataset)}", flush=True)
    if dataset.IsA("vtkUnstructuredGrid"):
        print("UnstructuredGrid", flush=True)
        block = dataset
    elif dataset.IsA("vtkMultiBlockDataSet"):
        print("MultiBlockDataSet", flush=True)
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
            child_info = input.GetSubsetDataInformation(0, child)  # rootChild)
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
        resultinfo(calculator1, ignored_keys, verbose)
        Delete(pointDatatoCellData)
        del pointDatatoCellData

        # Force a garbage collection
        collected = gc.collect()
        if verbose:
            print(f"meshinfo: Garbage collector: collected {collected} objects.")

        print("Data ranges without Air:", flush=True)
        selected_blocks = [block for block in blockdata.keys() if not "Air" in block]

        # TODO: need to restart from input and apply ... Calculator1
        extractBlock1 = ExtractBlock(registrationName="insert", Input=input)
        extractBlock1.Selectors = selected_blocks
        extractBlock1.UpdatePipeline()
        mergeBlocks1 = MergeBlocks(registrationName="MergeBlocks1", Input=extractBlock1)
        mergeBlocks1.UpdatePipeline()

        vol, statsdict = part(
            mergeBlocks1,
            "insert",
            fieldunits,
            ignored_keys,
            ureg,
            basedir,
            ComputeHisto,
            BinCount,
        )
        stats.append(statsdict)

        icsv = part_integrate(
            input, "insert", selected_blocks, basedir, merge=True, verbose=verbose
        )
        if verbose:
            print(f'insert: vol={vol}, ivol={icsv["AxiVol"].to_list()[0] * 2 * pi}')
        for key, value in statsdict.items():
            for array, avalue in value["Arrays"].items():
                if "Stats" in avalue:
                    (toolbox, physic, fieldname) = keyinfo(array)

                    symbol = fieldunits[fieldname]["Symbol"]
                    msymbol = symbol
                    if "mSymbol" in fieldunits[fieldname]:
                        msymbol = fieldunits[fieldname]["mSymbol"]
                    units = {fieldname: fieldunits[fieldname]["Units"]}
                    [in_unit, out_unit] = fieldunits[fieldname]["Units"]

                    M = [0] * 5
                    M[1] = icsv[f"{array}_moment1"].to_list()[0]
                    M[2] = icsv[f"{array}_moment2"].to_list()[0]
                    M[3] = icsv[f"{array}_moment3"].to_list()[0]
                    M[4] = icsv[f"{array}_moment4"].to_list()[0]
                    avalue["Stats"]["Mean"] = convert_data(units, M[1], fieldname)

                    if units[fieldname][0] != ureg.kelvin:
                        avalue["Stats"]["Standard Deviation"] = convert_data(
                            units, sqrt(abs(M[1] * M[1] - M[2])), fieldname
                        )
                    else:
                        avalue["Stats"]["Standard Deviation"] = sqrt(
                            abs(M[1] * M[1] - M[2])
                        )

                    for order in range(2, 5):
                        units = {f"M{order}": fieldunits[fieldname]["Units"]}
                        if in_unit != ureg.kelvin:
                            units[f"M{order}"] = [in_unit * order, out_unit * order]
                        else:
                            units[f"M{order}"] = [in_unit * order, in_unit * order]

                        avalue["Stats"][f"M{order}"] = convert_data(
                            units, M[order], f"M{order}"
                        )

                    # print(f"{array}: {avalue['Stats']}")

        stats.append(statsdict)

        # aggregate stats data
        createStatsTable([statsdict], "insert", fieldunits, basedir, verbose)

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

            vol, statsdict = part(
                extractBlock1,
                name,
                fieldunits,
                ignored_keys,
                ureg,
                basedir,
                ComputeHisto,
                BinCount,
            )
            stats.append(statsdict)
            sum_vol += vol

            icsv = part_integrate(
                input, name, [block], basedir, merge=False, verbose=verbose
            )
            if verbose:
                print(
                    f'name={name}: vol={vol}, ivol={icsv["Points_0"].to_list()[0] * 2 * pi}'
                )

            # print(f'insert: vol={vol}, ivol={icsv["AxiVol"].to_list()[0] * 2 * pi}')
            for key, value in statsdict.items():
                for array, avalue in value["Arrays"].items():
                    if "Stats" in avalue:
                        (toolbox, physic, fieldname) = keyinfo(array)

                        symbol = fieldunits[fieldname]["Symbol"]
                        msymbol = symbol
                        if "mSymbol" in fieldunits[fieldname]:
                            msymbol = fieldunits[fieldname]["mSymbol"]
                        units = {fieldname: fieldunits[fieldname]["Units"]}
                        [in_unit, out_unit] = fieldunits[fieldname]["Units"]

                        M = [0] * 5
                        M[1] = icsv[f"{array}_moment1"].to_list()[0]
                        M[2] = icsv[f"{array}_moment2"].to_list()[0]
                        M[3] = icsv[f"{array}_moment3"].to_list()[0]
                        M[4] = icsv[f"{array}_moment4"].to_list()[0]
                        avalue["Stats"]["Mean"] = convert_data(units, M[1], fieldname)

                        if units[fieldname][0] != ureg.kelvin:
                            avalue["Stats"]["Standard Deviation"] = convert_data(
                                units, sqrt(abs(M[1] * M[1] - M[2])), fieldname
                            )
                        else:
                            avalue["Stats"]["Standard Deviation"] = sqrt(
                                abs(M[1] * M[1] - M[2])
                            )

                        for order in range(2, 5):
                            units = {f"M{order}": fieldunits[fieldname]["Units"]}
                            if in_unit != ureg.kelvin:
                                units[f"M{order}"] = [in_unit * order, out_unit * order]
                            else:
                                units[f"M{order}"] = [in_unit * order, in_unit * order]

                            avalue["Stats"][f"M{order}"] = convert_data(
                                units, M[order], f"M{order}"
                            )

                        print(f"{array}: {avalue['Stats']}")

            stats.append(statsdict)
            Delete(extractBlock1)
            del extractBlock1

            # Force a garbage collection
            collected = gc.collect()
            if verbose:
                print(f"loopblock: Garbage collector: collected {collected} objects.")

            # aggregate stats data
            createStatsTable([statsdict], name, fieldunits, basedir, verbose)

        # check tvol == Sum(vol)
        print(
            f"reworked Total volume: tvol={tvol}, sum_vol={sum_vol}, error={abs(1-sum_vol/tvol)*100} %"
        )
        if abs(1 - sum_vol / tvol) > 1.0e-3:
            print(
                f"reworked Total volume != Sum(vol), tvol={tvol}, sum_vol={sum_vol}, error={abs(1-sum_vol/tvol)*100} %"
            )

        # aggregate stats data
        createStatsTable(stats, "total", fieldunits, basedir, verbose)

    if dataset.IsA("vtkUnstructuredGrid"):
        stats = []

        print("Data ranges:", flush=True)
        resultinfo(calculator1, ignored_keys, verbose)
        Delete(pointDatatoCellData)
        del pointDatatoCellData

        # Force a garbage collection
        collected = gc.collect()
        if verbose:
            print(f"meshinfo: Garbage collector: collected {collected} objects.")

        print("Data ranges without Air:", flush=True)
        selected_blocks = [block]

        vol, statsdict = part(
            calculator1,
            "insert",
            fieldunits,
            ignored_keys,
            ureg,
            basedir,
            ComputeHisto,
            BinCount,
        )
        stats.append(statsdict)

        icsv = part_integrate(
            input, "insert", selected_blocks, basedir, merge=True, verbose=verbose
        )
        if verbose:
            print(f'insert: vol={vol}, ivol={icsv["AxiVol"].to_list()[0] * 2 * pi}')
        for key, value in statsdict.items():
            for array, avalue in value["Arrays"].items():
                if "Stats" in avalue:
                    (toolbox, physic, fieldname) = keyinfo(array)

                    symbol = fieldunits[fieldname]["Symbol"]
                    msymbol = symbol
                    if "mSymbol" in fieldunits[fieldname]:
                        msymbol = fieldunits[fieldname]["mSymbol"]
                    units = {fieldname: fieldunits[fieldname]["Units"]}
                    [in_unit, out_unit] = fieldunits[fieldname]["Units"]

                    M = [0] * 5
                    M[1] = icsv[f"{array}_moment1"].to_list()[0]
                    M[2] = icsv[f"{array}_moment2"].to_list()[0]
                    M[3] = icsv[f"{array}_moment3"].to_list()[0]
                    M[4] = icsv[f"{array}_moment4"].to_list()[0]
                    avalue["Stats"]["Mean"] = convert_data(units, M[1], fieldname)

                    avalue["Stats"]["Standard Deviation"] = convert_data(
                        units, sqrt(abs(M[1] * M[1] - M[2])), fieldname
                    )

                    for order in range(2, 5):
                        units = {f"M{order}": fieldunits[fieldname]["Units"]}
                        if in_unit != ureg.kelvin:
                            units[f"M{order}"] = [in_unit * order, out_unit * order]
                        else:
                            units[f"M{order}"] = [in_unit * order, in_unit * order]

                        avalue["Stats"][f"M{order}"] = convert_data(
                            units, M[order], f"M{order}"
                        )

                    # print(f"{array}: {avalue['Stats']}")

        stats.append(statsdict)

        # aggregate stats data
        if ComputeStats:
            createStatsTable([statsdict], "insert", fieldunits, basedir, verbose)

    return input, blockdata, stats
