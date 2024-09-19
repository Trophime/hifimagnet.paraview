import gc

from paraview.simple import (
    CellSize,
    ExtractBlock,
    CellDatatoPointData,
    Calculator,
    MergeBlocks,
    Delete,
)
from paraview import servermanager as sm
from paraview.vtk.numpy_interface import dataset_adapter as dsa
from paraview.vtk.numpy_interface import algorithms as algs

from .method import convert_data, info, resultinfo
from .stats import resultStats, createStatsTable


def scaleField(input, key: str, nkey: str, AttributeType: str, factor: float):
    """scale a Field by 1/factor

    AttributeType:

    Args:
        input: paraview reader
        key (str): field name
        nkey (str): _description_
        AttributeType (str): 'Cell Data'
        factor (float): scale factor

    Returns:
        scaled paraview reader
    """
    calculator1 = Calculator(registrationName="Calculator1", Input=input)
    calculator1.AttributeType = AttributeType  # 'Cell Data'
    calculator1.ResultArrayName = f"{nkey}"
    calculator1.Function = f'{key}/{factor}")'

    calculator1.UpdatePipeline()
    return calculator1


def addField(input, key: str, nkey: str, AttributeType: str, factor: float):
    """add factor to Field

    Args:
        input: paraview reader
        key (str): field name
        nkey (str): _description_
        AttributeType (str): 'Cell Data'
        factor (float): add factor

    Returns:
        paraview reader with added field
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
    """compute r and theta component of a vector Field

    vr = ux * cos + uy * sin
    vt = ux* -sin + uy * cos

    with cos = x/sqrt(x²+y²) and sin = y/sqrt(x²+y²)

    x: coordsX
    y : coordsY

    ux: "thermo_electric.electric.current_density_X"

    Args:
        input: paraview reader
        key (str): field name
        nkey (str): _description_
        AttributeType (str): 'Point Data'

    Raises:
        RuntimeError: rectTocylField: key - unsupported AttributeType

    Returns:
        cylindric paraview reader
    """

    if AttributeType != "Point Data":
        raise RuntimeError(
            f"rectTocylField: {key} - unsupported AttributeType: {AttributeType}"
        )
    inputDataPoint = [field.Name for field in input.PointData]
    # print(
    #     f"rectTocylField: input PointData = {inputDataPoint}",
    #     flush=True,
    # )

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
    """create Norm of a vector

    Args:
        input: paraview reader
        key (str): field name
        nkey (str): _description_
        AttributeType (str): 'Cell Data'
        printed (bool, optional): Defaults to True.

    Returns:
        paraview reader
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
        input: paraview reader
        dim (int): geometry dimmension
        fieldunits (dict): dictionnary of field units
        ignored_keys (list[str]): list of ignored fields
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

    # rectTocyl: need CellDataToPointData before
    # for temperature add, for forces and densities norm, rescale
    cellDatatoPointData1 = CellDatatoPointData(
        registrationName="CellDatatoPointData", Input=input
    )

    # for vector
    print("Add Norm for vectors and RectToCyl:", flush=True)
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
                f"create {field.Name}ur and {field.Name}ut for {field.Name} PointData vector",
                flush=True,
            )
            calculator = rectTocylField(
                calculator, field.Name, field.Name, "Point Data"
            )

    print("Get mesh size", flush=True)
    cellsize = CellSize(calculator)  # input

    # set some params
    cellsize.ComputeLength = 0
    if dim == 2:
        cellsize.ComputeArea = 1
        cellsize.ComputeVolume = 0
        grandeur = "Area"
    elif dim == 3:
        cellsize.ComputeArea = 0
        cellsize.ComputeVolume = 1
        grandeur = "Volume"
    cellsize.ComputeVertexCount = 0
    cellsize.ComputeSum = 1
    # get params list
    if not printed:
        for prop in cellsize.ListProperties():
            print(f"cellsize: {prop}={cellsize.GetPropertyValue(prop)}", flush=True)

    # apply
    cellsize.UpdatePipeline()
    dataInfo = info(cellsize)

    dataset = sm.Fetch(cellsize)
    np_dataset = dsa.WrapDataObject(dataset)

    blockdata = {}
    # check dataset type
    if dataset.IsA("vtkUnstructuredGrid"):
        print("UnstructuredGrid", flush=True)
        block = dataset
    elif dataset.IsA("vtkMultiBlockDataSet"):
        print("MultiBlockDataSet", flush=True)
        block = dataset.GetBlock(0)

    grandeur_data = block.GetCellData().GetArray(grandeur)
    # print(f'block cellData[{grandeur}_data]: {grandeur_data}, type={type(grandeur_data)}', flush=True)
    AreaorVolume = block.GetFieldData().GetArray(grandeur)
    # print(f'block fieldData[{grandeur}]: {AreaorVolume}, type={type(AreaorVolume)}', flush=True)
    np_grandeur = np_dataset.FieldData[grandeur]
    # print(f'block fieldData[np_{grandeur}]: {np_grandeur}, type={type(np_grandeur)}, length={algs.shape(np_grandeur)}', flush=True)
    tvol = algs.sum(np_grandeur)
    vunits = fieldunits[grandeur]["Units"]
    mmdim = f"{vunits[1]:~P}"
    tvol_mmdim = convert_data(
        {grandeur: vunits},
        tvol,
        grandeur,
    )
    print(
        f"block fieldData[np_{grandeur}]: total={tvol_mmdim} {mmdim}, parts={algs.shape(np_grandeur)}",
        flush=True,
    )

    if dataset.IsA("vtkMultiBlockDataSet"):
        hierarchy = dataInfo.GetHierarchy()
        rootnode = hierarchy.GetRootNode()
        rootSelector = f"/{hierarchy.GetRootNodeName()}"
        blocks = hierarchy.GetNumberOfChildren(rootnode)
        print(f"Load blocks: {blocks}", flush=True)

        sum_vol = 0
        for i in range(blocks):
            child = hierarchy.GetChild(rootnode, i)
            name = hierarchy.GetNodeName(child)
            rootChild = f"{rootSelector}/{name}"
            child_info = input.GetSubsetDataInformation(0, child)  # rootChild)
            # print(f'block[{i}]: {name}, child={child}, childInfo: {child_info}', flush=True)

            nodes = child_info.GetNumberOfPoints()
            cells = child_info.GetNumberOfCells()
            vol = np_grandeur.Arrays[i][0]
            vol_mmdim = convert_data(
                {grandeur: vunits},
                vol,
                grandeur,
            )
            print(
                f"block[{i}]: {name}, nodes={nodes}, cells={cells}, vol={vol_mmdim} {mmdim}",
                flush=True,
            )

            blockdata[rootChild] = {
                "name": name,
                "nodes": nodes,
                "cells": cells,
                grandeur: vol,
            }

            sum_vol += vol

        # check tvol == Sum(vol)
        if abs(1 - sum_vol / tvol) > 1.0e-3:
            raise RuntimeError(
                f"Total {grandeur} != Sum({grandeur}), error={abs(1-sum_vol/tvol)}"
            )

        # Compute Stats
        stats = []

        print("Data ranges:", flush=True)
        datadict = resultinfo(cellsize, ignored_keys, verbose)

        """
        # To speed up Stats by blocks
        # !!! add Standard Deviation form Mean and M2 !!!
        # this would replace also createStatsTable

        print("Global Stats:", flush=True)

        for datatype in datadict:
            if datatype != "FieldData":
                AttributeMode = datadict[datatype]["AttributeMode"]
                TypeMode = datadict[datatype]["TypeMode"]
        
                statistics = DescriptiveStatistics(cellsize)
                statistics.VariablesofInterest = [
                    key for key in datadict["PointData"]["Arrays"] if not key in ignored_keys
                ]
                print(
                f"statistics.VariablesofInterest={statistics.VariablesofInterest}",
                    flush=True,
                )
                statistics.AttributeMode = AttributeMode
                spreadSheetView = CreateView("SpreadSheetView")
                descriptiveStatisticsDisplay = Show(
                    statistics, spreadSheetView, "SpreadSheetRepresentation"
                )
                spreadSheetView.Update()
                export = ExportView(
                    f"total-descriptivestats-ttt.csv",
                    view=spreadSheetView,
                    RealNumberNotation="Scientific",
                )
                csv = pd.read_csv(f"total-descriptivestats-ttt.csv")
                csv.rename(columns={"Block Name": "BlockName"}, inplace=True)
                csv = csv[(csv.BlockName != "Derived Statistics")]
                dropped_keys = [
                    "BlockName",
                    "Cardinality",
                    "Kurtosis",
                    "Skewness",
                    "Sum",
                    "Variance",
                ]
                csv.drop(columns=dropped_keys, inplace=True)
                
                # Compute Standard Deviation
                csv['Standard Deviation'] = sqrt(abs(csv["Mean"]**2 - csv["M2"]))
                    
                # Add Name 
                (nrows, ncols) = csv.shape
                csv['Name'] = [name for i in range(nrows)]

                # Reorder columns
                csv = csv[
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
                print(f"total stats: key={list(csv.keys())}", flush=True)
                print(
                    tabulate(
                        csv,
                        headers="keys",
                        tablefmt="psql", 
                        showindex=False
                    )
                )
                    
                # convert data to requested units, 
                # add units for Variable
                # set values with appropriate units

                # per blocks
                for i, block in enumerate(blockdata.keys()):
                    name = blockdata[block]["name"]
                    print(f"block[{i}]: extract {block}, name={name}", flush=True)
                    descriptiveStatisticsDisplay.BlockVisibilities = [block]
                    spreadSheetView.Update()
                    export = ExportView(
                        f"{name}-descriptivestats-ttt.csv",
                        view=spreadSheetView,
                        RealNumberNotation="Scientific",
                    )
                    csv = pd.read_csv(f"{name}-descriptivestats-ttt.csv")
                    csv.rename(columns={"Block Name": "BlockName"}, inplace=True)
                    csv = csv[(csv.BlockName != "Derived Statistics")]
                    dropped_keys = [
                        "BlockName",
                        "Cardinality",
                        "Kurtosis",
                        "Skewness",
                        "Sum",
                        "Variance",
                    ]
                    csv.drop(columns=dropped_keys, inplace=True)
                    
                    # Compute Standard Deviation
                    csv['Standard Deviation'] = sqrt(abs(csv["Mean"]**2 - csv["M2"]))
                    
                    # Add Name 
                    (nrows, ncols) = csv.shape
                    csv['Name'] = [name for i in range(nrows)]

                    # Reorder columns
                    csv = csv[
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
                    print(f"{name} stats: key={list(csv.keys())}", flush=True)
                    print(
                        tabulate(
                            csv,
                            headers="keys",
                            tablefmt="psql", 
                            showindex=False
                        )
                    )
                
                    # convert data to requested units, 
                    # add units for Variable
                    # set values with appropriate units

                Delete(descriptiveStatisticsDisplay)
                Delete(spreadSheetView)
                del spreadSheetView
                Delete(statistics)
                del statistics

        # concat csv per datatype
                

        exit(1)
        """

        print("Data ranges without Air:", flush=True)
        extractBlock1 = ExtractBlock(registrationName="insert", Input=cellsize)
        extractBlock1.Selectors = [
            block for block in blockdata.keys() if not "Air" in block
        ]
        extractBlock1.UpdatePipeline()
        Grandeurs = [blockdata[block][grandeur] for block in extractBlock1.Selectors]
        mergeBlocks1 = MergeBlocks(registrationName="MergeBlocks1", Input=extractBlock1)
        mergeBlocks1.UpdatePipeline()
        statsdict = resultStats(
            mergeBlocks1,
            "insert",
            dim,
            sum(Grandeurs),
            fieldunits,
            ignored_keys,
            ureg,
            basedir,
            histo=ComputeHisto,
            BinCount=BinCount,
            show=show,
            verbose=verbose,
        )
        if verbose:
            print(f"insert statsdict={statsdict}", flush=True)
        stats.append(statsdict)
        extractBlock1.UpdatePipeline()
        Delete(extractBlock1)
        del extractBlock1

        # Force a garbage collection
        collected = gc.collect()
        if verbose:
            print(f"Garbage collector: collected {collected} objects.", flush=True)

        # aggregate stats data
        createStatsTable([statsdict], "insert", fieldunits, basedir, ureg, verbose)

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
                    dim,
                    blockdata[block][grandeur],
                    fieldunits,
                    ignored_keys,
                    ureg,
                    basedir,
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
                    print(
                        f"Garbage collector: collected {collected} objects.", flush=True
                    )

                # aggregate stats data
                createStatsTable([statsdict], name, fieldunits, basedir, ureg, verbose)

            # aggregate stats data
            createStatsTable(stats, "total", fieldunits, basedir, ureg, verbose)

    elif dataset.IsA("vtkUnstructuredGrid"):
        stats = []

        print("Data ranges:", flush=True)
        datadict = resultinfo(cellsize, ignored_keys, verbose)

        statsdict = resultStats(
            cellsize,
            "insert",
            dim,
            tvol,
            fieldunits,
            ignored_keys,
            ureg,
            basedir,
            histo=ComputeHisto,
            BinCount=BinCount,
            show=show,
            verbose=verbose,
        )
        if verbose:
            print(f"insert statsdict={statsdict}", flush=True)
        stats.append(statsdict)

        # aggregate stats data
        if ComputeStats:
            createStatsTable([statsdict], "insert", fieldunits, basedir, ureg, verbose)

    return cellsize, blockdata, stats
