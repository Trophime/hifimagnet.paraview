import pytest
import gmsh
import re
import json
import numpy as np
import pandas as pd

from PIL import Image

from python_hifimagnetParaview.cli import init
from python_hifimagnetParaview.json import returnExportFields
from python_hifimagnetParaview.method import convert_data
from python_hifimagnetParaview.meshinfoAxi import meshinfo
from python_hifimagnetParaview.case2D.display2D import makeview
from python_hifimagnetParaview.caseAxi.methodAxi import create_dicts_fromjson


axi_cases = [
    (
        "./test/cases/cfpdes-thmagel-Axi-static-linear/Tore/np_1/cfpdes.exports/Export.case",
        "./test/models/cfpdes_Tore_axi_static/Toretest-cfpdes-thmagel_hcurl-Axi-sim.json",
    )
]
dim = 2
axis = True


@pytest.mark.parametrize("file,jsonfile", axi_cases)
def test_init(file, jsonfile):

    meshpath = file.replace("cfpdes.exports/Export.case", "cfpdes.mesh.path")
    with open(meshpath, "r") as f:
        meshfile = f.readlines()[0]

    meshfile = "./" + meshfile.split("hifimagnet.paraview/")[1]
    meshfile = re.sub(r"_p\d+.json", ".msh", meshfile)

    gmsh.initialize()
    gmsh.open(meshfile)

    # Paraview : sum of nodes for each marker
    #   -> doesn't account for duplicate nodes on boundary between markers,
    # so the same is done for gmsh in this test = sum(NodesForPhysicalGroup)
    groups = gmsh.model.getPhysicalGroups()
    NumberOfNodes = 0
    for g in groups:
        if g[0] == 2:
            nodestag, _ = gmsh.model.mesh.getNodesForPhysicalGroup(g[0], g[1])
            NumberOfNodes += len(nodestag)

    _, elementTags, _ = gmsh.model.mesh.getElements(dim=2)
    NumberOfElements = len(elementTags[0])
    (cwd, basedir, ureg, distance_unit, reader) = init(file)

    dataInfo = reader.GetDataInformation()
    assert (
        dataInfo.GetNumberOfPoints() == NumberOfNodes
    ), f"Number of Points: paraview:{dataInfo.GetNumberOfPoints()} != msh:{NumberOfNodes}"
    assert (
        dataInfo.GetNumberOfCells() == NumberOfElements
    ), f"Number of Cells: paraview:{dataInfo.GetNumberOfCells()} != msh:{NumberOfElements}"

    fieldtype = returnExportFields(jsonfile, basedir)
    fieldunits, ignored_keys = create_dicts_fromjson(
        fieldtype, ureg, distance_unit, basedir
    )

    with open(f"{basedir}/FieldType.json", "r") as jsonfile:
        fieldtypejson = json.load(jsonfile)

    assert (
        fieldtype.keys() == fieldtypejson.keys()
    ), f"fieldtype: {fieldtype.keys()} != json:{fieldtypejson.keys()}"

    with open(f"{basedir}/fieldunits.json", "r") as jsonfile:
        fieldunitsjson = json.load(jsonfile)

    assert (
        fieldunits.keys() == fieldunitsjson.keys()
    ), f"fieldunits: {fieldunits.keys()} != json:{fieldunitsjson.keys()}"

    assert len(fieldtype.keys()) <= len(
        fieldunits.keys()
    ), f"fieldtype:{len(fieldtype.keys())} > fieldunits:{len(fieldunits.keys())}"


def assert_images_equal(image_1: str, image_2: str):
    img1 = Image.open(image_1)
    img2 = Image.open(image_2)

    # Convert to same mode and size for comparison
    img2 = img2.convert(img1.mode)
    img2 = img2.resize(img1.size)

    sum_sq_diff = np.sum(
        (np.asarray(img1).astype("float") - np.asarray(img2).astype("float")) ** 2
    )

    if sum_sq_diff == 0:
        # Images are exactly the same
        pass
    else:
        normalized_sum_sq_diff = sum_sq_diff / np.sqrt(sum_sq_diff)
        assert (
            normalized_sum_sq_diff < 0.001
        ), f'{image_1.split("/")[-1]}: {normalized_sum_sq_diff} > 0.001'


@pytest.mark.parametrize("file,jsonfile", axi_cases)
def test_views(file, jsonfile):

    (cwd, basedir, ureg, distance_unit, reader) = init(file)

    fieldtype = returnExportFields(jsonfile, basedir)
    fieldunits, ignored_keys = create_dicts_fromjson(
        fieldtype, ureg, distance_unit, basedir
    )
    cellsize, blockdata, statsdict = meshinfo(
        reader, dim, fieldunits, ignored_keys, basedir, ureg, ComputeStats=False
    )

    for field in [
        "cfpdes.heat.temperature",
        "cfpdes.elastic.displacement",
        # "cfpdes.expr.B",
    ]:
        if field in list(cellsize.CellData.keys()):
            color = ["CELLS", field]
        if field in list(cellsize.PointData.keys()):
            color = ["POINTS", field]
        if field in list(cellsize.PointData.keys()) + list(cellsize.CellData.keys()):
            makeview(
                None,
                cellsize,
                blockdata,
                field,
                fieldunits,
                color,
                basedir,
                suffix="",
                addruler=False,
                background=False,
            )

        imageref = f"./test/Pictures/Axi/{field}.png"
        imagenew = f"{basedir}/views/{field}.png"
        assert_images_equal(imageref, imagenew)


### pour gros fichiers git lfs


@pytest.mark.parametrize("file,jsonfile", axi_cases)
def test_stats(file, jsonfile):

    (cwd, basedir, ureg, distance_unit, reader) = init(file)

    fieldtype = returnExportFields(jsonfile, basedir)
    fieldunits, ignored_keys = create_dicts_fromjson(
        fieldtype, ureg, distance_unit, basedir
    )
    cellsize, blockdata, statsdict = meshinfo(
        reader, dim, fieldunits, ignored_keys, basedir, ureg, ComputeStats=True
    )

    statsheat = pd.read_csv(
        f"{basedir}/stats/cfpdes.heat.temperature-descriptivestats-create.csv"
    )
    try:
        heatmeasures = pd.read_csv(
            f'{basedir.replace("cfpdes.exports/paraview.exports","heat.measures/values.csv")}'
        )
    except:
        heatmeasures = None

    if isinstance(heatmeasures, pd.DataFrame):
        units = {"temperature": fieldunits["temperature"]["Units"]}

        Feel_T_max = convert_data(
            units, heatmeasures["Statistics_Stat_T_max"].iloc[0], "temperature"
        )
        Feel_T_mean = convert_data(
            units, heatmeasures["Statistics_Stat_T_mean"].iloc[0], "temperature"
        )
        Feel_T_min = convert_data(
            units, heatmeasures["Statistics_Stat_T_min"].iloc[0], "temperature"
        )

        assert (
            abs(1 - Feel_T_max / statsheat["Maximum"].iloc[0]) < 0.001
        ), f'Tmax: abs(1-Feel:{Feel_T_max}/Paraview:{statsheat["Maximum"].iloc[0]}) > 0.001'
        assert (
            abs(1 - Feel_T_mean / statsheat["Mean"].iloc[0]) < 0.001
        ), f'Tmean: abs(1-Feel:{Feel_T_mean}/Paraview:{statsheat["Mean"].iloc[0]}) > 0.001'
        assert (
            abs(1 - Feel_T_min / statsheat["Minimum"].iloc[0]) < 0.001
        ), f'Tmin: abs(1-Feel:{Feel_T_min}/Paraview:{statsheat["Minimum"].iloc[0]}) > 0.001'

    statselastic = pd.read_csv(
        f"{basedir}/stats/cfpdes.expr.Vonmises-descriptivestats-create.csv"
    )
    try:
        elasticmeasures = pd.read_csv(
            f'{basedir.replace("cfpdes.exports/paraview.exports","elastic.measures/values.csv")}'
        )
    except:
        elasticmeasures = None

    if isinstance(elasticmeasures, pd.DataFrame):
        units = {"Vonmises": fieldunits["Vonmises"]["Units"]}

        Feel_VM_max = convert_data(
            units, elasticmeasures["Statistics_VonMises_Tore_max"].iloc[0], "Vonmises"
        )
        Feel_VM_mean = convert_data(
            units, elasticmeasures["Statistics_VonMises_Tore_mean"].iloc[0], "Vonmises"
        )
        Feel_VM_min = convert_data(
            units, elasticmeasures["Statistics_VonMises_Tore_min"].iloc[0], "Vonmises"
        )

        assert (
            abs(1 - Feel_VM_max / statselastic["Maximum"].iloc[0]) < 0.01
        ), f'VonMisesmax: abs(1-Feel:{Feel_VM_max}/Paraview:{statselastic["Maximum"].iloc[0]}) > 0.01'
        assert (
            abs(1 - Feel_VM_mean / statselastic["Mean"].iloc[0]) < 0.01
        ), f'VonMisesmean: abs(1-Feel:{Feel_VM_mean}/Paraview:{statselastic["Mean"].iloc[0]}) > 0.01'
        assert (
            abs(1 - Feel_VM_min / statselastic["Minimum"].iloc[0]) < 0.01
        ), f'VonMisesmin: abs(1-Feel:{Feel_VM_min}/Paraview:{statselastic["Minimum"].iloc[0]}) > 0.01'
