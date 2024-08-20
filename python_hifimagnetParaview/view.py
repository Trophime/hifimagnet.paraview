import numpy as np
import pandas as pd
import re

from paraview.simple import (
    Clip,
    Slice,
    WarpByVector,
    UpdatePipeline,
)
from paraview import servermanager as sm

from .method import getbounds, invert_convert_data


def deformed(input, factor: float = 1, printed: bool = True):
    """create deformed geometry

    Args:
        input: paraview reader
        factor (float, optional): factor of deformation. Defaults to 1.
        printed (bool, optional): Defaults to True.

    Returns:
        deformed paraview reader
    """
    print(f"deformed view: factor={factor}", flush=True)

    print(f"input: PointData={list(input.PointData.keys())}", flush=True)
    print(f"input: CellData={list(input.CellData.keys())}", flush=True)

    warpByVector1 = WarpByVector(registrationName="WarpByVector1", Input=input)
    warpByVector1.Vectors = ["POINTS", "elasticity.displacement"]
    if not warpByVector1.PointData.keys():
        warpByVector1.Vectors = ["POINTS", "cfpdes.elastic.displacement"]
    warpByVector1.ScaleFactor = factor
    # get params list
    if not printed:
        for prop in warpByVector1.ListProperties():
            print(
                f"warpByVector1: {prop}={warpByVector1.GetPropertyValue(prop)}",
                flush=True,
            )

    warpByVector1.UpdatePipeline()
    print(
        f"warpByVector1: PointData={list(warpByVector1.PointData.keys())}", flush=True
    )
    print(f"warpByVector1: CellData={list(warpByVector1.CellData.keys())}", flush=True)

    sm.Fetch(warpByVector1)
    UpdatePipeline()
    return warpByVector1


def makeclip(input, name: str, invert: bool = True, printed: bool = True):
    """_summary_

    Args:
        input: paraview reader
        name (str): clip name
        invert (bool, optional): invert the clip. Defaults to True.
        printed (bool, optional): Defaults to True.

    Returns:
        clipped paraview reader
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
            print(f"clip: {prop}={clip.GetPropertyValue(prop)}", flush=True)

    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip.HyperTreeGridClipper.Origin = [0.0, 0.0, 0.0]

    clip.UpdatePipeline()
    return clip


def makethetaclip(input, theta: float, invert: bool = True, printed: bool = True):
    """create a clip with a theta angle from 2D input dataset

    Args:
        input: paraview reader
        theta (float): theta angle for the clip
        invert (bool, optional): invert the clip. Defaults to True.
        printed (bool, optional): Defaults to True.

    Returns:
        clipped paraview reader
    """
    print(f"makeclip: theta={theta}°, invert={invert}", flush=True)
    radian = theta * np.pi / 180.0
    clip = Clip(registrationName="Cliptheta", Input=input)
    clip.ClipType = "Plane"
    clip.HyperTreeGridClipper = "Plane"

    # init the 'Plane' selected for 'ClipType'
    clip.ClipType.Origin = [0.0, 0.0, 0.0]
    clip.ClipType.Normal = [-np.sin(radian), np.cos(radian), 0.0]

    print(
        f"clip.ClipType.Normal = [{-np.sin(radian)}, {np.cos(radian)}, 0.0]", flush=True
    )
    clip.Invert = 1
    if invert:
        clip.Invert = 0

    # get params list
    if not printed:
        for prop in clip.ListProperties():
            print(f"clip: {prop}={clip.GetPropertyValue(prop)}", flush=True)

    # init the 'Plane' selected for 'HyperTreeGridClipper'
    clip.HyperTreeGridClipper.Origin = [0.0, 0.0, 0.0]

    clip.UpdatePipeline()
    return clip


def makeboxclip(input, name: str):
    """create a box clip from input dataset

    Args:
        input: paraview reader
        name (str): name of the clip

    Returns:
        clipped paraview reader
    """
    print(f"makeboxclip: name={name}", flush=True)

    clip = Clip(registrationName=name, Input=input)
    clip.ClipType = "Box"

    bounds = getbounds(input)
    xmin = bounds[0]
    xmax = bounds[1]
    ymin = bounds[2]
    ymax = bounds[3]
    zmin = bounds[4]
    zmax = bounds[5]
    dz = abs(zmax - zmin)

    # Properties modified on clip1.ClipType
    clip.ClipType.Position = [xmin + xmax, ymin + ymax, zmin + zmax - dz]
    clip.ClipType.Rotation = [0.0, 0.0, 135.0]
    clip.ClipType.Length = [abs(xmax - xmin), abs(ymax - ymin), 2 * dz]
    clip.Invert = 0

    clip.UpdatePipeline()
    return clip


def makecylinderslice(input, name: str, r: float, printed: bool = True):
    """create a cylinder slice from input dataset

    Args:
        input: paraview reader
        name (str): name of the slice
        r (float): r coordinates of the slice in m
        printed (bool, optional): Defaults to True.

    Returns:
        sliced paraview reader
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
            print(f"slice: {prop}={slice1.GetPropertyValue(prop)}", flush=True)

    slice1.UpdatePipeline()
    return slice1


def makeplaneOrOzslice(input, name: str, theta: float = 0.0):
    """create a OrOz plane slice from input dataset

    Args:
        input (_type_): paraview reader
        name (str): name of the slice
        theta (float, optional): angle of normal in degrees. Defaults to 0.0.

    Returns:
        sliced paraview reader
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
    """create a plane slice from input dataset

    Args:
        input (_type_): paraview reader
        name (str): name of slice
        z (float, optional): z coordinates of the plane in m. Defaults to 0.0.

    Returns:
        sliced paraview reader
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
    """create a sphere slice from input dataset

    Args:
        input: paraview reader
        name (str): name of slice
        radius (float): radius of slice
        r (float, optional): r coordinates in meter. Defaults to 0.
        theta (float, optional): angle in degree. Defaults to 0.
        z (float, optional): z coordinates in meter. Defaults to 0.

    Returns:
        sliced paraview reader
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


def setCamera(
    renderView,
    Position: tuple = None,
    Focal: tuple = None,
    Up: tuple = None,
    Angle: float = 30,
    pProjection: bool = True,
    roll: float = 0,
    elevation: float = 0,
    azimuth: float = 0,
):
    """adapt camera settings

    ref:
    https://docs.paraview.org/en/latest/ReferenceManual/customizingParaView.html#camera-settings
    https://docs.paraview.org/en/latest/Tutorials/ClassroomTutorials/pythonAndBatchParaViewAndPython.html#control-the-camera

    Args:
        renderView: render view
        Position (tuple, optional): _description_. Defaults to None.
        Focal (tuple, optional): _description_. Defaults to None.
        Up (tuple, optional): _description_. Defaults to None.
        Angle (float, optional): _description_. Defaults to 30.
        pProjection (bool, optional): _description_. Defaults to True.
        roll (float, optional): _description_. Defaults to 0.
        elevation (float, optional): _description_. Defaults to 0.
        azimuth (float, optional): _description_. Defaults to 0.
    """

    renderView.ResetCamera()
    camera = renderView.GetActiveCamera()

    if Position is not None:
        camera.SetPosition(Position[0], Position[1], Position[2])
        camera.Roll(roll)  # rotate around ??
        # or renderView.AdjustRoll(roll)

    if Focal is not None:
        camera.SetFocalPoint(Focal[0], Focal[1], Focal[2])

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

    # print(f"help={dir(camera)}", flush=True)


def rangeHisto(field: str, fieldname: str, fieldunits: dict, filename: str) -> tuple:
    """find custom range from histogram of field
    (removes data from extremities whose area/volume is less than 0.1% of total area/volume)

    Args:
        field (str): name of field
        fieldname (str): name of field without feelpp prefixes (ex: cfpes.expr.)
        fieldunits (dict): dict of units for fields
        filename (str): name of futur view file to find corresponding histogram

    Returns:
        tuple: new custom range (None if doesn't change/ hist doesn't exist)
    """

    histfile = filename.replace("views/", "histograms/insert-").replace(
        ".png", "-histogram-matplotlib.csv"
    )
    histfile = re.sub(r"-deformed_factor\d+", "", histfile)
    histfile = re.sub(r"-OrOz-theta=\d+deg", "", histfile)
    histfile = re.sub(r"-OxOy-z=\d+.\d+mm", "", histfile)
    try:
        df = pd.read_csv(histfile)
    except:
        print(f"No histogram {histfile} - cannot use custom range !")
        return None

    for col in df.columns.tolist():
        if "Fraction of total" in col:
            grandeur = col
        else:
            fieldcol = col

    artefact_values = []
    for index, row in df.iterrows():
        fraction = row[grandeur]
        if fraction < 0.1:
            artefact_values.append(index)
        else:
            break

    for index, row in df.iloc[::-1].iterrows():
        fraction = row[grandeur]
        if fraction < 0.1:
            artefact_values.append(index)
        else:
            break

    if artefact_values:
        df.drop(index=artefact_values, inplace=True)
    else:
        return None

    units = {field: fieldunits[fieldname]["Units"]}
    [r1, r2] = invert_convert_data(
        units, [df[fieldcol].iloc[0], df[fieldcol].iloc[-1]], field
    )

    print(f"Custom range({r1:.3g},{r2:.3g})")
    return (r1, r2)
