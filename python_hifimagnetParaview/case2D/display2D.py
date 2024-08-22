import os
import gc

from paraview.simple import (
    BoundingRuler,
    Text,
    Delete,
    CreateView,
    Show,
    GetScalarBar,
    GetColorTransferFunction,
    ColorBy,
    HideScalarBarIfNotNeeded,
    GetOpacityTransferFunction,
    SaveScreenshot,
    ExtractBlock,
)

from ..method import selectBlocks, keyinfo
from ..view import rangeHisto


def displayField(
    input,
    selectedblocks: list[str],
    field: str,
    fieldunits: dict,
    color,
    addruler: bool = True,
    renderView=None,
    filename: str = None,
    comment: str = None,
    polargrid: bool = False,
    printed: bool = True,
    excludeBlocks: bool = False,
    background: bool = False,
    customRangeHisto: bool = False,
):
    """display field in renderview

    Args:
        input: paraview reader
        selectedblocks (list[str]): list of markers of field
        field (str): field name
        fieldunits (dict): dict of field units
        color (_type_): color PointData or CellData
        addruler (bool, optional): add ruler to view. Defaults to True.
        renderView (optional): pre-existing renderview. Defaults to None.
        filename (str, optional): name and path of futur view file. Defaults to None.
        comment (str, optional): add comment. Defaults to None.
        polargrid (bool, optional): add polar grid to view. Defaults to False.
        printed (bool, optional): _description_. Defaults to True.
        excludeBlocks (bool, optional): field excluded blocks. Defaults to False.
        background (bool, optional): transparent background (& text black). Defaults to False.
        customRangeHisto (bool, optional): create custom range from field histogram. Defaults to False.

    Returns:
        renderView
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
        if background:
            textDisplay.Color = [0.0, 0.0, 0.0]

    display = Show(input, renderView)

    # TODO if args.field is not Magnetic something
    display.BlockSelectors = selectedblocks
    display.ColorArrayName = color
    if polargrid:
        display.PolarAxes.Visibility = 1
        # display.PolarAxes.MaximumAngle = 360.0
        if background:
            display.PolarAxes.PolarAxisColor = [0.0, 0.0, 0.0]
            display.PolarAxes.PolarArcsColor = [0.0, 0.0, 0.0]
            display.PolarAxes.LastRadialAxisColor = [0.0, 0.0, 0.0]
            display.PolarAxes.SecondaryPolarArcsColor = [0.0, 0.0, 0.0]
            display.PolarAxes.SecondaryRadialAxesColor = [0.0, 0.0, 0.0]
            display.PolarAxes.PolarAxisTitleColor = [0.0, 0.0, 0.0]
            display.PolarAxes.PolarAxisLabelColor = [0.0, 0.0, 0.0]
            display.PolarAxes.LastRadialAxisTextColor = [0.0, 0.0, 0.0]
            display.PolarAxes.SecondaryRadialAxesTextColor = [0.0, 0.0, 0.0]
        display.PolarAxes.Use2DMode = 0

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
                print(f"ruler: {prop}={ruler.GetPropertyValue(prop)}", flush=True)
        Show(ruler, renderView)  # Reset Camera

    resolution = [1600, 1200]
    renderView.ViewSize = resolution
    if not printed:
        for prop in renderView.ListProperties():
            print(f"renderView: {prop}={renderView.GetPropertyValue(prop)}", flush=True)
    renderView.OrientationAxesVisibility = 1
    display.SetScalarBarVisibility(renderView, True)
    display.RescaleTransferFunctionToDataRange(True, False)
    # get color transfer function/color map for 'thermo_electricheattemperature'
    field_name = field.replace(".", "")
    LUT = GetColorTransferFunction(field_name)
    LUT.ScalarRangeInitialized = 1.0

    if input.GetDataInformation().DataInformation.GetNumberOfUniqueBlockTypes() == 0:
        excludeBlocks = False
        print("only one block: excludeBlocks= False", flush=True)
    if excludeBlocks:
        print(f"Extract block: {selectedblocks}", flush=True)
        extractBlock1 = ExtractBlock(registrationName="insert", Input=input)
        extractBlock1.Selectors = selectedblocks
        extractBlock1.UpdatePipeline()
        if field in list(extractBlock1.CellData.keys()):
            arrayInfo = extractBlock1.CellData[field]
        if field in list(extractBlock1.PointData.keys()):
            arrayInfo = extractBlock1.PointData[field]
        r = arrayInfo.GetRange(0)
        LUT.RescaleTransferFunction(r[0], r[1])
        Delete(extractBlock1)
        del extractBlock1
        collected = gc.collect()

    # valid range for temperature but where do the range come from?
    # see post https://stackoverflow.com/questions/63028755/paraview-rescaling-colour-scheme-to-visible-data-in-range-in-python

    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    # LUT.ApplyPreset("Turbo", True)
    LUT.ApplyPreset("Rainbow Uniform", True)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(LUT, renderView)

    # get color legend/bar for LUT in view renderView1
    LUTColorBar = GetScalarBar(LUT)
    # LUTColorBar.Position = [0.9118075801749271, 0.01059135039717564]

    (toolbox, physic, fieldname) = keyinfo(field)

    symbol = fieldunits[fieldname]["Symbol"]
    msymbol = symbol
    if "mSymbol" in fieldunits[fieldname]:
        msymbol = fieldunits[fieldname]["mSymbol"]
    [in_unit, out_unit] = fieldunits[fieldname]["Units"]
    LUTColorBar.Title = rf"{msymbol} [{in_unit:~P}]"
    print(f"LUTColorBar.ComponentTitle={LUTColorBar.ComponentTitle}", flush=True)
    # LUTColorBar.ComponentTitle = ''

    if customRangeHisto:
        r = rangeHisto(field, fieldname, fieldunits, filename)
        if r:
            LUT.RescaleTransferFunction(r[0], r[1])

    # Properties modified on LUTColorBar
    if background:
        LUTColorBar.TitleColor = [0.0, 0.0, 0.0]
        LUTColorBar.LabelColor = [0.0, 0.0, 0.0]
    # LUTColorBar.TitleFontFamily = 'Courier'
    LUTColorBar.TitleBold = 1
    LUTColorBar.TitleItalic = 1
    # LUTColorBar.TitleShadow = 1
    LUTColorBar.TitleFontSize = 30
    LUTColorBar.HorizontalTitle = 1
    # LUTColorBar.LabelFontFamily = 'Courier'
    # LUTColorBar.LabelBold = 1
    # LUTColorBar.LabelItalic = 1
    # LUTColorBar.LabelShadow = 1
    LUTColorBar.LabelFontSize = 30
    LUTColorBar.ScalarBarThickness = 32
    LUTColorBar.ScalarBarLength = 0.9
    LUTColorBar.AutomaticLabelFormat = 0
    # LUTColorBar.LabelFormat = '%-#6.2g'
    LUTColorBar.RangeLabelFormat = "%-#6.3g"
    LUTColorBar.DataRangeLabelFormat = "%-#6.3e"
    LUTColorBar.DrawDataRange = 1
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
            TransparentBackground=background,
        )

    if comment is not None:
        Delete(textDisplay)
    Delete(display)
    return renderView


################################################################
# create a 2D view
def make2Dview(
    input,
    blockdata,
    field: str,
    fieldunits: dict,
    color,
    basedir: str,
    suffix: str = None,
    addruler: bool = False,
    printed: bool = True,
    background: bool = False,
    customRangeHisto: bool = False,
):
    """create a 2D view

    Args:
        input: paraview reader
        blockdata: blockdata from meshinfo
        field (str): field name
        fieldunits (dict): dict of field units
        color: color for PointData or CellData
        basedir (str): result directory
        suffix (str, optional): None or -deformed. Defaults to None.
        addruler (bool, optional): add ruler to view. Defaults to False.
        printed (bool, optional): Defaults to True.
        background (bool, optional): transparent background (& text black). Defaults to False.
        customRangeHisto (bool, optional):  create custom range from field histogram. Defaults to False.
    """
    os.makedirs(f"{basedir}/views", exist_ok=True)
    print(f"make2Dview: field={field}", end="")
    if suffix:
        print(f", suffix={suffix}", end="")
    print(flush=True)

    (toolbox, physic, fieldname) = keyinfo(field)

    excludeBlocks = False
    if fieldunits[fieldname]["Exclude"]:
        excludeBlocks = True
        print(f"Exclude blocks = {fieldunits[fieldname]['Exclude']}", flush=True)

    selectedblocks = selectBlocks(
        list(blockdata.keys()), fieldunits[fieldname]["Exclude"]
    )
    if selectedblocks:
        print(f"input.Selectors = {selectedblocks}", flush=True)

    filename = f"{basedir}/views/{field}.png"
    if suffix is not None:
        filename = f"{basedir}/views/{field}{suffix}.png"

    comm = ""
    if fieldunits["Current"]["Val"]:
        comm = f'I={fieldunits["Current"]["Val"]}'
    if fieldunits["B0"]["Val"]:
        comm = comm + f'\nB0={fieldunits["B0"]["Val"]}T'
    if fieldunits["Bbg"]["Val"]:
        comm = comm + f'\nBackground field: {fieldunits["Bbg"]["Val"]}'
    renderView = displayField(
        input,
        selectedblocks,
        field,
        fieldunits,
        color,
        addruler=addruler,
        filename=filename,
        comment=comm,
        excludeBlocks=excludeBlocks,
        background=background,
        customRangeHisto=customRangeHisto,
    )

    Delete(renderView)
    del renderView


def makeview(
    args,
    input,
    blockdata,
    field: str,
    fieldunits: dict,
    color,
    basedir: str,
    suffix: str = None,
    addruler: bool = False,
    printed: bool = True,
    background: bool = False,
    customRangeHisto: bool = False,
):
    """create views

    Args:
        args: options
        input: paraview reader
        blockdata: blockdata from meshinfo
        field (str): field name
        fieldunits (dict): dict of field units
        color: color for PointData or CellData
        basedir (str):  result directory
        suffix (str, optional):  None or -deformed. Defaults to None.
        addruler (bool, optional): add ruler to view. Defaults to False.
        printed (bool, optional): Defaults to True.
        background (bool, optional): transparent background (& text black). Defaults to False.
        customRangeHisto (bool, optional):  create custom range from field histogram. Defaults to False.
    """

    make2Dview(
        input,
        blockdata,
        field,
        fieldunits,
        color,
        basedir,
        suffix=suffix,
        addruler=addruler,
        background=background,
        customRangeHisto=customRangeHisto,
    )
