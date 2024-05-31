import os

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
)

from ..method import selectBlocks


def displayField(
    input,
    selectedblocks: list,
    field: str,
    fieldunits: dict,
    color,
    addruler: bool = True,
    renderView=None,
    filename: str = None,
    comment: str = None,
    polargrid: bool = False,
    printed: bool = True,
):
    """
    display field in renderview

    TODO: eventually add an annotation
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
        textDisplay.Color = [0.0, 0.0, 0.0]

    display = Show(input, renderView)

    # TODO if args.field is not Magnetic something
    display.BlockSelectors = selectedblocks
    display.ColorArrayName = color
    if polargrid:
        display.PolarAxes.Visibility = 1
        # display.PolarAxes.MaximumAngle = 360.0

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

    resolution = [1200, 1200]
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

    LUT.ApplyPreset("Rainbow Uniform", True)

    # LUT.RescaleTransferFunction(293.6058044433594, 397.88848876953125)
    # valid range for temperature but where do the range come from?
    # see post https://stackoverflow.com/questions/63028755/paraview-rescaling-colour-scheme-to-visible-data-in-range-in-python

    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    # LUT.ApplyPreset("Turbo", True)

    # Hide the scalar bar for this color map if no visible data is colored by it.
    HideScalarBarIfNotNeeded(LUT, renderView)

    # get color legend/bar for LUT in view renderView1
    LUTColorBar = GetScalarBar(LUT)
    # LUTColorBar.Position = [0.9118075801749271, 0.01059135039717564]

    keyinfo = field.split(".")
    print(f"keyinfo={keyinfo}", flush=True)
    if len(keyinfo) == 2:
        (physic, fieldname) = keyinfo
    elif len(keyinfo) == 3:
        (toolbox, physic, fieldname) = keyinfo
    else:
        raise RuntimeError(f"{field}: cannot get keyinfo as splitted char")
    symbol = fieldunits[fieldname]["Symbol"]
    msymbol = symbol
    if "mSymbol" in fieldunits[fieldname]:
        msymbol = fieldunits[fieldname]["mSymbol"]
    [in_unit, out_unit] = fieldunits[fieldname]["Units"]
    LUTColorBar.Title = rf"{msymbol} [{in_unit:~P}]"
    print(f"LUTColorBar.ComponentTitle={LUTColorBar.ComponentTitle}", flush=True)
    # LUTColorBar.ComponentTitle = ''

    # Properties modified on LUTColorBar
    LUTColorBar.TitleColor = [0.0, 0.0, 0.0]
    # LUTColorBar.TitleFontFamily = 'Courier'
    LUTColorBar.TitleBold = 1
    LUTColorBar.TitleItalic = 1
    # LUTColorBar.TitleShadow = 1
    LUTColorBar.TitleFontSize = 24
    LUTColorBar.HorizontalTitle = 1
    LUTColorBar.LabelColor = [0.0, 0.0, 0.0]
    # LUTColorBar.LabelFontFamily = 'Courier'
    # LUTColorBar.LabelBold = 1
    # LUTColorBar.LabelItalic = 1
    # LUTColorBar.LabelShadow = 1
    LUTColorBar.LabelFontSize = 20
    LUTColorBar.ScalarBarThickness = 32
    LUTColorBar.ScalarBarLength = 0.5
    LUTColorBar.AutomaticLabelFormat = 0
    # LUTColorBar.LabelFormat = '%-#6.2g'
    LUTColorBar.RangeLabelFormat = "%-#6.1f"
    LUTColorBar.DataRangeLabelFormat = "%-#5.2f"
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
            TransparentBackground=1,
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
):
    """
    create a 2D view
    """
    os.makedirs(f"{basedir}/views", exist_ok=True)
    print(f"make2Dview: field={field}", end="")
    if suffix:
        print(f", suffix={suffix}", end="")
    print(flush=True)
    print(f"blockdata={blockdata}", flush=True)

    keyinfo = field.split(".")
    print(f"keyinfo={keyinfo}", flush=True)
    if len(keyinfo) == 2:
        (physic, fieldname) = keyinfo
    elif len(keyinfo) == 3:
        (toolbox, physic, fieldname) = keyinfo
    else:
        raise RuntimeError(f"{field}: cannot get keyinfo as splitted char")
    print(f"Exclude blocks = {fieldunits[fieldname]['Exclude']}", flush=True)

    selectedblocks = selectBlocks(
        list(blockdata.keys()), fieldunits[fieldname]["Exclude"]
    )
    if selectedblocks:
        print(f"input.Selectors = {selectedblocks}", flush=True)

    filename = f"{basedir}/views/{field}.png"
    if suffix is not None:
        filename = f"{basedir}/views/{field}{suffix}.png"

    renderView = displayField(
        input,
        selectedblocks,
        field,
        fieldunits,
        color,
        addruler=addruler,
        filename=filename,
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
):
    make2Dview(
        input,
        blockdata,
        field,
        fieldunits,
        color,
        basedir,
        suffix=suffix,
        addruler=addruler,
    )
