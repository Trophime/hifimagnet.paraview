# Scripts for HiFiMagnet post-processing with Paraview

This directory contains scripts useful for post-processing feelpp results with Paraview

## `python_hifimagnetParaview.cli`:

* get range per PointData, CellData
* compute stats per PointData, CellData for insert
* compute histogram per PointData, CellData for insert
* display 3D/2D/Axi view
* display 2D OrOz view for theta in 

All data file are saved in csv format for other use.

Required
* `dimmension`: choose between 3D, 2D or Axi
* `file`: input case file (ex. Export.case)

Optional
* `--json`: 
    * give `feelpp` json file to detect exported fields
* `--views`: 
    * create views per PointData, CellData and save them to png
    * `--field`: select a field, by default get all fields
    * `--transparentBG`: enable transparent background on views
    * `--customRangeHisto`: enable custom range in views, recovered from histograms
    * `--deformedfactor`: select a deformation factor, by default 1
* `--stats`: 
    * compute stats per PointData, CellData per block (aka `feelpp` marker) 
* `--histos`: 
    * compute histogram per PointData, CellData per insert
    * `--bins`: select number of bins in histograms, by default 20
* `--plots`: 
    * create plots per PointData, CellData using given coordinates :
    * `--z`: 
    * `--theta`: 
    * `--r`: 
    * `--plotmarker`: choose marker for plots calculations
    * `--greyspace`: plot grey bar for holes (channels/slits) in plot 
    * `--show`: show plots

Optional specific to 2D:
* `--cliptheta`: select an angle to clip the geometry

Optional specific to 3D:
* `--channels`: enable creation of stl files for test-meshlib.py
* `--z`: with "--views", create a OxOy view at z
* `--theta`: with "--views", create OrOz views at theta=0/30/60/90/120/150deg

 
help

```bash
pvbatch -m python_hifimagnetParaview.cli --help
pvbatch -m python_hifimagnetParaview.cli 3D --help
pvbatch -m python_hifimagnetParaview.cli 2D --help
pvbatch -m python_hifimagnetParaview.cli Axi --help
```

examples

```bash
pvbatch -m python_hifimagnetParaview.cli 3D  ../../HL-31/test/hybride-Bh27.7T-Bb9.15T-Bs9.05T_HPfixed_BPfree/bmap/np_32/elasticity.exports/Export.case --plots --z -0.15 -0.1 -0.05 0 0.05 0.1 0.15  --r 1.94e-2 2.52e-2 3.17e-2 --save
pvbatch -m python_hifimagnetParaview.cli 2D  tmp/cfpdes-thmagel_hcurl-Axi-static-nonlinear/M9Bitters_18MW_laplace/gradH/Montgomery/Colebrook/np_16/cfpdes.exports/Export.case --views
pvbatch -m python_hifimagnetParaview.cli Axi  tmp/cfpdes-thmagel_hcurl-Axi-static-nonlinear/M9Bitters_18MW_laplace/gradH/Montgomery/Colebrook/np_16/cfpdes.exports/Export.case --stats --histos --json tmp/M9Bitters_18MW-cfpdes-thmagel_hcurl-nonlinear-Axi-sim.json
```

## running testsuite

Use paraview.simple for Python 3.10
```bash
export PYTHONPATH=/opt/paraview/lib/python3.10/site-packages/
```

To run tests located in test/
```bash
pytest
```


## `pv-statistics`:

* get range per PointData, CellData
* compute stats per PointData, CellData for insert
* compute histogram per PointData, CellData for insert
* display 3D view
* display 2D OrOz view for theta in 

All data file are saved in csv format for other use.

Optional
* `--views`: 
    * `--field`: select a field, by default get first PointData array
* `--stats`: 
    * compute stats per PointData, CellData per block (aka `feelpp` marker) 
* `--histos`: 
    * compute histogram per PointData, CellData per block (aka `feelpp` marker)
* `--plots`: 
    * `--z`: 
       * display 2D OxOy view for z in `args.z`
       * `--r`: display theta plot for r in `args.r` and z in `args.z`
    * `--save`: save plots  
 
examples

```bash
python statistics.py --help
python pv-statistics.py ../../HL-31/test/hybride-Bh27.7T-Bb9.15T-Bs9.05T_HPfixed_BPfree/bmap/np_32/thermo-electric.exports/Export.case
pvbatch pv-statistics.py ../../HL-31/test/hybride-Bh27.7T-Bb9.15T-Bs9.05T_HPfixed_BPfree/bmap/np_32/elasticity.exports/Export.case --z -0.15 -0.1 -0.05 0 0.05 0.1 0.15  --r 1.94e-2 2.52e-2 3.17e-2 --save
```

## `pv-statistics2D`:

* for 2D

examples

```bash
pvbatch pv-statistics2D.py M9Bitters_18MW_thmagel/M9Bi_18MW_elas_laplace_withoutTierod/gradH/Montgomery/Colebrook/np_1/cfpdes.exports/Export.case
```

All data file are saved in csv format for other use.

Optional
* `--views`: 
    * `--field`: select a field, by default get all PointData and CellData array
    * `--transparent`: make background transparent and fonts black
    * `--colormap`: use a given colormap 
* `--stats`: 
    * compute stats per PointData, CellData per block (aka `feelpp` marker) 
* `--histos`: 
    * compute histogram per PointData, CellData per block (aka `feelpp` marker)
* `--plots`: 
    * `--theta`: 
    * `--r`: 
    * `--save`: save plots  

## `pv-statisticsAxi`:

* for Axi

All data file are saved in csv format for other use.

Optional
* `--views`: 
    * `--field`: select a field, by default get all PointData and CellData array
    * `--transparent`: make background transparent and fonts black
    * `--colormap`: use a given colormap 
* `--stats`: 
    * compute stats per PointData, CellData per block (aka `feelpp` marker) 
* `--histos`: 
    * compute histogram per PointData, CellData per block (aka `feelpp` marker)
* `--plots`: 
    * `--z`: 
    * `--r`: 
    * `--save`: save plots  

examples

```bash
pvbatch pv-statisticsAxi.py M9Bitters_18MW_laplace/gradH/Montgomery/Colebrook/np_1/np_1/cfpdes.exports/Export.case
```


# Use with matplotlib

To view the plot:

```bash
python vonmises-vs-theta.py --file 'r=0.0194m-z=0.07m-1.csv' 'r=0.0194m-z=0.07m-0.csv' --key thermo_electric.heat.temperature --ylabel 'T [K]' --title 'Temperature in H1: r=xx, z=yy' --show
```

# Get estimation of Channel width

Need to install MeshLib:

```bash
python3 -m venv --system-site-packages meshlib-env
source ./meshlib-env/bin/activate
python3 -m pip install meshlib
```

Run:

```bash
python3 ./test-meshlib.py --help
python test-meshlib.py H*_Cu0.stl [--rfiles R[0-9]0.stl R1[0-3]0.stl] --deformed
```

To start the virtual env:

```bash
source ./meshlib-env/bin/activate
```

To exit the virtual env:

```bash
deactivate
```

References:

* [MeshLib](https://github.com/MeshInspector/MeshLib)
* see https://stackoverflow.com/questions/61159587/measure-distance-between-meshes

# To-do:
- [ ] monitor memory - for stats,histos improve by using selectblock instead of extractblock?
- [ ] create a cli with argsparse command stasts,plots,views
- [ ] split pv-statistics.py accordingly to cli commands (note add subparser option)
- [ ] Add field for Magnetostatics and ThermoElectric physics
- [ ] Add PlotOr and PlotOz
- [ ] Test with MPI
- [ ] Offset rendering - requires specific build options for paraview (see paraview downloads - healess ), '--force-offscreen-rendering'
- [ ] Add min/max to legend
- [ ] Change colormap for display
- [ ] Use eventually CustomRange for display
- [ ] Adapt for 2D
- [ ] Adapt for Axi: trouble with stats, need to rewrite this part - what about generate 3D from Axi?
- [ ] Display on selected points
- [ ] Export scene for VtkJs
- [ ] Streamline (see pv-thmagstreamline.py)
- [ ] Create Contour (see pv-contour.py)
- [ ] Support for tensor
- [ ] Create developed view of cylinder slice (see https://discourse.paraview.org/t/how-to-create-a-developed-view-of-a-cylinder-slice/14569, https://www.kitware.com/paraviews-python-programmable-filters-in-geophysics/, https://www.kitware.com/dataset-resampling-filters/)
- [ ] Compute spherical harmonics from Sphere slice ? see shtools?
- [ ] Create fieldunits, ignored_keys from json 
- [ ] Add tools to read:write:add to fieldunits json
- [ ] Create a CLI for post-processing operations by loading a json
- [ ] Run in client/server mode

- Open Questions:
  - How to discard matplotlib plots when line is not in insert?
  - Add statics per matplotlib?
   
<!-- example with pvpython and pvbatch
connect with server?

for Axi, create a 3D view and apply the rest?
ExtractSurface then RotationalExtrusion () then Transform (Rotation?) - non cree surface en 3D pour all
a tester par block
 -->

# With Python virtual env

For Linux/Mac Os X:

```bash
$ python3 -m venv --system-site-packages hifimagnet-env
$ source ./hifimagnet-env/bin/activate
$ pip install -r requirements.txt
```

# With Docker

Building docker images from python slim docker image.

* standard paraview:
  * `PYVER python` version to use (default: 3.10) 
  * `PV_VERSION_MAJOR` paraview major version (eg: 5.12)
  * `PV_VERSION_MINOR` paraview minor version (eg: 0)
  * `PV_FLAVOR`
* headless paraview:
  * `PV_FLAVOR` specify the paraview flavor to "package". Valid values are: osmesa-MPI, egl-MPI

To build the image:

```bash
docker build ...
```

# With Singularity

# References

- [feelpp](https://docs.feelpp.org/home/index.html)
- [HiFiMagnet](https://github.com/feelpp/hifimagnet)
- [Paraview](https://docs.paraview.org/en/latest/Tutorials/SelfDirectedTutorial/batchPythonScripting.html)

<!--
![Control Camera in Paraview](/assets/images/Paraview-camera.png)
--!>
