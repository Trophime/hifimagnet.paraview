# Scripts for 'HiFiMagnet' post-processing with 'Paraview'

This directory contains scripts useful for post-processing 'feelpp' results with 'Paraview'

pv-statistics:

* get range per PointData, CellData
* compute stats per PointData, CellData for insert
* compute histogram per PointData, CellData for insert
* display 3D view
* display 2D OrOz view for theta in 

All data file are saved in csv format for other use.

Optional
* select a field, by default get first PointData array
* compute stats per PointData, CellData per block (aka feelpp marker) 
* compute histogram per PointData, CellData per block (aka feelpp marker)
* display 2D OxOy view for z in args.z
* display theta plot for r in args.r and z in args.z
 
examples

```bash
python statistics.py --help
python pv-statistics.py ../../HL-31/test/hybride-Bh27.7T-Bb9.15T-Bs9.05T_HPfixed_BPfree/bmap/np_32/thermo-electric.exports/Export.case
pvbatch pv-statistics.py ../../HL-31/test/hybride-Bh27.7T-Bb9.15T-Bs9.05T_HPfixed_BPfree/bmap/np_32/elasticity.exports/Export.case --z -0.15 -0.1 -0.05 0 0.05 0.1 0.15  --r 1.94e-2 2.52e-2 3.17e-2 --save
```

# Use with 'matplotlib'

To view the plot:

```bash
python vonmises-vs-theta.py --file 'r=0.0194m-z=0.07m-1.csv' 'r=0.0194m-z=0.07m-0.csv' --key thermo_electric.heat.temperature --ylabel 'T [K]' --title 'Temperature in H1: r=xx, z=yy' --show
```

# Get estimation of Channel width

Need to install 'MeshLib':

```bash
python3 -m venv --system-site-packages meshlib-env
source ./meshlib-env/bin/activate
python3 -m pip install meshlib
```

Run:

```bash
python3 ./test-meshlib.py --help
python test-meshlib.py H*_Cu0.stl --deformed
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
- [ ] Add min/max to legend
- [ ] Change colormap for display
- [ ] Use eventually CustomRange for display
- [ ] Display on selected points
- [ ] Export scene for VtkJs
- [ ] Streamline (see pv-thmagstreamline.py)
- [ ] Create Contour (see pv-contour.py)
- [ ] Support for tensor
- [ ] Adapt for 2D
- [ ] Adapt for Axi: trouble with stats, need to rewrite this part - what about generate 3D from Axi?
- [ ] Test with MPI
- [ ] Offset rendering - requires specific build options for paraview
- [ ] Create developed view of cylinder slice (see https://discourse.paraview.org/t/how-to-create-a-developed-view-of-a-cylinder-slice/14569, https://www.kitware.com/paraviews-python-programmable-filters-in-geophysics/, https://www.kitware.com/dataset-resampling-filters/)
- [ ] Compute spherical harmonics from Sphere slice ?
- [ ] Create fieldunits, ignored_keys from json 
- [ ] Add tools to read:write:add to fieldunits json
- [ ] Create a CLI for post-processing operations by loading a json
- [ ] Run in client/server mode


<!-- example with pvpython and pvbatch
pvbatch '--force-offscreen-rendering'
test with mpi
connect with server?

for Axi, create a 3D view and apply the rest?
ExtractSurface then RotationalExtrusion () then Transform (Rotation?) - non cree surface en 3D pour all
a tester par block
 -->
# References

- [feelpp](https://docs.feelpp.org/home/index.html)
- [HiFiMagnet](https://github.com/feelpp/hifimagnet)
- [Paraview](https://docs.paraview.org/en/latest/Tutorials/SelfDirectedTutorial/batchPythonScripting.html)