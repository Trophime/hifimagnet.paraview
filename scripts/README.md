# Additionnal indivdual scripts

This folder regroups other python scripts for other display uses.

## profiles.py

Create plots from feelpp bitters Axi and 2D results with paraview for bitters to compare them with Ansys or Simple Model results 
=> is uses in conjonction with `display_results_v0.1.py`

### Options
dir: add the directory where you want the results (must contain subdirectory like cooling/heatcorrelation/friction/np_...)
--np: specify nb of cores for results
--coolings: specify all coolings you want to display
--heatcorrelations: specify all heat correlations you want to display
--frictions: specify all frictions you want to display
--rint: specify inner radius (if 2D give 1, if axi give 2, one for each bitter)
--rint: specify outer radius (if 2D give 1, if axi give 2, one for each bitter)
--toolbox: specify toolbox (cfpdes or thermo_electric)

axi: profiles for axi coord
    --measures: give measures to add to plot (can be multiple)
                (-for "display_results_v0.1.py plot_profiles", need "Jth","temperature"-)
                (-can also add B for profile on Z axis-)
    --Z: give one Z for each Bitter for 2ndtop & bottom profiles
    --Zmax: give Zmax for B profile on Z axis

2D: profiles for 2D coord
    --measures: give measures to add to plot (can be multiple)
                (-for "display_results_v0.1.py plot_profiles", need "temperature"-)
    --yaml_file: input yaml file of bitter for more precise boxplots
    --nb: specify number of boxplots (only used if yaml_file is not specified)
    --theta: specify angle of the slice (default pi/16)

### Examples
```bash
python profiles.py axi bitters_old/M9Bitters-frans --rint 0.2 0.343 --rext 0.34 0.5 --Z 0.22292 0.250615 --Zmax 0.15 --coolings gradH meanH --frictions Colebrook --measures Jth temperature B
python profiles.py 2D bitters_old/M9Bi-frans-nonlinear-ansys/ --yaml_file yaml/M9_Bi.yaml --rint 0.2 --rext 0.34 --coolings gradH meanH --frictions Colebrook
```


## profiles.py

Compare feelpp results with pupitre or Ansys or SimpleModel results
    * plot_measure: plot measures.csv after workflow [vs pupitre results]
    * plot_values: plot *.measures/values.csv after workflow [vs Simple Model results]
    * plot_profiles: plot profiles created by `profiles.py`

### Options:
  --mtype: bitter or helix
  --filter: magnet filter if measures from MSite
  --dir: add the directory where you want the results (must contain subdirectory like cooling/heatcorrelation/friction/np_...)
  --np: specify nb of cores for results
  --pupitrefile: txt file of exp
  --coolings: specify all coolings you want to display
  --heatcorrelations: specify all heat correlations you want to display
  --frictions: specify all frictions you want to display
  --oneI: IF RESULTS ARE NOT COMMISSIONING, ONLY ONE I
  --save: save plots & err tables
  --show: show plots
plot_measure: plot measures from measures.csv vs I & give relative error
  --measures: choose from ["Ucoil","B0","Power","flow","R(I)","MinTH","MaxTH","MeanTH","TH","Tout"] (can be multiple)
plot_values
  --measures: choose from ["Uw","DT","Flux","HeatCoeff","PowerH","statsTH","cf"] (can be multiple)
  --current: set specific currents to plot out of the hole commissioning
  --modeldir: specify external model directory with : [..].measures/values.csv
plot_profiles:
  files: specify simple model files containing temperature profiles
  --measures: choose from ["T","J"] (can be multiple)
  --parts: choose on which parts of the magnet the profile will be plot; choices=["2ndtop", "core", "2ndbot"] (can be multiple)
  --toolbox: specify toolbox used between ["cfpdes", "thermo_electric"]
  --hideSimpleModel: hide the SimpleModel profiles
  --boxplotdir: specify 2D measures directory
  --show2Dlines: show 2D profiles
  --Ansys: specify Ansys file containing temperature profiles

### Examples
```bash
python display_results_v0.1.py plot_measures --mtype bitter --dir bitters_new/M9Bitters_18MWcommi9june/ --coolings meanH gradH gradHZ --frictions Colebrook  --pupitrefile pupitre/M9_2023.06.09---13\:48\:53_from13\:54\:09to14\:00\:05.txt --save --show --measures "Ucoil" "B0" "Power" "flow" "R(I)" "MinTH" "MaxTH" "MeanTH" "TH" Tout
python display_results_v0.1.py plot_values --mtype bitter --dir bitters_new/M9Bitters_18MW25oct/ --coolings meanH gradH gradHZ --frictions Colebrook --oneI --pupitrefile pupitre/M9_2023.10.25---15\:50\:45_to16\:08\:42.txt --save --show --measures "Uw" "DT" "Flux" "HeatCoeff" "PowerH" "statsTH" "cf" --modeldir Exported_Results_Modified_Inner_Bitter/Exported_Results/
python display_results_v0.1.py plot_profiles '{"inner":"Exported_Results/20231127-LNCMI_Inner_Original-Simple_Model_Data.csv","outer":"Exported_Results/20231127-LNCMI_Outer_Original-Simple_Model_Data.csv"}' --coolings mean grad --frictions Constant --dir bitters_old/M9Bitters-frans --save --show --measures T J
python display_results_v0.1.py plot_profiles '{"inner":"Exported_Results/20231127-LNCMI_Inner_Original-Simple_Model_Data.csv"}' --part core --coolings meanH --frictions Colebrook --dir bitters_old/M9Bitters-frans --boxplotsdir bitters_old/M9Bi-frans-nonlinear-ansys/ --save --show --show2Dlines --Ansys Exported_Results_fix/20240129-LNCMI_Inner_Original-Cross-sectional_data.csv
```

## test-tin.py

View R(I) as function of Tin

### Options:
input_files: enter input file(s) (ex: ~/R_14Helices_Tinit=20degCelsius.csv)
--ikey: specify I key (starting from 1)
--list: list valid ikeys values (action="store_true")