"""
Show how to create distribution in Axi

Paraview:
Convert PointData to CellData
Cellsize
CellCenters with Vertex Cells option
create SpreadSheatView
export data to csv

see pv-statistics.py to automate this
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import hist
from math import pi

csv = pd.read_csv('Axi-B-cellcenters.csv')
keys = csv.columns.values.tolist()
print(f"read_csv: csv=Axi-B-cellcenters.csv, keys={keys}", flush=True)
    
csv['AxiSurf'] = csv['Points_0'] * csv['Area'] * (2*pi) / 1.e-9 
sum = csv["AxiSurf"].sum()
csv['AxiVol'] = csv['AxiSurf'] / sum
print(f'Sum(AxiSurf)={sum}')
print(f'Sum(AxiVol)={csv["AxiVol"].sum()}')
print(f'Stats(Bz)={csv["cfpdes.expr.B_1"].describe()}')

# from Bounding Box
r = 3
rG = r/2.
h = 2 * 1.81616
V = 2*pi*rG*r*h/1.e-9
print(f'V={V} mm³')
print(f'V-Sum()={abs(1-sum/V)} %')

# histo
hist(csv["cfpdes.expr.B_1"], bins=40, weights=csv['AxiVol'])

plt.xlabel('Bz [T]')
plt.ylabel('[%]')
plt.legend('Bz')
plt('Bz Distribution')
plt.grid()
plt.show()
