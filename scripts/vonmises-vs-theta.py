import argparse

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import MaxNLocator

from numpy import pi, arctan2

parser = argparse.ArgumentParser()
parser.add_argument(
    "--file", nargs="+", help="input txt file (ex. HL31_2018.04.13.txt)"
)
parser.add_argument(
    "--key",
    nargs="+",
    help="select key to be plotted (ex. elasticity.von_mises_criterions)",
)
parser.add_argument(
    "--legend",
    nargs="*",
    help="select key legend to be plotted (ex. von_mises_criterions)",
)

parser.add_argument("--ylabel", help="set ylabel", type=str, default="VonMises [Pa]")
parser.add_argument("--ylim", help="set ylim", type=str, nargs=2)
parser.add_argument("--yscale", help="set yscale", type=float, default="1")

parser.add_argument(
    "--title", help="set title", type=str, default="VonMises Distribution"
)
parser.add_argument(
    "--show",
    help="display graphs (requires X11 server active)",
    action="store_true",
)
parser.add_argument(
    "--save",
    help="save graphs",
    action="store_true",
)
parser.add_argument(
    "--output", help="set output filename", type=str, default="output.png"
)
args = parser.parse_args()

print(f"files: {args.file} (type: {type(args.file)}")

legends = []
if args.legend:
    legends = args.legend
    if len(legends) != len(args.file) * len(args.key):
        raise RuntimeError("there should be as many legend args as file args")
"""
else:
    for i, file in enumerate(args.file):
        for key in args.key:
            legends.append(f"{key} ({file})")
"""

ax = plt.gca()
for file in args.file:
    csv = pd.read_csv(file)  # , sep='\s+', engine='python', skiprows=1)
    keys = csv.columns.values.tolist()
    print(f'{file}: {keys}')

    for key in args.key:
        print(f"key={key}")
        keycsv = csv[["Points:0", "Points:1", "arc_length", key]]
        print(f'new keys: {keycsv.columns.values.tolist()}')
        print(f'x: {keycsv["Points:0"].describe()}')
        print(f'y: {keycsv["Points:1"].describe()}')
        
        # rename columns
        keycsv.rename(columns={"Points:0": "x", "Points:1": "y"}, inplace=True)
        print(f'new keys: {keycsv.columns.values.tolist()}')
        
        keycsv.to_csv(f"{file}-{key}.csv")
        

        # rescale columns to plot
        theta = arctan2(keycsv['y'].to_numpy(), keycsv['x'].to_numpy())  * 180 / pi
        # print(f'theta (type={type(theta)}, nelem={theta.shape}): {theta}')
        keycsv["theta"] =  keycsv.apply(
            lambda row: arctan2(row.y, row.x) * 180 / pi,
            axis=1,
        )
        # Sort the rows of dataframe by 'Name' column 
        keycsv = keycsv.sort_values(by = 'theta')
        # drop last line
        print(f'drop first line from {file}', flush=True)
        keycsv.drop(index=keycsv.index[-1],axis=0,inplace=True)
            
        # keycsv['theta'] = theta.tolist()
        keycsv[args.key] = keycsv[args.key] / args.yscale

        keycsv.plot(x="theta", y=args.key, marker='o', grid=True, ax=ax)

if args.legend:
    print('add legend')
    ax.legend(legends)
else:
    ax.get_legend().remove()

plt.xlabel('theta [deg]')
plt.ylabel(args.ylabel)

# x range
x0 = -180
x1 = 180
plt.xlim(x0, x1)

# y range
if args.ylim:
    y0 = float(args.ylim[0])
    y1 = float(args.ylim[1])
    plt.ylim(y0, y1)
    
ax.yaxis.set_major_locator(MaxNLocator(10))
plt.title(args.title)

filename = "output.png"
if args.output :
    filename = args.output
    
if args.save:
    plt.savefig(filename, dpi=300)
else:
    plt.show()
