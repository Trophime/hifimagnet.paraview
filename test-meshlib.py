"""
Compute the distance between stl
"""

import argparse
from natsort import natsorted
from natsort import natsort_keygen
import meshlib.mrmeshpy as mr

import pandas as pd
from tabulate import tabulate
import pint

parser = argparse.ArgumentParser()
parser.add_argument(
    "files", nargs="*", type=str, help="input case file (ex. H1-surf0.stl)"
)
parser.add_argument("--innerbore", type=float, help="set inner bore radius in m")
parser.add_argument("--outerbore", type=float, help="set outer bore radius in m")
parser.add_argument(
    "--deformed",
    help="compute distance for deformed",
    action="store_true",
)
args = parser.parse_args()
print(f"args: {args}")

# natsort files
files = natsorted(args.files)

shift = 1

# create an stl for inner and outer bore
if args.innerbore:
    shift = 0
    print(f"innerbore={args.innerbore}")
    # create a Stl for a cylindre with radius=innerbore
    # numpy-stl, salome, or simpler
    # add stl as 1st element of files
    # files.insert(0, istl)

if args.outerbore:
    print(f"outerbore={args.outerbore}")
    # create a Stl for a cylindre with radius=outerbore
    # numpy-stl, salome, or simpler
    # add stl as latest element of files
    # files.append(ostl)

# Original
index = []
original = []
deformed = []
expand = []
df = pd.DataFrame()

for i in range(len(files) - 1):
    mesh1 = mr.loadMesh(files[i])
    mesh2 = mr.loadMesh(files[i + 1])
    z = mr.findSignedDistance(mesh1, mesh2)
    # print(f"Channel{i+shift}: {z.signedDist} m", flush=True)  # around 0.8 mm
    index.append(f"Channel{i+shift}")
    original.append(z.signedDist / 1.0e-3)

# Deformed
if args.deformed:
    dfiles = []
    for i, file in enumerate(files):
        (filename, ext) = file.split("0.")
        # print(f"filename={filename}, ext={ext}")
        new_filename = f"{filename}-deformed0.{ext}"
        dfiles.append(new_filename)

    for i in range(len(dfiles) - 1):
        mesh1 = mr.loadMesh(dfiles[i])
        mesh2 = mr.loadMesh(dfiles[i + 1])
        z = mr.findSignedDistance(mesh1, mesh2)
        # print(f"Channel{i+shift}: {z.signedDist} m", flush=True)  # around 0.8 mm
        deformed.append(z.signedDist / 1.0e-3)
        expand.append((deformed[-1] / original[i] - 1) * 100.0)

    df = pd.DataFrame.from_dict(
        {
            "Name": index,
            "Original [mm]": original,
            "Deformed [mm]": deformed,
            "Expand [%]": expand,
        }
    )

else:
    df = pd.DataFrame.from_dict({"Name": index, "Original [mm]": original})

df.sort_values(by="Name", key=natsort_keygen())
print(tabulate(df, headers="keys", tablefmt="psql", showindex=False))
