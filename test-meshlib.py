"""
Compute the distance between stl
"""

import argparse
from natsort import natsorted
from natsort import natsort_keygen
import meshlib.mrmeshpy as mr

import copy
import pandas as pd
from tabulate import tabulate
import pint

parser = argparse.ArgumentParser()
parser.add_argument(
    "files", nargs="*", type=str, help="input case file (ex. H1-surf0.stl)"
)
parser.add_argument(
    "--rfiles", nargs="*", type=str, help="input case file (ex. R1-surf0.stl)"
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
rfiles = []
if args.rfiles:
    if len(args.rfiles) != len(args.files) - 1:
        raise RuntimeError(
            f"Expect to have {len(args.files)-1} rfiles: got {len(args.rfiles)}"
        )
    rfiles = natsorted(args.rfiles)
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
df = pd.DataFrame()

for i in range(len(files) - 1):
    mesh1 = mr.loadMesh(files[i])
    mesh2 = mr.loadMesh(files[i + 1])
    z = mr.findSignedDistance(mesh1, mesh2)
    # print(f"Channel{i+shift}: {z.signedDist} m", flush=True)  # around 0.8 mm
    index.append(f"Channel{i+shift}")
    original.append(z.signedDist / 1.0e-3)
# print(f"original={original} len={len(original)}")

if rfiles:
    roriginal = copy.deepcopy(original)
    # TODO R1/innerbore
    roriginal[0] = 0
    for i in range(1, len(files) - 2):
        mesh1 = mr.loadMesh(rfiles[i - 1])
        mesh2 = mr.loadMesh(rfiles[i + 1])
        z = mr.findSignedDistance(mesh1, mesh2)
        # print(
        #     f"Channel{i+shift} R{i-1}/R{i+1}: {z.signedDist} m", flush=True
        # )  # around 0.8 mm
        roriginal[i] = z.signedDist / 1.0e-3
    # TODO R13/outerbore
    roriginal[len(files) - 2] = 0
    # print(f"roriginal={roriginal} len={len(roriginal)}")

# Deformed
if args.deformed:
    deformed = []
    expand = []
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

    if rfiles:
        rdeformed = copy.deepcopy(deformed)
        rexpand = copy.deepcopy(expand)
        dfiles = []
        for i, file in enumerate(rfiles):
            (filename, ext) = file.split("0.")
            # print(f"filename={filename}, ext={ext}")
            new_filename = f"{filename}-deformed0.{ext}"
            dfiles.append(new_filename)

        # TODO R1/innerbore
        rdeformed[0] = 0
        rexpand[0] = 0
        for i in range(1, len(dfiles) - 2):
            mesh1 = mr.loadMesh(dfiles[i - 1])
            mesh2 = mr.loadMesh(dfiles[i + 1])
            z = mr.findSignedDistance(mesh1, mesh2)
            # print(
            #     f"deformed Channel{i+shift} R{i-1}/R{i+1}: {z.signedDist} m", flush=True
            # )  # around 0.8 mm
            rdeformed[i] = z.signedDist / 1.0e-3
            rexpand[i] = (rdeformed[i] / roriginal[i] - 1) * 100.0
        # TODO R13/outerbore
        rdeformed[len(files) - 2] = 0
        rexpand[len(files) - 2] = 0

    df = pd.DataFrame.from_dict(
        {
            "Name": index,
            "Original [mm]": [f"{val:.3f}" for val in original],
            "Deformed [mm]": [f"{val:.3f}" for val in deformed],
            "Expand [%]": [f"{val:.3f}" for val in expand],
        }
    )

    if rfiles:
        df["rOriginal [mm]"] = [f"{val:.3f}" for val in roriginal]
        df["rDeformed [mm]"] = [f"{val:.3f}" for val in rdeformed]
        df["rExpand [%]"] = [f"{val:.3f}" for val in rexpand]


else:
    df = pd.DataFrame.from_dict(
        {"Name": index, "Original [mm]": [f"{val:.3f}" for val in original]}
    )

df.sort_values(by="Name", key=natsort_keygen())
print(tabulate(df, headers="keys", tablefmt="psql", showindex=False))
df.to_csv("Channels-meshlib.csv")
