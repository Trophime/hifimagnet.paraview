import argparse, argcomplete
import os
import sys
import json
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image

import warnings

warnings.filterwarnings("ignore")


def options(description: str, epilog: str):
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument("--mdata", help="specify results data", type=json.loads)
    #'{"3D":"tmp/HL-31_HPfixed_BPfixed/gradH/Colebrook/np_32/elasticity.exports/paraview.exports","Axi":"tmp/M19061901_laplace_dilatation31k/gradH/Montgomery/Colebrook/np_16/cfpdes.exports/paraview.exports"}'
    parser.add_argument("--name", type=str, help="input name", default="")
    parser.add_argument(
        "--views", help="activate views calculations", action="store_true"
    )
    parser.add_argument(
        "--stats", help="activate stats calculations", action="store_true"
    )
    parser.add_argument(
        "--histos", help="activate histograms calculations", action="store_true"
    )
    parser.add_argument(
        "--plots", help="activate plots calculations", action="store_true"
    )
    parser.add_argument("--r", help="plots vs r", action="store_true")
    parser.add_argument("--z", help="plots vs z", action="store_true")
    parser.add_argument("--theta", help="plots vs theta", action="store_true")
    parser.add_argument(
        "--cooling",
        help="what coolings do you want",
        type=str,
        choices=["mean", "meanH", "grad", "gradH", "gradHZ"],
        default="grad",
    )
    parser.add_argument(
        "--heatcorrelation",
        help="what heat correlations do you want",
        type=str,
        choices=["Montgomery", "Dittus", "Colburn", "Silverberg"],
        default="Montgomery",
    )
    parser.add_argument(
        "--friction",
        help="what frictions do you want",
        type=str,
        choices=["Constant", "Blasius", "Filonenko", "Colebrook", "Swanee"],
        default="Constant",
    )

    return parser


def key_dataframe(dirs, geos):
    with open(f"{dirs[0]}/FieldType.json", "r") as jsonfile:
        data1 = json.load(jsonfile)
    with open(f"{dirs[1]}/FieldType.json", "r") as jsonfile:
        data2 = json.load(jsonfile)

    types1 = {}
    for key, value in data1.items():
        types1[value["Type"]] = key

    types = []
    key1 = []
    key2 = []
    for key, value in data2.items():
        if value["Type"] in types1.keys():
            types.append(value["Type"])
            key1.append(types1[value["Type"]])
            key2.append(key)

    dfkeys = pd.DataFrame()
    dfkeys["types"] = types
    dfkeys["key1"] = key1
    dfkeys["key2"] = key2
    excluded = ["Stress", "ElectricPotential"]
    for type in excluded:
        dfkeys = dfkeys[dfkeys["types"] != type]
    print(dfkeys, flush=True)
    print("\n")

    return dfkeys


def get_files_list(measure: str, dir: str, rowkey: str, need: str = None):
    files1 = []
    conditions = {
        "plots": lambda filename1: filename1.endswith(".csv")
        and rowkey in filename1
        and need in filename1,
        "views": lambda filename1: filename1.endswith(".png")
        and rowkey in filename1
        and "deformed" not in filename1,
        "histograms": lambda filename1: filename1.endswith(".png")
        and rowkey in filename1
        and filename1.startswith("insert-"),
    }

    for filename1 in os.listdir(f"{dir}/{measure}"):
        if measure in conditions and conditions[measure](filename1):
            files1.append(os.path.join(f"{dir}/{measure}", filename1))

    return files1


def set_ax(
    ax,
    xlabel: str,
    ylabel: str,
    legend: list[str],
    title: str,
    xticks: list[float] = None,
):
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.legend(legend, fontsize=18, loc="best")
    ax.set_title(
        title,
        fontsize=20,
    )
    ax.grid(True, linestyle="--")
    if xticks:
        ax.set_xticks(xticks, xticks)
        ax.set_xlim([xticks[0], xticks[-1]])
    ax.tick_params(axis="x", which="major", labelsize=15)
    ax.tick_params(axis="y", which="major", labelsize=15)


def filter_files(files, exclude_terms=None, include_term=None, unique_term="norm"):
    if exclude_terms:
        files = [x for x in files if all(term not in x for term in exclude_terms)]
    if len(files) != 1:
        filetemp = [x for x in files if unique_term in x]
        if filetemp:
            files = filetemp
    if include_term:
        files = [x for x in files if include_term in x]
    return files


def merge_images(files: list[str], savefile: str):
    images = [Image.open(x) for x in files]
    widths, heights = zip(*(i.size for i in images))

    total_width = sum(widths)
    max_height = max(heights)

    new_im = Image.new("RGB", (total_width, max_height))

    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset, 0))
        x_offset += im.size[0]

    new_im.save(savefile)


def main():
    parser = options("", "")
    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    units = {
        "Temperature": "[°C]",
        "Displacement": "[mm]",
        "MagneticField": "[T]",
        "ElectricPotential": "[V]",
        "CurrentDensity": "[A/mm²]",
        "VonMises": "[MPa]",
        "Stress": "[MPa]",
    }

    basedir = f"res_compare/{args.name}/{args.cooling}/{args.friction}"
    os.makedirs(basedir, exist_ok=True)
    geos = []
    dirs = []
    for geo, dir in args.mdata.items():
        geos.append(geo)
        dirs.append(dir)

    dfkeys = key_dataframe(dirs, geos)

    if args.plots:
        print("Compare plots:", flush=True)
        os.makedirs(f"{basedir}/plots", exist_ok=True)
        if args.r:
            need = "vs-r"
        elif args.z:
            need = "vs-z"
        elif args.theta:
            need = "vs-theta"
        files = {}
        for index, row in dfkeys.iterrows():
            files1 = get_files_list("plots", dirs[0], row["key1"], need)
            files2 = get_files_list("plots", dirs[1], row["key2"], need)
            if files1 and files2:
                files[row["types"]] = {geos[0]: files1, geos[1]: files2}

        # print(files)

        show = False
        xticks = {}
        for index, values in files.items():
            fig, ax = plt.subplots(figsize=(12, 8))
            legend = []
            print(f"    plot {index}-{need}", flush=True)
            for i, v in values.items():
                for file in v:
                    df = pd.read_csv(file, index_col=0)
                    df.plot(x=df.columns[0], y=df.columns[1], linewidth=3, ax=ax)
                    legend.append(f'{i} {file.split(f"{need}-")[1].replace(".csv","")}')

            if args.r:
                xlabel = "r [m]"
            elif args.z:
                xlabel = "z [m]"
            elif args.theta:
                xlabel = "theta [deg]"
            ylabel = rf"{index} {units[index]}"
            title = f"{args.name} {args.cooling} {args.friction} {index} {geos[0]} vs {geos[1]}"
            set_ax(ax, xlabel, ylabel, legend, title)

            if show:
                fig.show()
            else:
                fig.tight_layout()
                fig.savefig(
                    f"{basedir}/plots/{args.name}-{args.cooling}-{args.friction}-{index}-{geos[0]}-{geos[1]}-{need}.png",
                    dpi=300,
                )
            xticks[index] = [round(x, 4) for x in ax.get_xticks()]

        if args.r and args.theta:
            if geos[0] in ["3D", "2D"]:
                geotheta = geos[0]
                dirtheta = dirs[0]
                keytheta = "key1"
            elif geos[1] in ["3D", "2D"]:
                geotheta = geos[1]
                dirtheta = dirs[1]
                keytheta = "key2"

            for index, row in dfkeys.iterrows():
                filestheta = get_files_list(
                    "plots", dirtheta, row[keytheta], "vs-theta"
                )

                if filestheta and row["types"] in files:
                    files[row["types"]][geotheta] = filestheta

            for index, values in files.items():
                fig, ax = plt.subplots(figsize=(12, 8))
                legend = []
                print(f"    plot {index}-vs-r w/ Boxplot-vs-theta", flush=True)
                for i, v in values.items():
                    if i in ["3D", "2D"]:
                        for file in v:
                            df = pd.read_csv(file, index_col=0)
                            r = float(file.split("r=")[-1].split("mm-z=")[0])
                            measure = file.split("/")[-1].split("-vs-theta")[0]
                            df["r"] = [r] * len(df[measure])
                            df.boxplot(
                                column=measure,
                                by=["r"],
                                positions=[r * 1e-3],
                                widths=0.01,
                                boxprops=dict(linestyle="-", linewidth=1.5, color="r"),
                                flierprops=dict(
                                    linestyle="-", linewidth=1.5, color="r"
                                ),
                                medianprops=dict(
                                    linestyle="-", linewidth=1.5, color="r"
                                ),
                                whiskerprops=dict(
                                    linestyle="-", linewidth=1.5, color="r"
                                ),
                                capprops=dict(linestyle="-", linewidth=1.5, color="r"),
                                showfliers=True,
                                ax=ax,
                            )
                            if legend:
                                legend.append("_nolegend_")
                            else:
                                legend.append(
                                    f'{i} z={file.split(f"z=")[1].replace(".csv","")}'
                                )
                    else:
                        for file in v:
                            df = pd.read_csv(file, index_col=0)
                            df.plot(
                                x=df.columns[0], y=df.columns[1], linewidth=3, ax=ax
                            )
                            legend.append(
                                f'{i} {file.split(f"{need}-")[1].replace(".csv","")}'
                            )

                xlabel = "r [m]"
                ylabel = rf"{index} {units[index]}"
                title = f"{args.name} {args.cooling} {args.friction} {index} {geos[0]} vs {geos[1]}"
                set_ax(ax, xlabel, ylabel, legend, title, xticks[index])
                if show:
                    fig.show()
                else:
                    fig.tight_layout()
                    fig.savefig(
                        f"{basedir}/plots/{args.name}-{args.cooling}-{args.friction}-{index}-{geos[0]}-{geos[1]}-vs-r-boxplot-theta.png",
                        dpi=300,
                    )
        print("\n", flush=True)

    if args.views:
        os.makedirs(f"{basedir}/views", exist_ok=True)
        print("Compare views:", flush=True)
        need = ""
        if args.z:
            need = "OxOy"
        elif args.theta:
            need = "OrOz"
        for index, row in dfkeys.iterrows():

            files1 = get_files_list("views", dirs[0], row["key1"])
            files2 = get_files_list("views", dirs[1], row["key2"])

            if files1 and files2:
                print(f"    views for {row['types']}", flush=True)
                if not need:
                    # Exclude files with OrOz or OxOy in their name
                    files1 = filter_files(files1, exclude_terms=["OrOz", "OxOy"])
                    files2 = filter_files(files2, exclude_terms=["OrOz", "OxOy"])
                else:
                    # Keep files with need in their name in geo is 3D
                    files1 = filter_files(
                        files1, include_term=need if geos[0] == "3D" else None
                    )
                    files2 = filter_files(
                        files2, include_term=need if geos[1] == "3D" else None
                    )

                # Ensure the "norm" filter is applied if length of files is not 1
                files1 = filter_files(files1)
                files2 = filter_files(files2)
                for file1 in files1:
                    for file2 in files2:
                        savefile = f"{basedir}/views/{args.name}-{args.cooling}-{args.friction}-{row['types']}-{file1.split('/')[-1].replace('.png','')}-{file2.split('/')[-1].replace('.png','')}.png"
                        merge_images([file1, file2], savefile)
        print("\n", flush=True)

    if args.histos:
        os.makedirs(f"{basedir}/histograms", exist_ok=True)
        print("Compare histograms:", flush=True)
        for index, row in dfkeys.iterrows():
            files1 = get_files_list("histograms", dirs[0], row["key1"])
            files2 = get_files_list("histograms", dirs[1], row["key2"])
            if files1 and files2:
                print(f"    histogram for {row['types']}", flush=True)
                # Ensure the "norm" filter is applied if length of files is not 1
                files1 = filter_files(files1)
                files2 = filter_files(files2)

                savefile = f"{basedir}/histograms/{args.name}-{args.cooling}-{args.friction}-{row['types']}-{files1[-1].split('/')[-1].replace('.png','')}-{files2[-1].split('/')[-1].replace('.png','')}.png"
                merge_images(files1 + files2, savefile)
        # print("\n", flush=True)


if __name__ == "__main__":
    sys.exit(main())
