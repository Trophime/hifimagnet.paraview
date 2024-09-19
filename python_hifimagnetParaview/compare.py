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
    ## Example '{"HL-31-3D":{"geo":"3D",""dir":"tmp/HL-31_HPfixed_BPfixed/gradH/Colebrook/np_32/elasticity.exports/paraview.exports"},"HL-31-Axi":{"geo":"Axi","dir":"tmp/M19061901_laplace_dilatation31k/gradH/Montgomery/Colebrook/np_16/cfpdes.exports/paraview.exports"}}'
    parser.add_argument(
        "--name", type=str, help="input result directory name", default=""
    )
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
        default="",
    )
    parser.add_argument(
        "--heatcorrelation",
        help="what heat correlations do you want",
        type=str,
        choices=["Montgomery", "Dittus", "Colburn", "Silverberg"],
        default="",
    )
    parser.add_argument(
        "--friction",
        help="what frictions do you want",
        type=str,
        choices=["Constant", "Blasius", "Filonenko", "Colebrook", "Swanee"],
        default="",
    )

    return parser


def key_dataframe(dirs: list[str]) -> pd.DataFrame:
    """create a DataFrame, containing the key names corresponding to physical measures.
    Only takes the measures in common in all the directories.

    Args:
        dirs (list[str]): list of directories that are compared

    Returns:
        pd.DataFrame: translator measures<->key names
    """
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
        if (
            value["Type"] in types1.keys()
        ):  ## Type of measure must be present in both dir
            types.append(value["Type"])
            key1.append(types1[value["Type"]])
            key2.append(key)

    dfkeys = pd.DataFrame()
    dfkeys["types"] = types
    dfkeys["key1"] = key1
    dfkeys["key2"] = key2

    ## Option for a third directory
    if len(dirs) == 3:
        with open(f"{dirs[2]}/FieldType.json", "r") as jsonfile:
            data3 = json.load(jsonfile)
        key3 = []
        for key, value in data3.items():
            if value["Type"] in types1.keys():
                key3.append(key)
        dfkeys["key3"] = key3

    excluded = ["Stress", "Strain"]
    for type in excluded:
        dfkeys = dfkeys[dfkeys["types"] != type]
    print(dfkeys, flush=True)
    print("\n")

    return dfkeys


def get_files_list(measure: str, dir: str, rowkey: str, need: str = None) -> list[str]:
    """get the list of files to compare for one measure in a directory for:
           - views comparison
           - plots comparison
           - histos comparison

    Args:
        measure (str): Type of measure to compare
        dir (str): directory of the results
        rowkey (str): key name for the type
        need (str, optional): string needed if specified (ex: 'vs-r'). Defaults to None.

    Returns:
        list[str]: return list of files to compare
    """
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
    """Formatting for the plot

    Args:
        ax: plt ax
        xlabel (str): label for x axis
        ylabel (str): label for y axis
        legend (list[str]): legend list
        title (str): title of plot
        xticks (list[float], optional): ticks for x axis. Defaults to None.
    """
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    if "Temperature" in ylabel and ax.get_ylim()[0] < 0:
        ax.set_ylim(bottom=0)
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


def filter_files(
    files: list[str],
    exclude_terms: list[str] = None,
    include_term: list[str] = None,
    unique_term: str = "norm",
) -> list[str]:
    """filter files from the list of files to compare

    Args:
        files (list[str]): list of files to compare
        exclude_terms (list[str], optional): strings which exclude the file from the list . Defaults to None.
        include_term (list[str], optional): strings which include the file from the list. Defaults to None.
        unique_term (str, optional): terme that need to be in the filename. Defaults to "norm".

    Returns:
        list[str]: return filtered files list
    """
    if exclude_terms:
        files = [x for x in files if all(term not in x for term in exclude_terms)]
    if len(files) != 1:
        filetemp = [x for x in files if unique_term in x.split("/")[-1]]
        if filetemp or unique_term != "norm":
            files = filetemp
    if include_term:
        files = [x for x in files if include_term in x]
    return files


def merge_images(files: list[str], savefile: str):
    """merge image for views comparison

    Args:
        files (list[str]): list of images names
        savefile (str): name of the comparison result image
    """
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
    color_cycler = [
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#F0E442",
        "c",
        "m",
        "y",
        "k",
    ]
    linestyle_cycler = [
        "-",
        "--",
        "-.",
        ":",
        "-",
        "--",
        "-.",
        "-",
        "--",
        "-.",
        ":",
        "-",
        "--",
        "-.",
    ]
    units = {
        "Temperature": "[°C]",
        "Displacement": "[mm]",
        "MagneticField": "[T]",
        "MagneticPotential": "[Wb/mm]",
        "ElectricPotential": "[V]",
        "ElectricField": "[V/m]",
        "ElectricConductivity": "[S/mm]",
        "ForceLaplace": "[N/mm³]",
        "Q": "[W/mm³]",
        "CurrentDensity": "[A/mm²]",
        "VonMises": "[MPa]",
        "Stress": "[MPa]",
        "Strain": "[]",
        "HoopStress": "[MPa]",
        "HoopStrain": "[]",
        "Stress_T": "[MPa]",
    }

    basedir = "res_compare"
    cooling = ""
    for arg in [args.name, args.cooling, args.heatcorrelation, args.friction]:
        if arg:
            basedir += f"/{arg}"
            if arg != args.name:
                cooling += f"-{arg}"
    os.makedirs(basedir, exist_ok=True)

    ## extract geometries directories & names from mdata
    geos = []
    dirs = []
    names = []
    for name, val in args.mdata.items():
        geos.append(val["geo"])
        dirs.append(val["dir"])
        names.append(name)

    name = f"{names[0]}_vs_{names[1]}"
    if len(dirs) == 3:
        name = name + f"_vs_{names[2]}"

    ## create measures translation dataframe
    dfkeys = key_dataframe(dirs)

    if args.plots:
        print("Compare plots:", flush=True)
        os.makedirs(f"{basedir}/plots", exist_ok=True)
        ## add needed str to find the plots we want
        if args.r:
            need = "vs-r"
        elif args.z:
            need = "vs-z"
        elif args.theta:
            need = "vs-theta"

        ## select the files to compare
        files = {}
        for index, row in dfkeys.iterrows():
            files1 = get_files_list("plots", dirs[0], row["key1"], need)
            files2 = get_files_list("plots", dirs[1], row["key2"], need)
            if len(dirs) == 3:
                files3 = get_files_list("plots", dirs[2], row["key3"], need)

            if len(files1) > 1 and len(files2) > 1:
                for suffix in ["_r", "_t", "_z", "_ur", "_ut"]:
                    files1_ = filter_files(files1, unique_term=suffix)
                    files2_ = filter_files(files2, unique_term=suffix)

                    if files1_ and files2_:
                        files[f'{row["types"]}{suffix}'] = {
                            names[0]: files1_,
                            names[1]: files2_,
                        }

            ## Ensure the "norm" filter is applied if length of files is not 1 (i.e. multiple files for 1 measure)
            files1 = filter_files(files1)
            files2 = filter_files(files2)
            if len(dirs) == 3:
                files3 = filter_files(files3)

            if files1 and files2:
                files[row["types"]] = {names[0]: files1, names[1]: files2}
                if len(dirs) == 3 and files3:
                    files[row["types"]][names[2]] = files3

        # print(files)

        ## plot every files in the same ax
        show = False
        xticks = {}
        for index, values in files.items():
            fig, ax = plt.subplots(figsize=(12, 8))
            legend = []
            print(f"    plot {index}-{need}", flush=True)
            col = 0
            for i, v in values.items():
                for file in v:
                    df = pd.read_csv(file, index_col=0)
                    print(f"        file: {file.split('/')[-1]}")
                    if "Jth" in file.split("/")[-1]:
                        df[df.columns[1]] = abs(df[df.columns[1]])
                    df.plot(
                        x=df.columns[0],
                        y=df.columns[1],
                        linewidth=3,
                        ax=ax,
                        color=color_cycler[col],
                        linestyle=linestyle_cycler[col],
                    )
                    legend.append(
                        f'{i} {file.split("/")[-1].split(f"{need}-")[1].replace(".csv","")}'
                    )
                    col = col + 1
                    if col == len(color_cycler):
                        col = 0
            ## Format the ax
            real_index = (
                index.replace("_r", "")
                .replace("_t", "")
                .replace("_z", "")
                .replace("_ur", "")
                .replace("_ut", "")
            )
            if args.r:
                xlabel = "r [m]"
            elif args.z:
                xlabel = "z [m]"
            elif args.theta:
                xlabel = "theta [deg]"
            ylabel = rf"{index} {units[real_index]}"
            title = f"{name}{cooling} {index} {geos[0]} vs {geos[1]}"
            set_ax(ax, xlabel, ylabel, legend, title)

            ## display or save the plot
            if show:
                fig.show()
            else:
                fig.tight_layout()
                fig.savefig(
                    f"{basedir}/plots/{name}{cooling}-{index}-{geos[0]}-{geos[1]}-{need}.png",
                    dpi=300,
                )
            xticks[index] = [round(x, 4) for x in ax.get_xticks()]

        if args.r and args.theta:
            ## WIP
            ## if one of the result compared in 3D or 2D, make boxplot given r for theta
            if geos[0] in ["3D", "2D"]:
                geotheta = geos[0]
                nametheta = names[0]
                dirtheta = dirs[0]
                keytheta = "key1"
            elif geos[1] in ["3D", "2D"]:
                geotheta = geos[1]
                nametheta = names[1]
                dirtheta = dirs[1]
                keytheta = "key2"

            for index, row in dfkeys.iterrows():
                filestheta = get_files_list(
                    "plots", dirtheta, row[keytheta], "vs-theta"
                )
                ## Ensure the "norm" filter is applied if length of files is not 1 (i.e. multiple files for 1 measure)
                filestheta = filter_files(filestheta)
                if filestheta and row["types"] in files:
                    files[row["types"]][nametheta] = filestheta

            for index, values in files.items():
                fig, ax = plt.subplots(figsize=(12, 8))
                legend = []
                print(f"    plot {index}-vs-r w/ Boxplot-vs-theta", flush=True)
                toboxplot = []
                for i, v in values.items():
                    if args.mdata[i]["geo"] in ["3D", "2D"]:
                        toboxplot.append((i, v))

                    else:
                        for file in v:
                            df = pd.read_csv(file, index_col=0)
                            print(f"        file: {file.split('/')[-1]}")
                            if "Jth" in file.split("/")[-1]:
                                df[df.columns[1]] = abs(df[df.columns[1]])
                            df.plot(
                                x=df.columns[0], y=df.columns[1], linewidth=3, ax=ax
                            )
                            legend.append(
                                f'{i} {file.split("/")[-1].split(f"{need}-")[1].replace(".csv","")}'
                            )

                for i, v in toboxplot:
                    for file in v:
                        df = pd.read_csv(file, index_col=0)
                        print(f"        file: {file.split('/')[-1]}")
                        r = float(file.split("/")[-1].split("r=")[-1].split("mm-z=")[0])
                        measure = file.split("/")[-1].split("-vs-theta")[0]
                        df["r"] = [r] * len(df[measure])
                        df.boxplot(
                            column=measure,
                            by=["r"],
                            positions=[r * 1e-3],
                            widths=0.01,
                            whis=[0, 100],
                            boxprops=dict(linestyle="-", linewidth=1.5, color="r"),
                            flierprops=dict(linestyle="-", linewidth=1.5, color="r"),
                            medianprops=dict(linestyle="-", linewidth=1.5, color="r"),
                            whiskerprops=dict(linestyle="-", linewidth=1.5, color="r"),
                            capprops=dict(linestyle="-", linewidth=1.5, color="r"),
                            showfliers=False,
                            ax=ax,
                        )
                        if file == v[0]:
                            legend.append(
                                f'{i} z={file.split("/")[-1].split(f"z=")[1].replace(".csv","")}'
                            )
                        else:
                            legend.append("_nolegend_")

                real_index = (
                    index.replace("_r", "")
                    .replace("_t", "")
                    .replace("_z", "")
                    .replace("_ur", "")
                    .replace("_ut", "")
                )
                xlabel = "r [m]"
                ylabel = rf"{index} {units[real_index]}"
                title = f"{name}{cooling} {index} {geos[0]} vs {geos[1]}"
                set_ax(ax, xlabel, ylabel, legend, title, xticks[index])
                if show:
                    fig.show()
                else:
                    fig.tight_layout()
                    fig.savefig(
                        f"{basedir}/plots/{name}{cooling}-{index}-{geos[0]}-{geos[1]}-vs-r-boxplot-theta.png",
                        dpi=300,
                    )
        print("\n", flush=True)

    if args.views:
        ## compare views (merge images)
        os.makedirs(f"{basedir}/views", exist_ok=True)
        print("Compare views:", flush=True)

        ## choose a view for 3D
        need = ""
        if args.z:
            need = "OxOy"
        elif args.theta:
            need = "OrOz"

        ## select files to display
        for index, row in dfkeys.iterrows():
            files1 = get_files_list("views", dirs[0], row["key1"])
            files2 = get_files_list("views", dirs[1], row["key2"])
            if len(dirs) == 3:
                files3 = get_files_list("views", dirs[2], row["key3"])
                print(files3)

            if files1 and files2:
                print(f"    views for {row['types']}", flush=True)
                if not need:
                    ## Exclude files with OrOz or OxOy in their name
                    files1 = filter_files(files1, exclude_terms=["OrOz", "OxOy"])
                    files2 = filter_files(files2, exclude_terms=["OrOz", "OxOy"])
                    if len(dirs) == 3:
                        files3 = filter_files(files3, exclude_terms=["OrOz", "OxOy"])
                else:
                    ## Keep files with need in their name in geo is 3D
                    files1 = filter_files(
                        files1, include_term=need if geos[0] == "3D" else None
                    )
                    files2 = filter_files(
                        files2, include_term=need if geos[1] == "3D" else None
                    )
                    if len(dirs) == 3:
                        files3 = filter_files(
                            files3, include_term=need if geos[2] == "3D" else None
                        )

                ## Ensure the "norm" filter is applied if length of files is not 1 (i.e. multiple files for 1 measure)
                files1 = filter_files(files1)
                files2 = filter_files(files2)
                if len(dirs) == 3:
                    files3 = filter_files(files3)

                ## merge the images
                for file1 in files1:
                    # print(f"        file1: {file1.split('/')[-1]}")
                    for file2 in files2:
                        # print(f"            file2: {file2.split('/')[-1]}")
                        if len(dirs) == 3:
                            for file3 in files3:
                                print(f"                file3: {file3.split('/')[-1]}")
                                savefile = f"{basedir}/views/{name}{cooling}-{row['types']}-{file1.split('/')[-1].replace('.png','')}-{file2.split('/')[-1].replace('.png','')}-{file3.split('/')[-1].replace('.png','')}.png"
                                merge_images([file1, file2, file3], savefile)
                        else:
                            savefile = f"{basedir}/views/{name}{cooling}-{row['types']}-{file1.split('/')[-1].replace('.png','')}-{file2.split('/')[-1].replace('.png','')}.png"
                            merge_images([file1, file2], savefile)
        print("\n", flush=True)

    if args.histos:
        ## merge histogram pictures
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

                savefile = f"{basedir}/histograms/{name}{cooling}-{row['types']}-{files1[-1].split('/')[-1].replace('.png','')}-{files2[-1].split('/')[-1].replace('.png','')}.png"
                merge_images(files1 + files2, savefile)
        # print("\n", flush=True)


if __name__ == "__main__":
    sys.exit(main())
