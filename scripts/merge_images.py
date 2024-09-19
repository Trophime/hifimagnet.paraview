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
    parser.add_argument(
        "--file1",
        help="first file",
        type=str,
    )
    parser.add_argument(
        "--file2",
        help="second file",
        type=str,
    )
    parser.add_argument(
        "--name", type=str, help="input result directory name", default=""
    )

    return parser


def merge_images(files: list[str], savefile: str):
    """merge image for views comparison

    Args:
        files (list[str]): list of images names
        savefile (str): name of the comparison result image
    """
    images = [Image.open(x) for x in files]
    widths, heights = zip(*(i.size for i in images))

    max_height = max(heights)
    min_height = min(heights)

    for im in images:
        width, height = im.size
        ratio = width / height
        if height > min_height:
            new_width = int(ratio * min_height)
            im.thumbnail((new_width, min_height), Image.Resampling.LANCZOS)

    widths, heights = zip(*(i.size for i in images))
    total_width = sum(widths)
    # new_im = Image.new("RGB", (total_width, max_height))
    new_im = Image.new("RGB", (total_width, min_height))

    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset, 0))
        x_offset += im.size[0]

    new_im.save(savefile)


def main():
    parser = options("", "")
    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    os.makedirs(f"{args.name}", exist_ok=True)

    savefile = f"{args.name}/{args.file1.split('/')[-1].replace('.png','')}-{args.file2.split('/')[-1].replace('.png','')}.png"
    merge_images([args.file1, args.file2], savefile)
    # print("\n", flush=True)


if __name__ == "__main__":
    sys.exit(main())
