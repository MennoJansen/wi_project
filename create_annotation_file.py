"""Generates an annotation file with color strips for the respective classes
    
    Parameters
    ----------
    file : str
        path to csv containing all species in column 1 and all respective classes in column 2.
"""
import argparse
import os

COLORS = {
    "Eurotiomycetes": "#ce7e00",
    "Dothideomycetes": "#8fce00",
    "Lecanoromycetes": "#2986cc",
    "Leotiomycetes": "#6a329f",
    "Sordariomycetes": "#6dbbcc",
    "Pezizomycetes": "#f44336",
    "Xylonomycetes": "#990000",
    "Orbiliomycetes": "#744700",
    "Peziz_Incertae": "#000000",
    "Taph_Incertae": "#7bffbd",
    "Saccharomycetes": "#ecff0e",
    "Neolectomycetes": "#cc00ff",
    "Pneumocystidomycetes": "#8bc34a",
    "Schizosaccharomycetes": "#fad0c3",
    "Taphrinomycetes": "#d4cfb",
    "Coniocybomycetes": "#bedadc",
}


def create_annotation_file(file: str):
    """Creates an annotation file using input csv, color coding the species based on their class
    This file is named annotation.txt and will be placed in the same directory as this script

    Args:
        file (str): _description_
    """
    with open("annotation.txt", "w", encoding="utf-8") as ann_file:
        ann_file.write("DATASET_COLORSTRIP\n")
        ann_file.write("SEPARATOR COMMA\n")
        ann_file.write("DATASET_LABEL,classes\n")
        ann_file.write("COLOR,#ff0000\n")
        ann_file.write("COLOR_BRANCHES,1\n\n")

        ann_file.write("DATA\n\n")
        with open(file, encoding="utf-8") as csv_file:
            for line in csv_file.readlines():
                split_line = line.split(",")
                name = split_line[0]
                color = COLORS[split_line[1].strip()]
                ann_file.write(f"{name},{color}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help="Input file in CSV format, containing species in column 1 and classes in column 2",
    )

    args = parser.parse_args()

    if args.input is None:
        print("Input file not specified")
        exit()

    input_file = args.input

    assert os.path.exists(input_file), "Input file could not be found"

    create_annotation_file(input_file)
