"""File generating a heatmap for iTOL species tree based on the occurrence
of certain genes of the dhn-melanin biosynthetic pathway.

Colors used:

0 4thn1 #9A1663
1 4thn2 #E0144C
2 athn #FF5858
3 ywa1 #FF97C1
4 ayg13 #A9907E
5 ayg12 #E5D9B6
6 ayg11 #A4BE7B
7 yg1 #5F8D4E
8 cladea #285430
9 4hnr #94DDFF
10 3hnr #0099e0
11 arp2 #AA77FF
12 scd1 #FFF38C
13 arp1 #C0B236

Input is an info file, structured like this (on single line):
gene_name;species_name;protein_id;(scaffold, pos_1, pos_2, orientation);
scaffold_length;[(species_name|protein_id, gene_name, pos_1, pos_2, orientation)]

The last list is the group of genes on the same scaffold.

Example:
ayg13;Asparx1;203013;('scaffold_84', 27956, 29337, '+');
85777;[('Asparx1|203011', 'arp2', 25517, 26370, '-'),
('Asparx1|203013', 'ayg13', 27956, 29337, '+')]
"""

from collections import defaultdict
import argparse
import os

HEATMAP = defaultdict(list)
GENES = {
    "4thn1": 0,
    "4thn2": 1,
    "athn": 2,
    "ywa1": 3,
    "ayg13": 4,
    "ayg12": 5,
    "ayg11": 6,
    "yg1": 7,
    "cladea": 8,
    "4hnr": 9,
    "3hnr": 10,
    "arp2": 11,
    "scd1": 12,
    "arp1": 13,
}


def read_data(inf_file: str):
    """Reads data of info file, processes it to heatmap variable (dictionary of lists)

    Args:
        info_file (str): info file, structured as described in file docstring
    """
    with open(inf_file, "r", encoding="utf-8") as info:
        for line in info:
            split_line = line.split(";")
            if len(HEATMAP[split_line[1]]) == 0:
                HEATMAP[split_line[1]] = [0] * 14
            HEATMAP[split_line[1]][GENES[split_line[0]]] = 1

    print("Done reading info!")


def write_annotation(ann_file: str):
    """Creates an annotation file ready for iTOL"""
    with open(ann_file, "w", encoding="utf-8") as open_ann_file:
        open_ann_file.write(
            "DATASET_EXTERNALSHAPE\nSEPARATOR"
            " COMMA\nDATASET_LABEL,heatmap_genes\nCOLOR,#00ffFF\n"
            "FIELD_COLORS,#9a1663,#e0144c,#ff5858,#ff97c1,#a9907e,"
            "#e5d9b6,#a4be7b,#5f8d4e,#285430,#94ddff,#0099e0,#aa77ff,#fff38c,#c0b236\n"
            "FIELD_LABELS,4thn1,4thn2,athn,ywa1,ayg13,ayg12,"
            "ayg11,yg1,cladea,4hnr,3hnr,arp2,scd1,arp1\n"
            "LEGEND_TITLE,Heatmap Genes\n"
            "LEGEND_SHAPES,1,1,1,1,1,1,1,1,1,1,1,1,1,1\n"
            "LEGEND_COLORS,#9a1663,#e0144c,#ff5858,#ff97c1,#a9907e,"
            "#e5d9b6,#a4be7b,#5f8d4e,#285430,#94ddff,#0099e0,#aa77ff,#fff38c,#c0b236\n"
            "LEGEND_LABELS,4thn1,4thn2,athn,ywa1,ayg13,ayg12,ayg11,"
            "yg1,cladea,4hnr,3hnr,arp2,scd1,arp1\n"
            "SHAPE_SPACING,0\nSHAPE_TYPE,1\n\nDATA\n"
        )
        for key, value in HEATMAP.items():
            open_ann_file.write(f"{key}")
            for item in value:
                open_ann_file.write(f",{item}")
            open_ann_file.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--info",
        dest="info_file",
        help=(
            "Input is an info file, structured like this (on single"
            " line):\ngene_name;species_name;protein_id;(scaffold, pos_1, pos_2,"
            " orientation);scaffold_length;[(species_name|protein_id, gene_name, pos_1,"
            " pos_2, orientation)]\n\nThe last list is the group of genes on the same"
            " scaffold.\n\nExample\n:ayg13;Asparx1;203013;('scaffold_84', 27956, 29337,"
            " '+');85777;[('Asparx1|203011', 'arp2', 25517, 26370,"
            " '-'),('Asparx1|203013', 'ayg13', 27956, 29337, '+')]"
        ),
    )
    parser.add_argument(
        "-a",
        "--annotation",
        dest="annotation_file",
        help="Name and location of output annotation file",
    )
    args = parser.parse_args()
    if args.info_file is None:
        print("Info file not specified...")
        exit()
    if args.annotation_file is None:
        print("Annotation file not specified...")
        exit()

    info_file = args.info_file
    assert os.path.exists(info_file), "Info file not found..."
    annotation_file = args.annotation_file

    read_data(info_file)
    write_annotation(annotation_file)
