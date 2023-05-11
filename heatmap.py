from collections import defaultdict

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

# 0 4thn1 #9A1663
# 1 4thn2 #E0144C
# 2 athn #FF5858
# 3 ywa1 #FF97C1
# 4 ayg13 #A9907E
# 5 ayg12 #E5D9B6
# 6 ayg11 #A4BE7B
# 7 yg1 #5F8D4E
# 8 cladea #285430
# 9 4hnr #94DDFF
# 10 3hnr #0099e0
# 11 arp2 #AA77FF
# 12 scd1 #FFF38C
# 13 arp1 #C0B236


def read_data(info_file: str):
    with open(info_file, "r", encoding="utf-8") as info:
        for line in info:
            split_line = line.split(";")
            if len(HEATMAP[split_line[1]]) == 0:
                HEATMAP[split_line[1]] = [0] * 14
            HEATMAP[split_line[1]][GENES[split_line[0]]] = 1

    print("Done reading info!")


def write_annotation():
    with open(
        "locus/annotation_heatmap_genes.txt", "w", encoding="utf-8"
    ) as annotation_file:
        annotation_file.write(
            f"DATASET_EXTERNALSHAPE\nSEPARATOR"
            f" COMMA\nDATASET_LABEL,heatmap_genes\nCOLOR,#00ffFF\n"
            f"FIELD_COLORS,#9a1663,#e0144c,#ff5858,#ff97c1,#a9907e,#e5d9b6,#a4be7b,#5f8d4e,#285430,#94ddff,#0099e0,#aa77ff,#fff38c,#c0b236\n"
            f"FIELD_LABELS,4thn1,4thn2,athn,ywa1,ayg13,ayg12,ayg11,yg1,cladea,4hnr,3hnr,arp2,scd1,arp1\n"
            f"LEGEND_TITLE,Heatmap Genes\n"
            f"LEGEND_SHAPES,1,1,1,1,1,1,1,1,1,1,1,1,1,1\n"
            f"LEGEND_COLORS,#9a1663,#e0144c,#ff5858,#ff97c1,#a9907e,#e5d9b6,#a4be7b,#5f8d4e,#285430,#94ddff,#0099e0,#aa77ff,#fff38c,#c0b236\n"
            f"LEGEND_LABELS,4thn1,4thn2,athn,ywa1,ayg13,ayg12,ayg11,yg1,cladea,4hnr,3hnr,arp2,scd1,arp1\n"
            f"SHAPE_SPACING,0\nSHAPE_TYPE,1\n\nDATA\n"
        )
        for key, value in HEATMAP.items():
            annotation_file.write(f"{key}")
            for item in value:
                annotation_file.write(f",{item}")
            annotation_file.write("\n")


read_data("locus/info.txt")
write_annotation()
