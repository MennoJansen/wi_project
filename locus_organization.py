from collections import defaultdict
from typing import Self


LONGEST_DIST = defaultdict(int)


class Gene:
    def __init_main_gene(self, gene_name) -> str:
        if (
            gene_name == "4thn1"
            or gene_name == "4thn2"
            or gene_name == "athn"
            or gene_name == "ywa1"
        ):
            return "pks"
        if (
            gene_name == "ayg13"
            or gene_name == "ayg12"
            or gene_name == "ayg11"
            or gene_name == "yg1"
            or gene_name == "cladea"
        ):
            return "ayg"
        if gene_name == "3hnr" or gene_name == "4hnr" or gene_name == "arp2":
            return "rdt"
        if gene_name == "scd1" or gene_name == "arp1":
            return "scd"
        return "unknown"

    def __init_distance(self, identifier, group) -> int:
        counter = 1
        for item in group:
            if f"{item[1]}|{item[0]}" == identifier:
                return counter
            counter += 1

    def __init__(
        self,
        species_name: str,
        gene_name: str,
        prot_id: int,
        pos1: int,
        pos2: int,
        scaffold: str,
        orientation: str,
        group: list[tuple[str, str, int, int, str]],
    ):
        self.species_name = species_name
        self.gene_name = gene_name
        self.prot_id = prot_id
        self.identifier = f"{gene_name}|{species_name}|{prot_id}"
        self.main_gene = self.__init_main_gene(gene_name)

        self.pos1 = pos1
        self.pos2 = pos2
        self.scaffold = scaffold
        self.orientation = orientation
        self.group = group
        if len(self.group) > 1:
            self.group.sort(key=lambda x: x[2])
        self.drawn = False
        self.swap = False

        self.distance = self.__init_distance(self.identifier, self.group)

    def reverse_group(self):
        self.group.reverse()
        self.swap = True
        self.distance = self.__init_distance(self.identifier, self.group)

    def recount_gene_distance(self):
        self.distance = self.__init_distance(self.identifier, self.group)


class Species:
    def __init__(self, name: str):
        self.name = name
        self.genes = []
        self.gene_counts = [0, 0, 0, 0]


class Annotation_Line:
    def __init__(self, gene: Gene, file: str):
        self.gene = gene
        self.file = file
        self.dist_from_left = 0
        self.annotation_string = ""

    def set_distance_from_left(self, distance: int):
        self.dist_from_left = distance

    def generate_annotation_string(self):
        # length_of_group = len(self.gene.group)
        first_pos = self.dist_from_left - self.gene.distance
        shapes = [0] * first_pos
        for gene in self.gene.group:
            shape = get_shape(gene[4], self.gene.swap)
            color = get_color(gene[1])
            shapes.append((shape, color, gene[1]))

        annotation_string = f"{self.gene.species_name},[length_of_group]"

        current_dist = 10
        for shape in shapes:
            if shape == 0:
                # nothing, empty
                current_dist += 50 + 10
            else:
                line = (
                    f"{shape[0]}|{current_dist}|{current_dist + 50}|{shape[1]}|{shape[2]}"
                )
                annotation_string = f"{annotation_string},{line}"
                current_dist += 50 + 10
        self.annotation_string = f"{annotation_string}\n"
        if current_dist > LONGEST_DIST[self.file]:
            LONGEST_DIST[self.file] = current_dist


INFO_DATA2 = []
DRAWN2 = defaultdict(bool)
COUNTS = defaultdict(int)
GENES_LIST = []
SPECIES_GENES = defaultdict(list)
SPECIES = {}
ANNOTATION_LINES = []
ANNOTATION_FILES = defaultdict(set)


def load_data(info_file: str):
    global INFO_DATA
    global DRAWN
    DRAWN = defaultdict(bool)
    INFO_DATA = []
    with open(info_file, "r", encoding="utf-8") as info:
        for line in info:
            INFO_DATA.append(line.strip())


def load_data2(info_file: str):
    with open(info_file, "r", encoding="utf-8") as info:
        for line in info:
            split_line = line.split(";")
            eval_3 = eval(split_line[3])
            list_of_group_genes = []
            eval_5 = eval(split_line[5].strip())
            counter = 0
            for item in eval_5:
                if item[0] != f"{split_line[1]}|{split_line[2]}":
                    if counter != 0:
                        if item[2] - 50000 > eval_5[counter - 1][2]:
                            list_of_group_genes.append(item)
                    else:
                        list_of_group_genes.append(item)
                else:
                    if counter != 0:
                        if item[2] - 50000 > eval_5[counter - 1][2]:
                            list_of_group_genes.pop()
                            list_of_group_genes.append(item)
                    else:
                        list_of_group_genes.append(item)
            new_gene = Gene(
                split_line[1],
                split_line[0],
                split_line[2],
                eval_3[1],
                eval_3[2],
                eval_3[0],
                eval_3[3],
                list_of_group_genes,
            )
            if not new_gene.species_name in SPECIES:
                SPECIES[new_gene.species_name] = Species(new_gene.species_name)
            SPECIES[new_gene.species_name].genes.append(new_gene)

    print("Done reading data!")


def create_annotation_line(
    main_gene: tuple[str, str, tuple[str, int, int, str]],
    species: str,
    genes: list[tuple[str, int, int, str]],
    max_dist_from_left: dict[str, int],
) -> tuple[str, int]:
    # Check if swapping is needed
    for gene in genes:
        swap = False
        if check_gene(gene) == "pks":
            if gene[4] == "-":
                swap = True

    # reverse gene list if pks is in reverse
    if swap:
        genes.reverse()

    # get pks distance from left
    pks_distance = 1
    for gene in genes:
        if check_gene(gene) == "pks":
            break
        else:
            pks_distance += 1

    # make shape at right distance for all genes
    max_length = 0
    annotation_entry = f"{species},[length_of_group]"
    count = 0
    for gene in genes:
        distance = (
            50 + count * 150 + ((max_dist_from_left[main_gene[0]] - pks_distance) * 150)
        )
        type_gene = check_gene(gene)
        if type_gene == "pks":
            line = f"TR|{distance}|{distance+100}|{get_color(gene[1])}|{gene[1]}"
        else:
            shape = get_shape(gene[4], swap)
            color = get_color(type_gene)
            line = f"{shape}|{distance}|{distance+100}|{color}|{gene[1]}"
        max_length = distance
        annotation_entry = f"{annotation_entry},{line}"
        count += 1
        DRAWN[gene[0]] = True
    annotation_entry = f"{annotation_entry}\n"

    return (annotation_entry, max_length + 150)


def replace_placeholder_lengths(
    annotation_entries: list[str], max_length: int
) -> list[str]:
    final_annotation = []
    for item in annotation_entries:
        final_annotation.append(item.replace("[length_of_group]", str(max_length)))
    return final_annotation


def get_dist_from_left(
    all_genes: list[list[tuple[str, int, int, str]]]
) -> dict[str, int]:
    max_dist = defaultdict(int)
    for genes in all_genes:
        dist = 1
        for gene in genes:
            if check_gene(gene) == "pks":
                if gene[4] == "-":
                    dist = len(genes) - dist
                if dist > max_dist[gene[1]]:
                    max_dist[gene[1]] = dist
            else:
                dist += 1
    return max_dist


def get_shape(orientation: str, swap: bool) -> str:
    if orientation == "+":
        if swap:
            return "TL"
        return "TR"
    if swap:
        return "TR"
    return "TL"


def get_color(type_gene: str) -> str:
    match type_gene:
        case "3hnr":
            return "#0099e0"
        case "4hnr":
            return "#94DDFF"
        case "arp2":
            return "#AA77FF"
        case "scd1":
            return "#FFF38C"
        case "arp1":
            return "#C0B236"
        case "ayg13":
            return "#A9907E"
        case "ayg12":
            return "#E5D9B6"
        case "ayg11":
            return "#A4BE7B"
        case "yg1":
            return "#5F8D4E"
        case "cladea":
            return "#285430"
        case "4thn1":
            return "#9A1663"
        case "4thn2":
            return "#E0144C"
        case "athn":
            return "#FF5858"
        case "ywa1":
            return "#FF97C1"


def check_gene(gene: tuple[str, str, int, int, str]) -> str:
    if (
        gene[1] == "4thn1"
        or gene[1] == "4thn2"
        or gene[1] == "athn"
        or gene[1] == "ywa1"
    ):
        return "pks"
    if (
        gene[1] == "ayg13"
        or gene[1] == "ayg12"
        or gene[1] == "ayg11"
        or gene[1] == "yg1"
        or gene[1] == "cladea"
    ):
        return "ayg"
    if gene[1] == "3hnr" or gene[1] == "4hnr" or gene[1] == "arp2":
        return "rdt"
    if gene[1] == "scd1" or gene[1] == "arp1":
        return "scd"
    return "unknown"


def count_genes():
    for key, item in SPECIES.items():
        counts = [0, 0, 0, 0]
        # 0 = pks, 1 = ayg, 2 = rdt and 3 = scd
        for gene in item.genes:
            if gene.main_gene == "pks":
                counts[0] += 1
            if gene.main_gene == "ayg":
                counts[1] += 1
            if gene.main_gene == "rdt":
                counts[2] += 1
            if gene.main_gene == "scd":
                counts[3] += 1
        item.gene_counts = counts


def check_gene2(gene: str) -> str:
    if gene == "4thn1" or gene == "4thn2" or gene == "athn" or gene == "ywa1":
        return "pks"
    if (
        gene == "ayg13"
        or gene == "ayg12"
        or gene == "ayg11"
        or gene == "yg1"
        or gene == "cladea"
    ):
        return "ayg"
    if gene == "3hnr" or gene == "4hnr" or gene == "arp2":
        return "rdt"
    if gene == "scd1" or gene == "arp1":
        return "scd"
    return "unknown"


def calc_distance_from_left():
    for file in ANNOTATION_FILES:
        max_distance = 0
        for item in ANNOTATION_LINES:
            if item.file == file:
                if item.gene.distance > max_distance:
                    max_distance = item.gene.distance

        for item in ANNOTATION_LINES:
            item.set_distance_from_left = max_distance


def generate_annotation_file(info_file: str, output_path: str):
    annotation_dict = defaultdict(list)
    # load_data(info_file)
    load_data2(info_file)
    count_genes()

    # 1 PKS, 2 AYG, 3 RDT, 4 SCD
    current_gene = [
        "4thn1",
        "4thn2",
        "ywa1",
        "athn",
        "ayg13",
        "ayg12",
        "ayg11",
        "yg1",
        "cladea",
        "3hnr",
        "4hnr",
        "arp2",
        "scd1",
        "arp1",
    ]
    counter = 0
    while counter < 13:
        for species in SPECIES.values():
            for gene in species.genes:
                if gene.gene_name == current_gene[counter]:
                    file_counter = 0
                    while gene.drawn is False:
                        key = f"{current_gene[counter]}_{file_counter}"
                        draw = False
                        if key in ANNOTATION_FILES:
                            if species.name in ANNOTATION_FILES[key]:
                                file_counter += 1
                            else:
                                ANNOTATION_FILES[key].add(species.name)
                                draw = True
                        else:
                            ANNOTATION_FILES[key] = {species.name}
                            draw = True
                        if draw is True:
                            if gene.orientation == "-":
                                gene.reverse_group()
                            gene.drawn = True
                            for gene2 in species.genes:
                                for item in gene.group:
                                    if gene2.identifier == f"{item[1]}|{item[0]}":
                                        gene2.drawn = True
                            annotation_line = Annotation_Line(gene, key)
                            ANNOTATION_LINES.append(annotation_line)
        counter += 1

    # calculate and set distance from left for each file
    left_distances = defaultdict(int)
    for line in ANNOTATION_LINES:
        if left_distances[line.file] < line.gene.distance:
            left_distances[line.file] = line.gene.distance

    # make lines from line object and write to file
    for line in ANNOTATION_LINES:
        # if len(line.gene.group) > 1:
        line.set_distance_from_left(left_distances[line.file])
        line.generate_annotation_string()

    # prepare output files with headers
    for file in ANNOTATION_FILES.keys():
        with open(
            f"locus/annotation_{file}.txt", "w", encoding="utf-8"
        ) as prep_output_file:
            prep_output_file.write(
                "DATASET_DOMAINS\nSEPARATOR"
                f" COMMA\nDATASET_LABEL,{file}\nCOLOR,#00ff00\nHEIGHT_FACTOR,2.6\nWIDTH,{str(LONGEST_DIST[file])}\n\nDATA\n"
            )

    # write to output files
    for line in ANNOTATION_LINES:
        to_write = line.annotation_string.replace(
            "[length_of_group]", str(LONGEST_DIST[line.file])
        )
        with open(
            f"locus/annotation_{line.file}.txt", "a", encoding="utf-8"
        ) as output_file:
            output_file.write(to_write)

    all_genes = []
    for item in INFO_DATA:
        split_item = item.split(";")
        all_genes.append(eval(split_item[5]))

    max_dist = get_dist_from_left(all_genes)
    # DEAL WITH DOUBLE COPIES OF SAME PKS GENE THAT ARE SPLIT
    max_group_length = defaultdict(int)
    for item in INFO_DATA:
        split_item = item.split(";")
        if not DRAWN[f"{split_item[1]}|{split_item[2]}"]:
            if check_gene(("", split_item[0], 0, 0, "")) == "pks":
                line, group_length = create_annotation_line(
                    (split_item[0], split_item[1], split_item[3]),
                    split_item[1],
                    eval(split_item[5]),
                    max_dist,
                )
                if group_length > max_group_length[split_item[0]]:
                    max_group_length[split_item[0]] = group_length
                annotation_dict[split_item[0]].append(line)
    for key in annotation_dict.keys():
        annotation_dict[key] = replace_placeholder_lengths(
            annotation_dict[key], max_group_length[key]
        )
        with open(
            f"{output_path}annotation_{key}.txt", "w", encoding="utf-8"
        ) as output:
            output.write(
                "DATASET_DOMAINS\nSEPARATOR"
                f" COMMA\nDATASET_LABEL,{key}\nCOLOR,#00ff00\nHEIGHT_FACTOR,2.6\n\nDATA\n"
            )
            for line in annotation_dict[key]:
                output.write(line)
    print("Done!")


generate_annotation_file("/processing/jgi/locus/info.txt", "/processing/jgi/locus/")
