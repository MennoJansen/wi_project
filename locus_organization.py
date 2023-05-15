"""File for generating locus organization annotation file (iTOL).
In this file, it only works for the dhn-melanin biosynthetic pathway,
with the specified paralogues.

    Returns multiple annotation files (iTOL) in specified output directory.
"""
from collections import defaultdict
from ast import literal_eval

LONGEST_DIST = defaultdict(int)


class Gene:
    """Gene class, representing one single gene of one single species"""

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
        """Reverses order of group of genes and sets a bool to remember it is reversed"""
        self.group.reverse()
        self.swap = True
        self.distance = self.__init_distance(self.identifier, self.group)


class Species:
    """Species class with name, list of genes and counts of main genes"""

    def __init__(self, name: str):
        self.name = name
        self.genes = []
        self.gene_counts = [0, 0, 0, 0]


class AnnotationLine:
    """Annotation line class, allows generation of annotation strings"""

    def __init__(self, gene: Gene, file: str):
        self.gene = gene
        self.file = file
        self.dist_from_left = 0
        self.annotation_string = ""

    def set_distance_from_left(self, distance: int):
        """Set distance from left"""
        self.dist_from_left = distance

    def generate_annotation_string(self):
        """Generates a line for annotation file"""
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


SPECIES = {}
ANNOTATION_LINES = []
ANNOTATION_FILES = defaultdict(set)


def load_data2(info_file: str):
    """Loads data from info file, grouping genes correctly and
    adding all the generated Gene objects to Species objects"""
    with open(info_file, "r", encoding="utf-8") as info:
        for line in info:
            split_line = line.split(";")
            eval_3 = literal_eval(split_line[3])
            list_of_group_genes = []
            eval_5 = literal_eval(split_line[5].strip())
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


def get_shape(orientation: str, swap: bool) -> str:
    """Returns the shape of the gene, based on the
    orientation of scaffold and if the gene was reversed

    Args:
        orientation (str): orientation of scaffold
        swap (bool): if scaffold reversed

    Returns:
        str: shape of gene for annotation
    """
    if orientation == "+":
        if swap:
            return "TL"
        return "TR"
    if swap:
        return "TR"
    return "TL"


def get_color(type_gene: str) -> str:
    """Returns color of gene for visualization
    Colors here are based on different paralogues for different genes

    Args:
        type_gene (str): paralogue

    Returns:
        str: hex color code
    """
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


def count_genes():
    """Counts the amount of main genes."""
    for _, item in SPECIES.items():
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


def generate_annotation_file(info_file: str, output_path: str):
    """Generates annotation file

    Args:
        info_file (str): Info file for input
        output_path (str): path where output annotation files will be placed
    """
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
                            annotation_line = AnnotationLine(gene, key)
                            ANNOTATION_LINES.append(annotation_line)
        counter += 1

    # calculate and set distance from left for each file
    left_distances = defaultdict(int)
    for line in ANNOTATION_LINES:
        if left_distances[line.file] < line.gene.distance:
            left_distances[line.file] = line.gene.distance

    # make lines from line object and write to file
    for line in ANNOTATION_LINES:
        line.set_distance_from_left(left_distances[line.file])
        line.generate_annotation_string()

    # prepare output files with headers
    for file in ANNOTATION_FILES.keys():
        with open(
            f"{output_path}annotation_{file}.txt", "w", encoding="utf-8"
        ) as prep_output_file:
            prep_output_file.write(
                "DATASET_DOMAINS\nSEPARATOR"
                f" COMMA\nDATASET_LABEL,{file}\nCOLOR,#00ff00\n"
                f"HEIGHT_FACTOR,2.6\nWIDTH,{str(LONGEST_DIST[file])}\n\nDATA\n"
            )

    # write to output files
    for line in ANNOTATION_LINES:
        to_write = line.annotation_string.replace(
            "[length_of_group]", str(LONGEST_DIST[line.file])
        )
        with open(
            f"{output_path}annotation_{line.file}.txt", "a", encoding="utf-8"
        ) as output_file:
            output_file.write(to_write)
    print("Done!")


generate_annotation_file("/processing/jgi/locus/info.txt", "/processing/jgi/locus/")
