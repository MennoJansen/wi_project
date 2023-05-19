import zipfile
import os
from operator import itemgetter


def get_species_list(input_file: str, output_file: str):
    """Gets species list from list of nodes of gene tree

    Args:
        input_file (str): list of nodes of gene tree
        output_file (str): list of species
    """
    species = []
    with open(input_file, "r", encoding="utf-8") as f_in:
        for line in f_in:
            if line.startswith(">") and not line.startswith(">CHARAC"):
                specie = line.split("|")[2]
                species.append(specie)

    with open(output_file, "w", encoding="utf-8") as f_out:
        for specie in species:
            f_out.write(f"{specie}\n")


def get_diff_species(input_file_1: str, input_file_2: str, output_file: str):
    """Given two species lists, returns the difference between the two sets

    Args:
        input_file_1 (str): species list 1
        input_file_2 (str): species list 2
        output_file (str): difference of species list 1 and 2
    """
    species_1 = []
    species_2 = []
    with open(input_file_1, "r", encoding="utf-8") as f1_in:
        for line in f1_in:
            species_1.append(line.strip())
    with open(input_file_2, "r", encoding="utf-8") as f2_in:
        for line in f2_in:
            species_2.append(line.strip())

    species_3 = set(species_1).difference(species_2)

    with open(output_file, "w", encoding="utf-8") as f_out:
        for item in species_3:
            f_out.write(f"{item}\n")


def get_seqs(fasta_file: str, species_file: str, output_file: str):
    species = []
    with open(species_file, "r", encoding="utf-8") as species_in:
        for line in species_in:
            species.append(line.strip())

    seqs = {}
    with open(fasta_file, "r", encoding="utf-8") as fasta_in:
        seq = ""
        name = ""
        add_seq = False
        for line in fasta_in:
            if line.startswith(">"):
                if seq != "":
                    seqs[name] = seq
                    seq = ""
                    name = ""
                    add_seq = False
                split_line = line.split("|")
                if line.startswith(">CHARACTER") and len(split_line) <= 2:
                    name = (
                        line.split("|")[1].replace("[", "_").replace("]", "_").strip()
                    )
                    if name in species:
                        add_seq = True
                else:
                    split_line = line.split("|")
                    name = f"{split_line[2]}|{split_line[3]}"
                    if name in species:
                        add_seq = True
                continue
            if add_seq:
                seq += line.strip()
    with open(output_file, "w", encoding="utf-8") as f_out:
        for key, value in seqs.items():
            name = key.split("|")[0]
            f_out.write(f">{key}\n")
            f_out.write(f"{value}\n")


def get_gff_file(annotation_file: str, gff_path: str):
    """Get gff file if no gff3 file exists"""
    with zipfile.ZipFile(annotation_file) as zip_file:
        for path in zip_file.namelist():
            file_dir, file_name = path.split("/")
            if not os.path.exists(f"{gff_path}{file_dir}.gff3") and not os.path.exists(
                f"{gff_path}{file_dir}.gff"
            ):
                print(f"{file_dir} - {file_name}")


def get_locations(
    gene_tree_list_file: str,
    gff_directory: str,
    output_file: str,
):
    """With the input of a file containing a list of hits from the gene trees, structured like this:
    ayg13|EUROTIO|jgi|Talpro1|286554|e_gw1.1.1905.1

    Gets the locations of all these hits from gff or gff3 files,
    and puts them in an output file structured like this:
    ayg13;Talpro1;286554;('scaffold_1', 3769008, 3770491, '-')

    In parantheses, it is (scaffold, lowest_pos, highest_pos, orientation)

    All the inputs that are not found, are reported in error.txt,
    with either wrong protein id or outdated version.
    So for CocheC5, it has multiple versions, and the gene trees
    used for the internship were dated, so the location could
    not be found in the new version of the genome. This should
    not be an issue if you create new gene trees.

    Args:
        gene_tree_list_file (str): _description_
        gff_directory (str): _description_
        output_file (str): _description_
    """
    with open(gene_tree_list_file, "r", encoding="utf-8") as gene_tree_list:
        for line in gene_tree_list:
            split_line = line.split("|")
            species_name = ""
            prot_name = ""
            location = ("", 0, 0, "")
            gff_3 = False
            # Get species name and prot name
            if not line.startswith("\n"):
                gene = split_line[0]
                species_name = split_line[3]
                prot_id = split_line[4]
                print(f"Running for species: {gene}|{species_name}")
                if os.path.exists(f"{gff_directory}{species_name}.gff3"):
                    # GFF3 file
                    gff_3 = True
                    prot_name = f"proteinId={split_line[4]}"
                else:
                    # GFF file
                    gff_3 = False
                    prot_name = split_line[5].strip()
                lower_prot_name = prot_name.lower()
                # Get location from gff3/gff file
                try:
                    if gff_3:
                        with open(
                            f"{gff_directory}{species_name}.gff3", "r", encoding="utf-8"
                        ) as gff_file:
                            for line in gff_file:
                                if (
                                    f"{prot_name};" in line
                                    or f"{prot_name}\n" in line
                                    or f"{lower_prot_name};" in line
                                    or f"{lower_prot_name}\n" in line
                                ) and "Parent" in line:
                                    split_line = line.split()
                                    num_1 = int(split_line[3])
                                    num_2 = int(split_line[4])
                                    if num_1 < num_2:
                                        location = (
                                            split_line[0],
                                            num_1,
                                            num_2,
                                            split_line[6],
                                        )
                                    else:
                                        location = (
                                            split_line[0],
                                            num_2,
                                            num_1,
                                            split_line[6],
                                        )
                    else:
                        with open(
                            f"{gff_directory}{species_name}.gff", "r", encoding="utf-8"
                        ) as gff_file:
                            for line in gff_file:
                                if (
                                    f'"{prot_name}";' in line
                                    or f'"{prot_name}"\n' in line
                                ):
                                    lowest = 999999999999999999999
                                    highest = 0
                                    split_line = line.split()
                                    num_1 = int(split_line[3])
                                    num_2 = int(split_line[4])
                                    lowest = min(num_1, num_2, lowest)
                                    highest = max(num_1, num_2, highest)
                                    location = (
                                        split_line[0],
                                        lowest,
                                        highest,
                                        split_line[6],
                                    )
                    if location[0] == "":
                        with open("error.txt", "a", encoding="utf-8") as error:
                            error.write(
                                f"{species_name}:{prot_name};wrong_protein_id\n"
                            )
                        continue
                except FileNotFoundError:
                    with open("error.txt", "a", encoding="utf-8") as error:
                        error.write(f"{species_name}:{prot_name};wrong_version\n")
                    continue
            if location[0] == "":
                print(f"Could not find location for following line: {line}")
            else:
                with open(output_file, "a", encoding="utf-8") as o_file:
                    o_file.write(f"{gene};{species_name};{prot_id};{location}\n")


def get_scaffold_info(
    locations_file: str, gff_path: str, fasta_path: str, info_file: str
):
    """Using the output of the get_locations method, this method gathers more information about the clustering of the genes and some scaffold information using the gff3 files or fasta files (containing the whole genome of the species). It is considered a cluster when two or more genes are located on the same scaffold. Filtering for distance between the genes is done at a later stage.

    The output is written into an info file, structured like this:
    ayg13;Asparx1;203013;('scaffold_84', 27956, 29337, '+');85777;[('Asparx1|203011', 'arp2', 25517, 26370, '-'), ('Asparx1|203013', 'ayg13', 27956, 29337, '+')]

    It looks similar to the locations output, but additionnaly, we have a scaffold length (in this case 85777) and cluster information (an arp2 gene and an ayg13 gene on the same scaffold)
    """
    items = []
    gene_locations = {}
    scaffold_lengths = {}
    with open(locations_file, "r", encoding="utf-8") as locations:
        for line in locations:
            split_line = line.split(";")
            gene = split_line[0]
            name = split_line[1]
            prot_id = split_line[2]
            print(f"Running for {gene}:{name}")
            location = eval(split_line[3])
            scaffold = location[0]
            start = location[1]
            end = location[2]
            orientation = location[3]

            # get cluster organization
            key = f"{name}|{scaffold}"
            identifier = f"{name}|{prot_id}"
            gene_loc = (identifier, gene, start, end, orientation)
            if key in gene_locations:
                current_gene_locs = gene_locations[key]
                current_gene_locs.append(gene_loc)
                gene_locations[key] = current_gene_locs
            else:
                gene_locations[key] = [gene_loc]

            # Get scaffold length
            length = 0

            gff3_file = f"{gff_path}{name}.gff3"
            if os.path.exists(gff3_file):
                with open(gff3_file, "r", encoding="utf-8") as gff3:
                    for line in gff3:
                        if line.startswith("##"):
                            split_line = line.split()
                            if split_line[1] == scaffold:
                                length = split_line[3]
            else:
                fasta_file = f"{fasta_path}{name}.fasta"
                next_line = False
                with open(fasta_file, "r", encoding="utf-8") as fasta:
                    for line in fasta:
                        if next_line:
                            length = len(line)
                            next_line = False
                        if line == f">{scaffold}\n":
                            next_line = True
            scaffold_lengths[key] = length

        # Sort clusters on position on genome
        for key in gene_locations.keys():
            copy_list = gene_locations[key]
            sorted_list = sorted(copy_list, key=itemgetter(1))
            gene_locations[key] = sorted_list

    # Prepare new info for output
    with open(locations_file, "r", encoding="utf-8") as locations:
        for line in locations:
            split_line = line.split(";")
            name = split_line[1]
            location = eval(split_line[3])
            scaffold = location[0]
            key = f"{name}|{scaffold}"
            items.append(
                f"{line.strip()};{scaffold_lengths[key]};{gene_locations[key]}\n"
            )

    # write output to file
    with open(info_file, "w", encoding="utf-8") as output:
        for item in items:
            output.write(item)


# Below you can find how to use the locus organization functions
"""
get_locations(
    "/processing/jgi/locus/output.txt",
    "/drive/gff/",
    "/processing/jgi/locus/locations.txt",
)


get_scaffold_info(
    "/processing/jgi/locus/locations.txt",
    "/drive/gff/",
    "/processing/jgi/nucl/",
    "/processing/jgi/locus/info.txt",
)
"""
