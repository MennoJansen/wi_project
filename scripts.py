import argparse
import zipfile
import os
import datetime
import dateutil.parser
import re
from operator import itemgetter


def get_species_list(input_file: str, output_file: str):
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


def remove_end(input_file: str, output_file: str):
    species = []
    with open(input_file, "r", encoding="utf-8") as f_in:
        for line in f_in:
            species.append(line.split("|")[0])

    with open(output_file, "w", encoding="utf-8") as f_out:
        for item in species:
            f_out.write(f"{item}\n")


def get_name(input_file: str, output_file: str):
    species = []
    with open(input_file, "r", encoding="utf-8") as f_in:
        for line in f_in:
            split_line = line.split("|")
            if not split_line[0].startswith("\n"):
                if split_line[0].startswith("CHARACTER") and len(split_line) <= 3:
                    name = split_line[1].replace(" ", "_")
                    species.append(f"{name}")
                else:
                    name = split_line[2].replace(" ", "_")
                    species.append(f"{name}|{split_line[3]}")

    with open(output_file, "w", encoding="utf-8") as f_out:
        for item in species:
            f_out.write(f"{item}\n")


def get_seqs(fasta_file: str, species_file: str, output_file: str):
    species = []
    with open(species_file, "r", encoding="utf-8") as species_in:
        for line in species_in:
            species.append(line.strip())

    seqs = dict()
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


def add_classes(classes_file: str, names_file: str, output_file: str):
    classes = dict()
    output = dict()
    with open(classes_file, "r", encoding="utf-8") as classes_in:
        for line in classes_in:
            split_line = line.split(",")
            classes[split_line[0]] = split_line[1].strip()

    with open(names_file, "r", encoding="utf-8") as names_in:
        for line in names_in:
            split_line = line.split("|")
            try:
                output[line.strip()] = classes[split_line[0]]
            except KeyError:
                output[line.strip()] = "None"

    with open(output_file, "w", encoding="utf-8") as f_out:
        for key, value in output.items():
            f_out.write(f"{key},{value}\n")


def get_fasta_gff(hmm_file: str, busco: str, gff_file: str, fasta_file: str):
    prot_name = ""
    with open(hmm_file, "r", encoding="utf-8") as f_hmm:
        for line in f_hmm:
            if busco in line:
                prot_name = line.split()[0].split("|")[3]

    scaffold = ""
    orientation = ""
    min_num = 999999999999999
    max_num = 0
    if prot_name != "":
        with open(gff_file, "r", encoding="utf-8") as f_gff:
            for line in f_gff:
                if prot_name in line:
                    split_line = line.split()
                    num_1 = int(split_line[3])
                    num_2 = int(split_line[4])
                    if num_1 < num_2:
                        min_num = min(num_1, min_num)
                        max_num = max(num_2, max_num)
                    scaffold = split_line[0]
                    orientation = split_line[6]

    seq = ""
    if scaffold != "":
        with open(fasta_file, "r", encoding="utf-8") as f_fasta:
            next_line = False
            for line in f_fasta:
                if scaffold in line:
                    next_line = True
                    continue
                if next_line:
                    seq = line[min_num:max_num]
                    if orientation == "+":
                        seq = seq.upper().strip()
                else:
                    tab = str.maketrans("ACTG", "TGAC")
                    seq = seq.upper().translate(tab).strip()[::-1]
    return seq


def get_gff_file(annotation_file: str, gff_path: str):
    """Get gff file if no gff3 file exists"""
    with zipfile.ZipFile(annotation_file) as zip_file:
        for path in zip_file.namelist():
            file_dir, file_name = path.split("/")
            if not os.path.exists(f"{gff_path}{file_dir}.gff3") and not os.path.exists(
                f"{gff_path}{file_dir}.gff"
            ):
                print(f"{file_dir} - {file_name}")


def remove_duplicate_lines(input_file: str, output_file: str):
    lines = []
    with open(input_file, "r", encoding="utf-8") as i_file:
        for line in i_file.readlines():
            lines.append(line.split(",")[0].split("|")[0])
    unique_lines = set(lines)

    with open(output_file, "w", encoding="utf-8") as o_file:
        for line in unique_lines:
            o_file.write(f"{line}\n")


def get_locations(
    gene_tree_list_file: str,
    gff_directory: str,
    output_file: str,
):
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
                # if "Settu3" not in species_name:
                #    continue
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
                        with open("locus/error.txt", "a", encoding="utf-8") as error:
                            error.write(
                                f"{species_name}:{prot_name};wrong_protein_id\n"
                            )
                        continue
                except FileNotFoundError as e:
                    with open("locus/error.txt", "a", encoding="utf-8") as error:
                        error.write(f"{species_name}:{prot_name};wrong_version\n")
                    continue
            if location[0] == "":
                print(f"Could not find location for following line: {line}")
            else:
                with open(output_file, "a", encoding="utf-8") as o_file:
                    o_file.write(f"{gene};{species_name};{prot_id};{location}\n")


def filter_date_species(file_list_file: str, filter_date: datetime.datetime):
    dates = {}
    with open(file_list_file, "r", encoding="utf-8") as list_file:
        for line in list_file:
            split_line = line.split(":")
            name = split_line[0].strip()
            second_split_line = line.split(";")
            third_split_line = second_split_line[1].split()
            date = f"{third_split_line[1]} {third_split_line[2]} {third_split_line[5]}"
            parsed_date = dateutil.parser.parse(date)
            if name not in dates.keys():
                if parsed_date < filter_date:
                    dates[name] = parsed_date
                    # print(f"{name} | {dates[name]}")
                # dates[name] =
            else:
                if dates[name] < parsed_date and parsed_date < filter_date:
                    dates[name] = parsed_date
    fungi = {}
    with open("fungies.txt", "r", encoding="utf-8") as fungi_file:
        for line in fungi_file:
            split_line = line.split(",")
            name = split_line[1].replace('"', "").strip()
            id = split_line[2]
            fungi[name] = id

    with open("filtered_species.txt", "w", encoding="utf-8") as filter_file:
        for item in dates:
            filter_file.write(f"{fungi[item]}\n")


def rename_nodes(tree_file: str, new_file_name: str):
    p = re.compile(r"[A-Z]*\|jgi\|([A-Za-z0-9_]*)\|\d*\|[A-Za-z0-9_.-]*:")
    # print(p.match("DOTHIDEO|jgi|Melpu1|411912|estExt_fgenesh1_pm.C_7090002:").group(1))

    with open(tree_file, "r", encoding="utf-8") as tree:
        for line in tree:
            new_line = p.sub(r"\g<1>:", line)

            with open(new_file_name, "w", encoding="utf-8") as output_tree:
                output_tree.write(new_line)


def edit_bootstrap(tree_file: str, new_file_name: str):
    p = re.compile(r"(\d{1,3}\.?\d?)\/\d{1}\.?\d{0,3}\/\d{1,3}:(\d{1}.\d*)")

    with open(tree_file, "r", encoding="utf-8") as tree:
        for line in tree:
            new_line = p.sub(r":\g<1>", line)

            with open(new_file_name, "w", encoding="utf-8") as output_tree:
                output_tree.write(new_line)


def filter_lists(filtered_list_file: str, unfiltered_list_file: str, output_file: str):
    unfiltered_items = {}
    with open(unfiltered_list_file, "r", encoding="utf-8") as unfiltered_list:
        for item in unfiltered_list:
            split_item = item.split("|")
            if split_item[0] != "CHARACTERIZED":
                key = f"{split_item[0]}|{split_item[1]}|{split_item[2]}|{split_item[3]}"
                unfiltered_items[key] = split_item[4]

    filtered_items = {}
    with open(filtered_list_file, "r", encoding="utf-8") as filtered_list:
        for item in filtered_list:
            split_item = item.split("|")
            if split_item[1] != "CHARACTERIZED":
                name = split_item[3].replace(" ", "_")
                key = f"{split_item[1]}|{split_item[2]}|{name}|{split_item[4]}"
                filtered_items[f"{split_item[0]}|{key}"] = unfiltered_items[key]

    with open(output_file, "w", encoding="utf-8") as output:
        for key, value in filtered_items.items():
            output.write(f"{key}|{value}")


def get_scaffold_info(locations_file: str, info_file: str):
    items = []
    gene_locations = {}
    scaffold_lengths = {}
    fasta_path = "/processing/jgi/nucl/"
    with open(locations_file, "r", encoding="utf-8") as locations:
        max_length = 0
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
            if key.startswith("Pengl1"):
                print("break")
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

            gff3_file = f"/drive/gff/{name}.gff3"
            if os.path.exists(gff3_file):
                with open(gff3_file, "r", encoding="utf-8") as gff3:
                    for line in gff3:
                        if line.startswith("##"):
                            split_line = line.split()
                            if split_line[1] == scaffold:
                                length = split_line[3]
            else:
                fasta_file = f"/processing/jgi/nucl/{name}.fasta"
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


"""
get_species_list(
    "/processing/jgi/AYG1_trimal_all.fasta", "/processing/jgi/ayg1_species.txt"
)
"""

"""
get_diff_species(
    "/processing/jgi/jgi_unique_species.txt",
    "/processing/jgi/gene_species.txt",
    "/processing/jgi/hit_species_correct.txt",
)
"""
"""
remove_end(
    "/processing/jgi/ayg1_hit_species.txt", "/processing/jgi/ayg1_hit_species2.txt"
)
"""
# get_name("/processing/jgi/SCD/scd1_species.txt", "/processing/jgi/SCD/scd1_names.txt")
"""
get_seqs(
    "/processing/jgi/SCD/SCD1_blast_all.fasta",
    "/processing/jgi/SCD/scd1_names.txt",
    "/processing/jgi/SCD/scd1_seqs.fasta",
)
"""
"""
add_classes(
    "/processing/jgi/classes.txt",
    "/processing/jgi/SCD/scd1_names.txt",
    "/processing/jgi/SCD/scd1_classes.txt",
)
"""
# get_gff_file("/home/menno/Downloads/annotation.zip", "/drive/gff/")
# remove_duplicate_lines("/processing/jgi/jgi_hits.csv", "jgi_unique_species.txt")

# For Locus organization
"""
get_locations(
    "/processing/jgi/ayg1_labels.txt",
    "/drive/gff/",
    "/processing/jgi/locations.txt",
)
"""

get_locations(
    "/processing/jgi/locus/output.txt",
    "/drive/gff/",
    "/processing/jgi/locus/locations.txt",
)


get_scaffold_info(
    "/processing/jgi/locus/locations.txt", "/processing/jgi/locus/info.txt"
)

"""
filter_date_species(
    "/processing/jgi/file_list.txt", dateutil.parser.parse("Jan 12 2021")
)
"""
# For Notung
"""
rename_nodes("subtrees/ayg_cladea.txt", "subtrees/r_ayg_cladea.txt")
edit_bootstrap("subtrees/r_ayg_cladea.txt", "subtrees/p_ayg_cladea.txt")
"""
"""
filter_lists(
    "/processing/jgi/locus/filtered_list_genes_names.txt",
    "/processing/jgi/locus/unfiltered_list_genes.txt",
    "/processing/jgi/locus/output.txt",
)
"""
