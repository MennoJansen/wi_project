import argparse


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


"""
get_species_list(
    "/processing/jgi/AYG1_trimal_all.fasta", "/processing/jgi/ayg1_species.txt"
)
"""

"""
get_diff_species(
    "/processing/jgi/ayg1_hit_species2.txt",
    "/processing/jgi/species.txt",
    "/processing/jgi/missing_species_speciestree.txt",
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
add_classes(
    "/processing/jgi/classes.txt",
    "/processing/jgi/SCD/scd1_names.txt",
    "/processing/jgi/SCD/scd1_classes.txt",
)
