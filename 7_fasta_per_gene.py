""" This script creates a fasta file per gene
It retrieves the gene ID from the HMMsearch output
This gene ID is used to retrieve the location info (scaffold, start, end) from GFF files
This location info is used to extract the corresponding sequence
Input: .txt file with genes you want to do, HMMsearch output, GFF files and fasta files
Output: Fasta file per gene
"""

# Import packages
import os
import re
import logging

# Open relevant directories
gff_path = "/drive/gff/"
g_fasta_path = "/processing/jgi/nucl/"
hmmsearch_path = "/processing/jgi/hmm_output5/"
output_path = "/processing/jgi/output_fastas_per_gene_yanfang/"

gff_files = os.listdir(gff_path)
g_fasta_files = os.listdir(g_fasta_path)
hmmsearch_files = os.listdir(hmmsearch_path)

# Regex for hmm output files
p = re.compile(".*_out_(.*).*")

# Open .txt files with genes you need the sequence from
# genes = open('/processing/jgi/hmm_intersection/buscos.txt').readlines()
genes = open("/processing/jgi/buscos_left.txt", encoding="utf-8").readlines()
counter = 1
counter_end = len(genes)

logging.basicConfig(
    filename="7_fasta_per_gene.log",
    encoding="utf-8",
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
)
logging.getLogger().addHandler(logging.StreamHandler())

# Open species you need the sequence from
species = open("specieslist_1324.txt", encoding="utf-8").readlines()
s_counter = 1
s_counter_end = len(species)
species_list = []
for i in species:
    species_list.append(i.strip())


def get_gene_name(hmmsearch_output, gene):
    target_name = "not_found"
    for line in hmmsearch_output:
        if not line.startswith("#"):
            split_line = line.split("-")
            query_name = split_line[1].strip()
            if query_name == gene:
                target_name = split_line[0].strip()
                id_mrna = f"ID={target_name};"
                return id_mrna
    return target_name


def get_location(gff_output, gene_name):
    scaffold, start, end, orientation = "0", "0", "0", "0"
    for line in gff_output:
        split_gene_name = "".join(gene_name.split("_"))
        if gene_name in line or split_gene_name in line:
            scaffold = line.split()[0]
            start = line.split()[3]
            end = line.split()[4]
            orientation = line.split()[6]
            break

    if int(end) < int(start):
        new_end = start
        new_start = end
        start = new_start
        end = new_end
    if gene_name != "not_found" and start == "0":
        logging.warning("WARNING!")
    start = str(int(start) - 1)

    fasta_header = ">" + scaffold + ":" + start + "-" + end
    return fasta_header, orientation


def get_location_gff(gff_output, gene_name):
    scaffold, start, end, orientation = "0", "0", "0", "0"
    first = True
    for line in gff_output:
        if gene_name in line:
            if first is True:
                scaffold = line.split()[0]
                start = line.split()[3]
                end = line.split()[4]
                orientation = line.split()[6]
                first = False
            else:
                end = line.split()[4]
        elif not gene_name in line and first is False:
            break

    if int(end) < int(start):
        new_end = start
        new_start = end
        start = new_start
        end = new_end

    if gene_name != "not_found" and start == "0":
        logging.warning("WARNING!")
    start = str(int(start) - 1)

    return scaffold, orientation, start, end


def get_sequence(g_fasta_output, fasta_header, output_file, orientation, species):
    for i in range(0, len(g_fasta_output) - 1, 2):
        if g_fasta_output[i].strip() == fasta_header:
            output_file.write(f">{species}\n")

            # Take the complement sequence
            # + 1 in indexing is needed due to a +1 shift in counting
            if orientation == "+":
                output_file.write(g_fasta_output[i + 1].upper())
                return
            elif orientation == "-":
                tab = str.maketrans("ACTG", "TGAC")
                output_file.write(
                    g_fasta_output[i + 1].upper().translate(tab).strip()[::-1] + "\n"
                )
                return
            else:
                print('other orientation than "+" and "-"')


def get_sequence_gff(
    g_fasta_output, output_file, scaffold, orientation, start, end, species
):
    next_line = False
    for line in g_fasta_output:
        if next_line:
            seq = line[int(start) : int(end)]
            if orientation == "+":
                output_file.write(seq.upper().strip() + "\n")
                return
            else:
                tab = str.maketrans("ACTG", "TGAC")
                output_file.write(seq.upper().translate(tab).strip()[::-1] + "\n")
                return
        if scaffold in line:
            next_line = True
            output_file.write(f">{species}\n")


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
                split_line = line.split()
                gff_prot_name = split_line[9].split('"')[1]
                if prot_name == gff_prot_name:
                    num_1 = int(split_line[3])
                    num_2 = int(split_line[4])
                    if num_1 < num_2:
                        min_num = min(num_1, min_num)
                        max_num = max(num_2, max_num)
                    else:
                        min_num = min(num_2, min_num)
                        max_num = max(num_1, max_num)
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
                    break
    return seq


def create_fasta_per_gene(gene):
    global counter, s_counter
    logging.info("Running for gene: %s (%s/%s)", gene, counter, counter_end)
    counter = counter + 1
    s_counter = 1
    seq = ""
    for item in species_list:
        logging.info("Running for species %s (%s/%s)", item, s_counter, s_counter_end)
        s_counter += 1
        output_file_path = f"{output_path}{gene}.txt"
        if not os.path.exists(f"{gff_path}{item}.gff3"):
            seq = get_fasta_gff(
                f"{hmmsearch_path}hmm_out_{item}.txt",
                gene,
                f"{gff_path}{item}.gff",
                f"{g_fasta_path}{item}.fasta",
            )
            if seq == "":
                output_file = open(output_file_path, "a", encoding="utf-8")
                output_file.write(f">{item}\n\n")
                output_file.close()
            else:
                output_file = open(output_file_path, "a", encoding="utf-8")
                output_file.write(f">{item}\n")
                output_file.write(f"{seq}\n")
                output_file.close()
        else:
            with open(
                f"{hmmsearch_path}hmm_out_{item}.txt", "r", encoding="utf-8"
            ) as hmmsearch_output:
                gene_name = get_gene_name(hmmsearch_output, gene)
                if gene_name == "not_found":
                    output_file = open(output_file_path, "a", encoding="utf-8")
                    output_file.write(f">{item}\n\n")
                    output_file.close()
                else:
                    with open(
                        f"{gff_path}{item}.gff3", "rt", encoding="utf-8"
                    ) as gff_output:
                        get_loc = get_location_gff(gff_output, gene_name)
                        scaffold, orientation, start, end = (
                            get_loc[0],
                            get_loc[1],
                            get_loc[2],
                            get_loc[3],
                        )
                        with open(
                            f"{g_fasta_path}{item}.fasta", "r", encoding="utf-8"
                        ) as g_fasta_output:
                            g_fasta_output = g_fasta_output.readlines()
                            output_file = open(output_file_path, "a", encoding="utf-8")
                            get_sequence_gff(
                                g_fasta_output,
                                output_file,
                                scaffold,
                                orientation,
                                start,
                                end,
                                item,
                            )
                            output_file.close()


if __name__ == "__main__":
    for item in genes:
        item = item.strip()
        create_fasta_per_gene(item)
