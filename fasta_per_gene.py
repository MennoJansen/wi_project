""" This script creates a fasta file per gene
It retrieves the gene ID from the HMMsearch output
This gene ID is used to retrieve the location info (scaffold, start, end) from GFF files
This location info is used to extract the corresponding sequence
Input: .txt file with GENES you want to do, HMMsearch output, GFF files and fasta files
Output: Fasta file per gene
"""

# Import packages
import os
import logging
import argparse
from dataclasses import dataclass


@dataclass
class Paths:
    """Contains paths to .gff files, .fasta files, .hmm files and output files"""

    gff: str
    fasta: str
    hmm: str
    out: str

    def __init__(self, gff, fasta, hmm, out):
        self.gff = gff
        self.fasta = fasta
        self.hmm = hmm
        self.out = out


@dataclass
class GeneLocation:
    """Contains location attributes of gene"""

    scaffold: str
    orientation: str
    start: str
    end: str

    def __init__(self, scaffold, orientation, start, end):
        self.scaffold = scaffold
        self.orientation = orientation
        self.start = start
        self.end = end


def get_gene_name(hmmsearch_output: list[str], gene: str) -> str:
    """gets gene name or ID in select species

    Args:
        hmmsearch_output (list[str]): hmm output of single species hmm search
        gene (str): busco gene name

    Returns:
        str: gene id of species
    """
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


def get_location_gff(gff_output: list[str], gene_name: str) -> GeneLocation:
    """Gets the location of the gene from the gff file

    Args:
        gff_output (list[str]): gff output of file
        gene_name (str): gene name

    Returns:
        GeneLocation: gene location with scaffold, orientation, start and end
    """
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

    return GeneLocation(scaffold, orientation, start, end)


def get_sequence_gff(
    g_fasta_output: list[str],
    output_file: str,
    gene_location: GeneLocation,
    species: str,
) -> None:
    """_summary_

    Args:
        g_fasta_output (list[str]): _description_
        output_file (str): _description_
        gene_location (GeneLocation): _description_
        species (str): _description_
    """
    next_line = False
    for line in g_fasta_output:
        if next_line:
            seq = line[int(gene_location.start) : int(gene_location.end)]
            if gene_location.orientation == "+":
                output_file.write(seq.upper().strip() + "\n")
                return
            else:
                tab = str.maketrans("ACTG", "TGAC")
                output_file.write(seq.upper().translate(tab).strip()[::-1] + "\n")
                return
        if gene_location.scaffold in line:
            next_line = True
            output_file.write(f">{species}\n")


def get_fasta_gff(hmm_file: str, busco: str, gff_file: str, fasta_file: str) -> str:
    """Retrieves sequences data of gene in fasta file

    Args:
        hmm_file (str): hmm file
        busco (str): busco gene name
        gff_file (str): gff file
        fasta_file (str): fasta file

    Returns:
        str: sequence data of gene
    """
    prot_name = ""
    with open(hmm_file, "r", encoding="utf-8") as f_hmm:
        for line in f_hmm:
            if busco in line:
                prot_name = line.split()[0].split("|")[3]

    gene_location = GeneLocation("", "", 99999999999, 0)
    gene_location.start = 999999999999999
    gene_location.end = 0
    if prot_name != "":
        with open(gff_file, "r", encoding="utf-8") as f_gff:
            for line in f_gff:
                split_line = line.split()
                gff_prot_name = split_line[9].split('"')[1]
                if prot_name == gff_prot_name:
                    num_1 = int(split_line[3])
                    num_2 = int(split_line[4])
                    if num_1 < num_2:
                        gene_location.start = min(num_1, gene_location.start)
                        gene_location.end = max(num_2, gene_location.end)
                    else:
                        gene_location.start = min(num_2, gene_location.start)
                        gene_location.end = max(num_1, gene_location.end)
                    gene_location.scaffold = split_line[0]
                    gene_location.orientation = split_line[6]

    seq = ""
    if gene_location.scaffold != "":
        with open(fasta_file, "r", encoding="utf-8") as f_fasta:
            next_line = False
            for line in f_fasta:
                if gene_location.scaffold in line:
                    next_line = True
                    continue
                if next_line:
                    seq = line[gene_location.start : gene_location.end]
                    if gene_location.orientation == "+":
                        seq = seq.upper().strip()
                    else:
                        tab = str.maketrans("ACTG", "TGAC")
                        seq = seq.upper().translate(tab).strip()[::-1]
                    break
    return seq


def create_fasta_per_gene(gene: str, species_list: list[str], paths: Paths):
    """Generates fasta files as output, for each specie the specified
    gene is looked up, and the sequence data is added to the output file.

    Args:
        gene (str): the gene for which the program needs to find the sequence data in each species
        species_list (list[str]): List of species
        paths (Paths): Paths to .gff files, .fasta files, .hmm files and output files
    """
    species_counter = 1
    species_counter_end = len(species_list)
    seq = ""
    for item in species_list:
        logging.info(
            "Running for SPECIES %s (%s/%s)", item, species_counter, species_counter_end
        )
        species_counter += 1
        output_file_path = f"{paths.out}{gene}.txt"
        if not os.path.exists(f"{paths.gff}{item}.gff3"):
            seq = get_fasta_gff(
                f"{paths.hmm}hmm_out_{item}.txt",
                gene,
                f"{paths.gff}{item}.gff",
                f"{paths.fasta}{item}.fasta",
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
                f"{paths.hmm}hmm_out_{item}.txt", "r", encoding="utf-8"
            ) as hmmsearch_output:
                gene_name = get_gene_name(hmmsearch_output, gene)
                if gene_name == "not_found":
                    output_file = open(output_file_path, "a", encoding="utf-8")
                    output_file.write(f">{item}\n\n")
                    output_file.close()
                else:
                    with open(
                        f"{paths.gff}{item}.gff3", "rt", encoding="utf-8"
                    ) as gff_output:
                        gene_location = get_location_gff(gff_output, gene_name)
                        with open(
                            f"{paths.fasta}{item}.fasta", "r", encoding="utf-8"
                        ) as g_fasta_output:
                            g_fasta_output = g_fasta_output.readlines()
                            output_file = open(output_file_path, "a", encoding="utf-8")
                            get_sequence_gff(
                                g_fasta_output,
                                output_file,
                                gene_location,
                                item,
                            )
                            output_file.close()


def fasta_gene_generator(
    genes: list[str],
    species: list[str],
    paths: Paths,
):
    """Calls create_fasta_per_gene for each gene, adds logging

    Args:
        genes (list[str]): list of genes
        species (list[str]): list of species
        paths (Paths): contains paths to .gff, .fasta, .hmm and output files
    """
    gene_counter = 1
    gene_counter_end = len(genes)
    for gene in genes:
        gene = gene.strip()
        logging.info(
            "Running for gene: %s (%s/%s)", gene, gene_counter, gene_counter_end
        )
        create_fasta_per_gene(gene, species, paths)


def main():
    """Wrapper for main"""
    parser = argparse.ArgumentParser()

    parser.add_argument("species", dest="species", help="List of species")
    parser.add_argument(
        "gff_path", dest="gff", help="Path to folder containing .gff or .gff3 files"
    )
    parser.add_argument(
        "fasta_path", dest="fasta", help="Path to folder containing .fasta files"
    )
    parser.add_argument(
        "hmm_path", dest="hmm", help="Path to folder containing .hmm files"
    )
    parser.add_argument("out_path", dest="out", help="Path to output folder")
    parser.add_argument("genes", "--genes", dest="genes", help="List of genes")

    args = parser.parse_args(
        [
            "species.txt",
            "/drive/gff/",
            "/processing/jgi/nucl/",
            "/processing/jgi/hmm_output5/",
            "/processing/jgi/out_fastas_per_gene_yanfang/",
            "/processing/jgi/buscos_left.txt",
        ]
    )

    assert os.path.exists(args.species), "No species file found"
    assert os.path.exists(args.gff), "Invalid path to gff folder"
    assert os.path.exists(args.fasta), "Invalid path to fasta folder"
    assert os.path.exists(args.hmm), "Invalid path to hmm folder"
    assert os.path.exists(args.out), "Invalid path to output folder"
    assert os.path.exists(args.genes), "No genes file found"

    paths = Paths(args.gff, args.fasta, args.hmm, args.out)

    with open(args.genes, encoding="utf-8") as genes_list_file:
        genes = genes_list_file.readlines()

    logging.basicConfig(
        filename="7_fasta_per_gene.log",
        encoding="utf-8",
        level=logging.DEBUG,
        format="%(asctime)s %(message)s",
    )
    logging.getLogger().addHandler(logging.StreamHandler())

    with open(args.species, encoding="utf-8") as species_file:
        species_list = []
        for i in species_file:
            species_list.append(i.strip())

    fasta_gene_generator(genes, species_list, paths)


if __name__ == "__main__":
    main()
