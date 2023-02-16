"""Concats fasta files using a species list and a genes list to also add
empty sequences when a species is missing a certain gene

Looks for input files named: trimal_{gene}.fasta

Species names have to match those in input files

"""
import os
import argparse


def concat_fastas(
    input_path: str, output_path: str, species_file: str, buscos_file: str
):
    """Concats fasta files using a species list and a genes list to also add
    empty sequences when a species is missing a certain gene

    Args:
        input_path (str): path to input folder containing fasta files
        output_path (str): path to output file which will be the concatenated fasta files
        species_file (str): path to species file containing all species separated by new lines
        buscos_file (str): path to genes file containing all genes separated by new lines
    """
    new_fasta = dict()
    with open(species_file, encoding="utf-8") as all_species:
        for species in all_species:
            species = species.strip()
            new_fasta[species] = ""

    with open(buscos_file, encoding="utf-8") as genes:
        total_length = 0
        for gene in genes:
            # First, scan gene fasta file, save sequences in new_fasta dict
            gene_fasta_file = f"{input_path}trimal_{gene.strip()}.fasta"
            len_name = "undefined"
            print("Running for gene " + gene.strip())
            with open(gene_fasta_file, encoding="utf-8") as gene_fastas:
                name = "undefined"
                for line in gene_fastas:
                    # If line found is a specie name, set in var name
                    # Else, we are in a sequence of a specie, so we add this sequence
                    # to the existing sequence of this specie
                    if line.startswith(">"):
                        name = line.split(">")[
                            1
                        ].strip()  # >Aspcost1\n becomes Aspcost1
                        len_name = name
                    else:
                        new_fasta[name] = new_fasta[name] + line.strip()

            # Second, add empty sequences for species where the BUSCO gene is not present
            seq_length = len(new_fasta[len_name]) - total_length
            total_length = len(new_fasta[len_name])
            print(total_length)
            for name, seq in new_fasta.items():
                if len(seq) != total_length:
                    empty_seq = "-" * seq_length
                    new_fasta[name] = new_fasta[name] + empty_seq

    with open(output_path, "w", encoding="utf-8") as output_file:
        for name, seq in new_fasta.items():
            output_file.write(">" + name + "\n")
            output_file.write(seq + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help=(
            "Input path to directory containing trimmed and aligned .fasta files to be"
            " concatenated"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="Output path to file containing concatenated sequences",
    )
    parser.add_argument(
        "-s",
        "--species",
        dest="species",
        help="Path to file containing list of species, separated by newlines",
    )
    parser.add_argument(
        "-g", "--genes", dest="genes", help="Path to file containing genes to concat"
    )

    args = parser.parse_args(
        [
            "-i",
            "/processing/jgi/output_fastas_per_gene_v3/",
            "-o",
            "/processing/jgi/output_fastas_per_gene_v3/cat_trimal.fasta",
            "-s",
            "/processing/jgi/species.txt",
            "-g",
            "/processing/jgi/buscos.txt",
        ]
    )

    if args.input is None:
        print("Input path not specified")
        exit()
    if args.output is None:
        print("Output path not specified")
        exit()
    if args.species is None:
        print("Species path not specified")
        exit()
    if args.genes is None:
        print("Genes path not specified")
        exit()

    __input_path = args.input
    __output_path = args.output
    __species_path = args.species
    __genes_path = args.genes

    assert os.path.exists(__input_path), "Input path cannot be found"
    assert os.path.exists(__species_path), "Species path cannot be found"
    assert os.path.exists(__genes_path), "Genes path cannot be found"

    concat_fastas(__input_path, __output_path, __species_path, __genes_path)
