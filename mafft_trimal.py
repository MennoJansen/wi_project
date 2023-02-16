"""Uses output from script process_fasta_file, generates aligned, trimmed .fasta files for each of the BUSCO genes
It requires mafft and trimal to be installed on the system
"""
import os
import subprocess
import argparse


def mafft_trimal(path: str):
    """Uses mafft and trimal to generate aligned and trimmed .fasta files

    Args:
        path (str): input path to dir containing input .txt files (in fasta format)
    """
    files = os.listdir(path)

    for file in files:
        if file.endswith(".txt"):
            oldfile = path + file
            gene_name = str(file.split(".txt")[0])
            new_file_mafft = f"{path}mafft_{gene_name}.fasta"
            new_file_trimal = f"{path}trimal_{gene_name}.fasta"
            print("Running for " + oldfile)
            if not os.path.exists(new_file_mafft):
                command_1 = (
                    "mafft --preservecase --auto " + oldfile + " > " + new_file_mafft
                )
                subprocess.run(command_1, shell=True, check=False)
            if not os.path.exists(new_file_trimal):
                command_2 = (
                    "trimal -in "
                    + new_file_mafft
                    + " -out "
                    + new_file_trimal
                    + " -automated1"
                )
                subprocess.run(command_2, shell=True, check=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="input_path",
        help="Input path containing all the files resulting from process_fasta_file.py",
    )

    args = parser.parse_args(["-i", "/processing/jgi/output_fastas_per_gene_v3/"])

    if args.input_path is None:
        print("No input path specified")
        exit()

    input_path = args.input_path

    assert os.path.exists(input_path), "Input path could not be found"

    mafft_trimal(input_path)
