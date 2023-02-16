"""Writes a list of species based on output files.
For this it uses a specific regex only working in my case,
you should probably make a new regex for your case.

My case is like this *_out_{ID}_gff_*. It takes the id
and writes it in a file
"""
import os
import re
import argparse

PARSER = re.compile(".*_out_(.*)_gff_.*")


def write_species_list(species_path: str):
    """Writes a list of species (species.txt) based on output files plus a regex matching the IDs

    Args:
        species_path (str): path to dir containing the species specific files
    """
    result = os.listdir(species_path)

    with open("species.txt", "w", encoding="utf-8") as file:
        for i in result:
            file.write(PARSER.match(i).group(1) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        help="Input path to folder containing specific species files",
    )

    args = parser.parse_args("-i", "hmm_output5/")

    if args.input is None:
        print("Input path not specified")
        exit()

    __input_path = args.input

    assert os.path.exists(__input_path), "Input path could not be found"

    write_species_list(__input_path)
