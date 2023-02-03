"""Processes a .fasta file so the sequences are in one line per entry

Parameters
----------
input : str
    Input folder containing .fasta files needing to be processed
"""

import argparse
import os

def process_fasta(i_file : str):
    """Method processing .fasta file to move sequences into one line

    Parameters
    ----------
    i_file : str
        .fasta file that needs to be processed
    """
    lines = []
    with open(i_file, encoding="utf-8") as fin:
        seq = ''
        for line in fin.readlines():
            if line.startswith('>'):
                if seq != '':
                    lines.append(f"{seq}\n")
                lines.append(line)
                seq = ''
            else:
                seq+= line.strip()

    with open(i_file, 'w', encoding="utf-8") as fout:
        fout.writelines(lines)

#process_fasta("/drive/fasta/Albra1.fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", 
                        help="Input directory containing .fasta files needing to be processed")

    args = parser.parse_args()

    if args.input is None:
        print("Input directory not specified...")
        exit()
    input_dir = args.input

    assert os.path.exists(input_dir), "input_dir could not be found"

    for file in os.listdir(input_dir):
        process_fasta(f"{input_dir}{file}")
