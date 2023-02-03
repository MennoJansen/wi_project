"""Small script that removes duplicate lines from a file

Parameters
----------
file : str
    Path to file with duplicate lines

"""
import argparse
import os
from collections import OrderedDict

def remove_duplicate_lines(i_file : str):
    """Removes duplicate lines from file and overwrites same file

    Parameters
    ----------
    file : str
        Filepath of file with duplicate lines
    """
    with open(i_file, encoding="utf-8") as fin:
        lines = (line.rstrip() for line in fin)
        unique_lines = OrderedDict.fromkeys((line for line in lines if line))

    with open(i_file, 'w', encoding="utf-8") as out:
        for key in unique_lines.keys():
            out.write(key + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", dest= "file", help= "File with duplicate lines")

    args = parser.parse_args()

    if args.file is None:
        print("file not specified")
        exit()

    file = args.file
    assert os.path.exists(file), "file could not be found"

    print(f"Removing duplicate lines from: {file}")
    remove_duplicate_lines(file)
    print("Done! Exiting...")
    