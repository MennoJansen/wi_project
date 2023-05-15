"""Removes duplicate files with different extensions
Input is directory, extension, duplicate extension
Example, removes .gff files if a .gff3 file is found with the same name
"""
import os
import argparse


def remove_duplicate_files_diff_ext(i_folder: str, ext1: str, ext2: str):
    """Removes duplicate files with different file extensions
    For example, remove .gff files if a .gff3 file is found with the same name

    Args:
        i_folder (str): folder containing duplicate files
        ext1 (str): extension of files to be kept
        ext2 (str): extension of files to be removed
    """
    for file in os.listdir(i_folder):
        if file.endswith(ext1):
            file_name = file.split(".")[0]
            dup_file = f"{i_folder}{file_name}{ext2}"
            if os.path.exists(dup_file):
                os.remove(dup_file)
                print(f"Deleted file: {dup_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", dest="input_folder", help="Folder containing duplicate files"
    )
    parser.add_argument(
        "-e",
        "--ext",
        dest="extension",
        help="Extension of files to be kept",
    )
    parser.add_argument(
        "-d",
        "--dup-ext",
        dest="dup_extension",
        help="Extension of duplicate files to be removed",
    )

    args = parser.parse_args()

    if args.input_folder is None:
        print("Folder not specified")
        exit()
    input_folder = args.input_folder
    assert os.path.exists(input_folder), "Input folder could not be found"
    if args.extension is None:
        print("Extension not specified")
        exit()
    extension = args.extension
    if not extension.startsWith("."):
        extension = f".{extension}"

    if args.dup_extension is None:
        print("Duplicate extension not specified")
        exit()
    dup_extension = args.dup_extension
    if not dup_extension.startsWith("."):
        dup_extension = f".{dup_extension}"

remove_duplicate_files_diff_ext(input_folder, extension, dup_extension)
