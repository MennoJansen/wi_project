"""Extracts .gff.gz and .gff3.gz files from a zip file if it is missing.
Missing counts as an id where a .fasta file is found but a .gff or .gff3
is not found

Parameters
----------
f_dir : str
    Directory containing the .fasta files
g_dir : str
    Directory containing the .gff and .gff3 files
g_zip : str
    .zip file containing .gff and .gff3 files

"""

import os
import argparse
import zipfile
import gzip
import shutil


def get_file_from_zip(name: str, g_zip: str, g_dir: str):
    """Method for getting .gff.gz or .gff3.gz from .zip file
    Also unpacks and renames the file to [id].gff or [id].gff3

    Parameters
    ----------
    name : str
        Name of specie
    g_zip : str
        .zip file containing .gff and .gff3 files
    g_dir : str
        Directory for .gff and .gff3 files
    """
    with zipfile.ZipFile(g_zip) as z_file:
        for z_name in z_file.namelist():
            z_dir, file = z_name.split("/")
            split_file = file.split(".")
            extension = f".{split_file[1]}.{split_file[2]}"
            # print(f"{z_dir} = {name}")
            if z_dir == name:
                if (
                    "_GeneCatalog_" in file
                    or "filtered_proteins" in file
                    or "filtered_genes"
                    or "FilteredModels" in file
                ) and ".gff3" in file:
                    print(f"Found file [{file}], extracting it to {g_dir}")
                    g_file = f"{g_dir}{name}{extension}"

                    with open(g_file, "wb") as n_file:
                        n_file.write(z_file.read(z_name))

                    with gzip.open(g_file, "rb") as f_in:
                        if "filtered_proteins" in file:
                            ext = split_file[3]
                        elif "filtered_genes" in file:
                            ext = split_file[2]
                        else:
                            ext = split_file[1]
                        with open(f"{g_dir}{name}.{ext}", "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(g_file)
                    if ".gff3" in file:
                        break


def check_missing_files(f_dir: str, g_dir: str, g_zip: str):
    """Method checking for missing .gff or .gff3 files
    Calls get_file_from_zip to retrieve the missing files

    Parameters
    ----------
    f_dir : str
        Directory containing .fasta files
    g_dir : str
        Directory containing .gff and .gff3 files
    g_zip : str
        .zip file containing .gff and .gff3 files
    """
    counter = 0
    for file in os.listdir(f_dir):
        if file.endswith(".fasta"):
            name = file.split(".")[0]
            if not os.path.exists(f"{g_dir}{name}.gff3"):
                print(f"No gff file found for {file}")
                print("Trying to retrieve file from zip")
                get_file_from_zip(name, g_zip, g_dir)
                counter += 1
    print(f"{counter} missing files!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--fdir", dest="f_dir", help="Directory containing .fasta files"
    )
    parser.add_argument(
        "-g", "--gdir", dest="g_dir", help="Directory containing .gff and .gff3 files"
    )
    parser.add_argument(
        "-z", "--zip", dest="g_zip", help="Zip file containing .gff and .gff3 files"
    )

    args = parser.parse_args(
        [
            "-f",
            "/processing/jgi/nucl/",
            "-g",
            "/drive/gff/",
            "-z",
            "/home/menno/Downloads/annotation.zip",
        ]
    )

    if args.f_dir is None:
        print("f_dir not specified")
        exit()
    if args.g_dir is None:
        print("g_dir not specified")
        exit()
    if args.g_zip is None:
        print("g_zip not specified")
        exit()
    __f_dir = args.f_dir
    __g_dir = args.g_dir
    __g_zip = args.g_zip

    assert os.path.exists(__f_dir), "f_dir could not be found"
    assert os.path.exists(__g_dir), "g_dir could not be found"
    assert os.path.exists(__g_zip), "g_zip could not be found"
    assert __g_zip.endswith(".zip"), "g_zip is not a .zip file"

    check_missing_files(__f_dir, __g_dir, __g_zip)
