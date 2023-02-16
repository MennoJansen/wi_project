"""Given [x] protein sequences in [y] .fasta files, runs hmmsearch
over every file using a given model, resulting in [y] hmmsearch output files.

Input names should be formatted as follows:
ID_gff_proteins.fasta

Output names are formatted as follows:
hmm_out_ID_gff_proteins.txt

Parameters
----------
pdir : str
    Folder containing protein sequences in .fasta file
hmm_model : str
    Model that hmmer uses for hmmsearch
odir : str
    Folder for the outputs of the hmmsearches
"""
import os
import subprocess
import argparse

# check if file is ok
def check_file(file):
    """Checks if hmm_output file is finished by looking in the last line for '# [ok]'

    Parameters
    ----------
    file : str
        hmmer output file to scan

    Returns
    -------
    __return__ : bool
        True if hmmsearch was done, otherwise False
    """
    with open(file, "rb") as c_file:
        try:
            c_file.seek(-2, os.SEEK_END)
            while c_file.read(1) != b"\n":
                c_file.seek(-2, os.SEEK_CUR)
        except OSError:
            c_file.seek(0)
        last_line = c_file.readline().decode()
        if last_line.startswith("# [ok]"):
            return True
        return False


failed_hmms = []


def hmm_search(i_dir: str, hmm_file: str, o_dir: str):
    """Performs hmmsearch on all .fasta files in [pdir] using [hmm_file] as model

    Defaults values include threshold of 1e-10 for the E-value, 12 threads
    TODO Work default values into parameters

    Parameters
    ----------
    pdir : str
        Path to folder with protein .fasta files
    hmm_file : str
        File path to hmm model
    odir : str
        Path to output directory
    """
    total = 0
    counter = 1
    # Scan for how many files:
    for files in os.walk(i_dir):
        for file in files[2]:
            if file.endswith(".fasta"):
                total += 1
    for files in os.walk(i_dir):
        for file_name in files[2]:
            name = file_name.split(".")[0]
            print(f"hmmsearch for {file_name} ({counter}/{total})")
            counter += 1
            output_file = f"{o_dir}hmm_out_{name}.txt"
            if not os.path.exists(output_file) or not check_file(output_file):
                subprocess.run(
                    [
                        "hmmsearch",
                        "--cpu",
                        str(12),
                        "-E",
                        str(0.0000000001),
                        "--tblout",
                        output_file,
                        hmm_file,
                        i_dir + file_name,
                    ],
                    stdout=subprocess.DEVNULL,
                    check=False,
                )
            if not os.path.exists(o_dir + "hmm_out_" + name + ".txt"):
                failed_hmms.append(name)

    print("Hmmsearch failed for these entries:")
    print(failed_hmms)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="i_dir",
        help="Path to folder with sequences in .fasta files",
    )
    parser.add_argument("-m", "--model", dest="hmm_model", help="Filepath to hmm model")
    parser.add_argument("-o", "--out", dest="odir", help="Path to output directory")

    # args = parser.parse_args()
    args = parser.parse_args(
        [
            "-i",
            "/processing/jgi/proteins/",
            "-m",
            "fungal_models.hmm",
            "-o",
            "/processing/jgi/hmm_output5/",
        ]
    )
    cdir = os.getcwd()

    if args.i_dir is None:
        print("-i not specified")
        exit()
    if args.hmm_model is None:
        print("-m not specified")
        exit()
    if args.odir is None:
        print("-o not specified")
        exit()

    __i_dir = args.i_dir
    __hmm_model = args.hmm_model
    __odir = args.odir

    assert os.path.exists(__i_dir), "input could not be found"
    assert os.path.exists(__hmm_model), "hmm_model could not be found"
    assert os.path.exists(__odir), "odir could not be found"
    hmm_search(__i_dir, __hmm_model, __odir)
