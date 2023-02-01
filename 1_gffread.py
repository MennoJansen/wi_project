"""Creates files containing protein sequences for all genes for different species

Nucleotide .fasta files and .gff3 files should have the same names, protein .fasta
will receive the same name:
Abc123.fasta (n) + Abc123.gff3 -> Abc123.fasta (p)

gdir : str
    Location of the .gff3 files
fdir : str
    Location of the nucleotide fasta sequences in .fasta format
pdir : str
    New location of resulting protein sequences in .fasta format
"""

import os
import subprocess
import argparse

def gffread(gdir : str, fdir : str, pdir : str):
    """Gets the protein sequence using nucleotide sequence and .gff3 file, writes to .fasta file

    This uses the program gffread, tested on version v0.12.7
    More info on gffread: https://github.com/gpertea/gffread

    gdir : str
        Location of the .gff3 files
    fdir : str
        Location of the nucleotide fasta sequences in .fasta format
    pdir : str
        New location of resulting protein sequences in .fasta format
    
    """
    for root, dirs, files in os.walk(gdir):
        for file_name in files:
            if file_name.endswith(".gff3"):
                name = file_name.split('.')[0]
                path_to_gff = gdir + file_name 
                path_to_fasta = fdir + name + ".fasta"
                path_to_proteins_fasta = pdir + name + "_gff_proteins.fasta"

                #check if file exists as fasta)
                if os.path.exists(path_to_fasta):
                    #check if not already done
                    if not os.path.exists(path_to_proteins_fasta):
                        subprocess.run(['gffread', '-y', path_to_proteins_fasta, '-g', path_to_fasta, path_to_gff])
                else:
                    print("Fasta file for " + name + " does not exist!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-g", "--gffdir", dest = "gdir", help="Dir containing .gff3 files")
    parser.add_argument("-f", "--fastadir", dest = "fdir", help="Dir containing .fasta nucleotide sequences")
    parser.add_argument("-p", "--proteindir", dest = "pdir", help="Output dir for protein sequences (.fasta)")

    args = parser.parse_args()

    cdir = os.getcwd() + "/"

    if args.gdir is None:
        print("gdir not specified")
        exit()
    if args.fdir is None:
        print("fdir not specified")
        exit()
    if args.pdir is None:
        print("pdir not specified")
        exit()

    gdir = cdir + args.gdir
    fdir = cdir + args.fdir
    pdir = cdir + args.pdir
    if not gdir.endswith("/"):
        gdir += "/"
    if not fdir.endswith("/"):
        fdir += "/"
    if not pdir.endswith("/"):
        pdir += "/"
    print(f"Location of .gff files: {gdir}")
    print(f"Location of .fasta files (nucleotide sequences): {fdir}")
    print(f"Location for .fasta files (protein sequences): {pdir}")
    if input("Are locations okay? [y/n]: ") == "y":
        print("Running program...")
        gffread(gdir, fdir, pdir)
    print("Exiting program...")
    exit()

