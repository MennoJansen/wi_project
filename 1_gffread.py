import sys
import os
import subprocess

cdir = os.getcwd()
#gdir = cdir + "/gff/"
#fdir = cdir + "/fasta/"
gdir = "/drive/gff/"
fdir = "/drive/fasta/"
pdir = cdir + "/proteins/"

for root, dirs, files in os.walk(gdir):
    for file_name in files:
        #use only gff3 files
        if file_name.endswith(".gff3"):
            name = file_name.split('.')[0]
            path_to_gff = gdir + file_name #file name is with .gff3 extension
            path_to_fasta = fdir + name + ".fasta"
            path_to_proteins_fasta = pdir + name + "_gff_proteins.fasta"

            #check if file exists as fasta)
            if os.path.exists(path_to_fasta):
                #check if not already done
                if not os.path.exists(path_to_proteins_fasta):
                    print(['gffread', '-y', path_to_proteins_fasta, '-g', path_to_fasta, path_to_gff])
                    subprocess.run(['gffread', '-y', path_to_proteins_fasta, '-g', path_to_fasta, path_to_gff])
            else:
                print("Fasta file for " + name + " does not exist!")
