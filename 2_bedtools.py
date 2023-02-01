import sys
import os
import subprocess

cdir = os.getcwd()
gdir = "/drive/gff/"
fdir = "/drive/fasta/"
ndir = cdir + "/nucl/"

for root, dirs, files in os.walk(fdir):
    for file_name in files:
        #use only fasta files
        if file_name.endswith(".fasta"):
            name = file_name.split('.')[0]
            path_to_gff3 = gdir + name + ".gff3" #file name is with .gff3 extension
            path_to_gff = gdir + name + ".gff"
            path_to_fasta = fdir + file_name
            path_to_n_fasta = ndir + name + "_nucl.fasta"

            #check if file exists as gff3
            if os.path.exists(path_to_gff3):
                if not os.path.exists(path_to_n_fasta):
                    #print(['bedtools', 'getfasta', "-fi", path_to_fasta, '-bed', path_to_gff3, '-fo', path_to_n_fasta])
                    subprocess.run(['bedtools', 'getfasta', "-fi", path_to_fasta, '-bed', path_to_gff3, '-fo', path_to_n_fasta])
            else:
                if os.path.exists(path_to_gff):
                    print("No GFF3 or GFF file for " + name)
                else:
                    print("GFF3 file for " + name + " does not exist! Use GFF file!")
            
