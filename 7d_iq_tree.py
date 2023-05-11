import os
import subprocess

path = "/processing/jgi/output_fastas_per_gene_v3/"

files = os.listdir(path)

for file in files:
    if file.endswith(".fasta"):
        oldfile = path + file
        newfile = file.split(".fasta")[0].split("trimal_")[1]
        command = "iqtree2 -s " + oldfile + " -m MF -T AUTO --prefix iq_" + newfile
        subprocess.run(command, shell=True, check=False)
