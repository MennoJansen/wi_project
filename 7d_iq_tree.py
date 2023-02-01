import os
import sys
import subprocess

path = "/processing/jgi/output_fastas_per_gene_v3/"

files = os.listdir(path)

for file in files:
    if file.endswith(".phy"):
        oldfile = path + file
        newfile = file.split(".phy")[0].split("trim_")[1]
        command = "iqtree2 -s " + oldfile + " -m MF -T 6 --prefix iq_" + newfile
        subprocess.run(command, shell=True)