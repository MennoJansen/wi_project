import os
import sys

path = "/processing/jgi/output_fastas_per_gene_v3/"

files = os.listdir(path)

for file in files:
    if file.endswith(".txt"):
        with open(path + file) as ofile:
            for line in ofile.readlines():
                if not line.startswith('>') and 'n' in line:
                    print("Lowercase in " + file)
                    print("Line :" + line)