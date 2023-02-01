import sys
import os
import fileinput

path = "/processing/jgi/output_fastas_per_gene_v3/"

files = os.listdir(path)

for file in files:
    if file.endswith(".txt"):
        for line in fileinput.input(path + file, inplace=True):
            if line.strip() == '':
                print("XXXXXXXXXXXXXXXXXXXXXXXXXXX\n", end='')
            else:
                print(line, end='')
