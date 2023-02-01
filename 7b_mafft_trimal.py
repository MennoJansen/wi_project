import os
import sys
import subprocess

path = "/processing/jgi/output_fastas_per_gene_v3/"

files = os.listdir(path)

for file in files:
    if file.endswith(".txt"):
        oldfile = path + file
        newfile1 = path + "mafft_" + str(file.split(".txt")[0]) + ".fasta"
        newfile2 = path + "trimal_" + str(file.split(".txt")[0]) + ".fasta"
        print("Running for " + oldfile)
        if not os.path.exists(newfile1):
            command_1 = "mafft --preservecase --auto " + oldfile + " > " + newfile1
            subprocess.run(command_1, shell=True)
        if not os.path.exists(newfile2):
            command_2 = "trimal -in " + newfile1 + " -out " + newfile2 + " -automated1"            
            subprocess.run(command_2, shell=True)
        #subprocess.run(["trimal", "-in", "mafft" + file, "-out", "trim_" + str(file.split(".txt")) + ".phy", "-automated1", "-phylip"])
        #QC trimal -seqoverlap 80 -resoverlap 0.75