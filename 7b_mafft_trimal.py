import os
import sys
import subprocess

path = "/processing/jgi/output_fastas_per_gene_v3/"

files = os.listdir(path)

for file in files:
    if file.endswith(".txt"):
        oldfile = path + file
        gene_name = str(file.split(".txt")[0])
        new_file_mafft = f"{path}mafft_{gene_name}.fasta"
        new_file_trimal = f"{path}trimal_{gene_name}.fasta"
        print("Running for " + oldfile)
        if not os.path.exists(new_file_mafft):
            command_1 = (
                "mafft --preservecase --auto " + oldfile + " > " + new_file_mafft
            )
            subprocess.run(command_1, shell=True)
        if not os.path.exists(new_file_trimal):
            command_2 = (
                "trimal -in "
                + new_file_mafft
                + " -out "
                + new_file_trimal
                + " -automated1"
            )
            subprocess.run(command_2, shell=True)
        # subprocess.run(["trimal", "-in", "mafft" + file, "-out", "trim_" + str(file.split(".txt")) + ".phy", "-automated1", "-phylip"])
        # QC trimal -seqoverlap 80 -resoverlap 0.75
