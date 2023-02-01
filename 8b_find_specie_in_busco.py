import sys
import os

busco_path = "/processing/jgi/output_fastas_per_gene_v3/"
hmm_output_path = "/processing/jgi/hmm_output2/"
buscos_list_file = "buscos.txt"
buscos_len = 60

name = input("Enter name of species: ")
counter = 0 
counter_2 = 0
hmm_buscos = []
file_buscos = []


with open(hmm_output_path + "hmm_out_" + name + "_gff_proteins.txt") as hmm_file:
    for line in hmm_file:
        with open(buscos_list_file) as blf:
            for busco in blf:
                busco = busco.strip()               
                if busco in line:
                    hmm_buscos.append(busco)
                    counter = counter + 1
                    print("Busco " + busco.strip() + " is present in species hmm output")

with open(buscos_list_file) as blf:
    for busco in blf:
        with open(busco_path + busco.strip() + ".txt") as busco_file:
            for line in busco_file:
                if name.strip() in line:
                    file_buscos.append(busco.strip())
                    counter_2 = counter_2 + 1
                    print("Name of specie found in busco file: " + busco)

print(set(hmm_buscos) - set(file_buscos))
print("Busco present in (" + str(counter) + "/" + str(buscos_len) + ")")
print("Name present in (" + str(counter_2) + "/" + str(buscos_len) + ")")