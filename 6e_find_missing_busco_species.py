import sys
import os
import re

input = "hmm_intersection/buscos.txt"
input_hmms = "hmm_output2/"
output = "hmm_intersection/missing.txt"

p = re.compile('.*_out_(.*)_gff_.*')

# Import buscos from a txt file
# Returns list of buscos from file
def import_buscos(busco_file:str):
    buscos = []
    with open(busco_file) as bfile:
        for line in bfile.readlines():
            busco = line.strip()
            buscos.append(busco)
    return buscos

# Find buscos from list of buscos in hmm output file
# Returns list of buscos which are not found in hmm output file
def find_buscos(hmm_output:str, buscos:list):
    not_found_buscos = list(buscos)
    with open(hmm_output) as hfile:
        for line in hfile.readlines():
            if not line.startswith('#'):
                busco = line.split('-')[1].strip()
                if busco in not_found_buscos:
                    not_found_buscos.remove(busco)
    return not_found_buscos

# Write not found buscos to file
def write_to_file(not_found_buscos:list, specie:str, output:str):
    with open(output, 'a') as ofile:
        ofile.write(specie + " has to following missing BUSCOs: " + str(not_found_buscos) + '\n')

buscos = import_buscos(input)

files = os.listdir(input_hmms)
for file in files:
    print("running for " + file)
    name = p.match(file).group(1)
    not_found_buscos = find_buscos(input_hmms + file, buscos)
    if len(not_found_buscos) > 0:
        write_to_file(not_found_buscos, name, output)