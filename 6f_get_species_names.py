import os
import sys
import re

p = re.compile('.*_out_(.*)_gff_.*')

species_path = "hmm_output2/"
result = os.listdir(species_path)

with open("species.txt", 'w') as file:
    for i in result:
        file.write(p.match(i).group(1) + '\n')
