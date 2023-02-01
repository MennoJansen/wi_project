import sys
import os
from collections import OrderedDict

file = "/drive/gff/Exoaq1.gff3"
with open(file) as fin:
    lines = (line.rstrip() for line in fin)
    unique_lines = OrderedDict.fromkeys((line for line in lines if line))

#print(unique_lines.keys())
with open(file, 'w') as out:
    for key in unique_lines.keys():
        #print(key)
        out.write(key + '\n')
