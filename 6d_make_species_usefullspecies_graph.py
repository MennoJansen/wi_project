# This script determines in how many species specific genes are present and occur once per species
# Input: HMMsearch output
# Output: tsv file with the number of selected gene and in how many species this gene is present

# Import packages
import gzip
import os

# Relevant directories
input_path = 'hmm_output2/'
output_path = 'hmm_intersection/'
result = os.listdir(input_path)

# Variables
lower_bound = 5
upper_bound = 10

# Open files
long_list = open(output_path + 'all_presence_once_good.txt').readlines()
unique_long_list = list(set(long_list))
output_file = open(output_path + 'accurate_bigger_table1_' + str(lower_bound) + '_' + str(upper_bound) + '.txt', 'w+')

def count_genes(long_list, unique_long_list, threshold):
    intersection_result = []
    for gene in unique_long_list:
        if long_list.count(gene) >= threshold:
            intersection_result.append(gene.strip())
    return intersection_result

def count_species(file, gene_list, counter):
    with open(input_path + file, 'rt') as f:

        # genes_in_species is a list that contains the genes that are present in this species
        # incl duplicates
        genes_in_species = []

        for line in f:
            if not line.startswith('#'):
                gene_name = line.split('-')[1].strip()
                genes_in_species.append(gene_name.strip())
        
        genes_in_species_without_duplicates = list(set(genes_in_species))

        # check if all 928 genes are in this hmmsearch output file of this specific species
        result = all(elem in genes_in_species_without_duplicates for elem in gene_list)
        print(gene_list)
        print(genes_in_species_without_duplicates)
        if result:
            counter += 1
        return counter

# Save filenames from hmmsearch output in files
files = []
for i in result:
    # Four files are excluded
   # if i != 'predicted_genes_aa_Aspoch1.txt.gz' and i != 'predicted_genes_aa_Aspnom1.txt.gz' and \
            #i != 'predicted_genes_aa_Xylcur1.txt.gz' and i != 'predicted_genes_aa_Aspcal1.txt.gz':
    files.append(i)

# Create tsv file with the number of genes and the corresponding occurrence in species
for i in range(lower_bound, upper_bound, 1):
    print("Counting for lowerbound of " + str(i))
    counter = 0
    genes_above_threshold = count_genes(long_list, unique_long_list, i)
    for file in files:
        counter = count_species(file, genes_above_threshold, counter)
    output_file.write(str(i) + ' ' + str(len(genes_above_threshold)) + '  ' + str(counter) + '\n')

output_file.close()
