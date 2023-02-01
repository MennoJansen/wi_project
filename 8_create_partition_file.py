# This code creates the partition file which can be used by IQtree to build a tree based on multiple trees
# Input: .txt file with species names, trimmed fasta files per gene
# Output: new fasta file with all genes and all species + partition file

# Import packages
import os
import gzip

# Open directories
trimal_files_path = 'output_fastas_per_gene_v3/'
output_path = 'output_partition_files_v3/'


#new_fasta = open(output_path + 'new_fasta.fasta', 'w+')

species_file = open('species.txt')
genes_file = 'buscos.txt'

data = os.listdir(trimal_files_path)

# Create partition file
partition_file = open(output_path + 'partition_file.nex', 'w+')
partition_file.write('#nexus\n')
partition_file.write('begin sets;\n')

# Variables
gene_names = []
big_fasta_dict = dict()

# Create dict with species as keys without values
for species in species_file:
    big_fasta_dict[species.strip()] = ''

# Get gene names
for i in data:
    if ".fasta.gz" in i:
        gene_names.append(i.split('.')[0].split('_')[1])

start = 1
counter = 0
charpartition_mine = '  charpartition mine = '


with open(genes_file) as gf:
    for gene_name in gf:
        gene_name = gene_name.strip()
        partition_file.write('  charset ' + gene_name + ' = trim_' + gene_name + '.phy: *'  + ';' + '\n')
        with open("iq_" + gene_name + ".iqtree") as IQtree_output:
            for line in IQtree_output:
                if 'Best-fit model according to BIC:' in line:
                    best_model_BIC = line.split(': ')[1]
                    charpartition_mine += str(best_model_BIC.strip()) + ':' + gene_name + ', '

'''
# Concatenate sequence of all genes per species
for gene_name in gene_names:
    counter += 1
    with gzip.open(trimal_files_path + 'trimal_' + gene_name + '.fasta', 'rt') as trimal_of_gene:
        for line in trimal_of_gene:
            if line.startswith('>'):
                trimal_species = line.split('>')[1].strip()
            else:
                current_sequence = big_fasta_dict[trimal_species]
                new_sequence = current_sequence + line.strip()
                big_fasta_dict[trimal_species] = new_sequence

        length_of_gene = len(big_fasta_dict[trimal_species])
        end = length_of_gene
        partition_file.write('  charset ' + gene_name + ' = ' + str(start) + '-' + str(end) + ';' + '\n')
        start = end + 1
'''
'''    IQtree_output = open(IQtree_files_path + gene_name + '.iqtree')
    for line in IQtree_output:
        if 'Best-fit model according to BIC:' in line:
            best_model_BIC = line.split(': ')[1]
            charpartition_mine += str(best_model_BIC.strip()) + ':' + gene_name + ', '
'''
# Finish partition file
# Remove the last ', ' of charpartition
charpartition_mine = charpartition_mine[:-2]

# Add ';' to charpartition
charpartition_mine += ';'

partition_file.write(charpartition_mine + '\n')
partition_file.write('end;')

'''
# Save new fasta file
for species in species_file:
    species = species.strip()
    header = '>' + species
    sequence = big_fasta_dict[species]
'''