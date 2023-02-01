# This script creates a fasta file per gene
# It retrieves the gene ID from the HMMsearch output
# This gene ID is used to retrieve the location info (scaffold, start, end) from GFF files
# This location info is used to extract the corresponding sequence
# Input: .txt file with genes you want to do, HMMsearch output, GFF files and fasta files
# Output: Fasta file per gene

# Import packages
import gzip
import os
import re
import logging
from timeit import default_timer as timer

# Open relevant directories
gff_path = '/drive/gff/'
g_fasta_path = '/processing/jgi/nucl/'
hmmsearch_path = '/processing/jgi/hmm_output2/'
output_path = '/processing/jgi/output_fastas_per_gene_v3/'

gff_files = os.listdir(gff_path)
g_fasta_files = os.listdir(g_fasta_path)
hmmsearch_files = os.listdir(hmmsearch_path)

#Regex for hmm output files
p = re.compile('.*_out_(.*)_gff_.*')

# Open .txt files with genes you need the sequence from
#genes = open('/processing/jgi/hmm_intersection/buscos.txt').readlines()
genes = open('/processing/jgi/buscos_left.txt').readlines()
counter = 1
counter_end = len(genes)

logging.basicConfig(filename='7_fasta_per_gene.log', encoding='utf-8', level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.getLogger().addHandler(logging.StreamHandler())

# Open species you need the sequence from
species = open('species.txt').readlines()
s_counter = 1
s_counter_end = len(species)
species_list = []
for i in species:
    species_list.append(i.strip())

def get_gene_name(hmmsearch_output, gene):
    target_name = "not_found"
    for line in hmmsearch_output:
        if not line.startswith('#'):
            split_line = line.split('-')
            query_name = split_line[1].strip()
            if query_name == gene:
                target_name = split_line[0].strip()
                id_mrna = f'ID={target_name};'
                return id_mrna

    '''if 'mRNA_' in target_name:
        gene_ID = target_name.split('RNA_')[1]
    elif 'mRNA' in target_name and 'mRNA_' not in target_name:
        gene_ID = target_name.split('RNA')[1]
    else:
        #print('ERROR: OTHER TARGET NAME: ', target_name)
        return target_name

    gene_name = 'ID=gene_' + gene_ID + ';'
    '''
    return target_name
    

def get_location(gff_output, gene_name):
    scaffold, start, end, orientation = '0', '0', '0', '0'
    for line in gff_output:
        split_gene_name = ''.join(gene_name.split('_'))
        if gene_name in line or split_gene_name in line:
            scaffold = line.split()[0]
            start = line.split()[3]
            end = line.split()[4]
            orientation = line.split()[6]
            break

    if int(end) < int(start):
        new_end = start
        new_start = end
        start = new_start
        end = new_end
    if gene_name != "not_found" and start == '0':
        logging.warning("WARNING!")
    start = str(int(start) - 1)

    fasta_header = '>' + scaffold + ':' + start + '-' + end
    #print(fasta_header)
    return fasta_header, orientation

def get_sequence(g_fasta_output, fasta_header, output_file, orientation, species):
    for i in range(0, len(g_fasta_output) - 1, 2):
        if g_fasta_output[i].strip() == fasta_header:
            output_file.write('>' + species + '\n')

            # Take the complement sequence
            # + 1 in indexing is needed due to a +1 shift in counting
            if orientation == '+':
                output_file.write(g_fasta_output[i + 1].upper())
                return
            elif orientation == '-':
                complement = ''
                tab = str.maketrans('ACTG', 'TGAC')
                output_file.write(g_fasta_output[i + 1].upper().translate(tab).strip()[::-1] + '\n')
                '''
                seq = g_fasta_output[i + 1].strip()
                seq = seq.replace('A', 't').replace('C', 'g').replace('T', 'a').replace('G', 'c')
                seq = seq.upper()
                seq = seq[::-1]

                output_file.write(seq + '\n')
                '''
                '''
                for nucleotide in g_fasta_output[i + 1]:
                    if nucleotide == 'A':
                        complement += 'T'
                    elif nucleotide == 'C':
                        complement += 'G'
                    elif nucleotide == 'T':
                        complement += 'A'
                    elif nucleotide == 'G':
                        complement += 'C'
                reverse_complement = complement[::-1]
                output_file.write(reverse_complement.strip() + '\n')
                '''
                return
            else:
                print('other orientation than "+" and "-"')
    

def create_fasta_per_gene(gene):
    global counter, s_counter
    logging.info("Running for gene: %s (%s/%s)", gene, counter, counter_end)
    counter = counter + 1
    s_counter = 1
    # Get gene ID
    for hmmsearch_file in hmmsearch_files:
        with open(hmmsearch_path + hmmsearch_file, 'rt') as hmmsearch_output:
            #species = hmmsearch_file.split('_genes_aa_')[1].split('.txt')[0]
            species = p.match(hmmsearch_file).group(1)
            '''
            if "Clasph" in species:
                print("Found it!")
                print("breakpoint")
            else:
                continue
            '''
            
            #logging.INFO("\rRunning for species: " + species + " (" + str(s_counter) + "/" + str(s_counter_end) + ")                ", end='')
            logging.info("Running for species %s (%s/%s)", species, s_counter, s_counter_end)
            s_counter+=1
            if species in species_list:
                # Get gene name
                gene_name = get_gene_name(hmmsearch_output, gene)
                if gene_name == "not_found":
                    output_file = open(output_path + gene + '.txt', 'a')
                    output_file.write('>' + species + '\n\n')
                    output_file.close()

                # Retrieve location info of gene
                with open(gff_path + species + ".gff3", 'rt') as gff_output:
                    get_loc = get_location(gff_output, gene_name)
                    fasta_header = get_loc[0]
                    orientation = get_loc[1]

                    # Extract sequence
                    with open(g_fasta_path + species + "_nucl.fasta", 'rt') as g_fasta_output :
                        g_fasta_output = g_fasta_output.readlines()                                      
                        output_file = open(output_path + gene + '.txt', 'a')
                        get_sequence(g_fasta_output, fasta_header, output_file, orientation, species)
                        output_file.close()
#create_fasta_per_gene("408391at4751")
                                        
for gene in genes:
    gene = gene.strip()
    create_fasta_per_gene(gene)
