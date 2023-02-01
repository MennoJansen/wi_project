import os
import sys

input_fastas = "/processing/jgi/output_fastas_per_gene_v3/"
output_fasta = "/processing/jgi/output_fastas_per_gene_v3/cat_trimal.fasta"

species_file = "/processing/jgi/species.txt"
buscos = "/processing/jgi/buscos.txt"

fasta_files = os.listdir(input_fastas)

new_fasta = dict()
with open(species_file) as all_species:
    for species in all_species:
        species = species.strip()
        new_fasta[species] = ''

with open(buscos) as genes:
    total_length = 0
    for gene in genes:
        # First, scan gene fasta file, save sequences in new_fasta dict
        gene_fasta_file = input_fastas + "trimal_" + gene.strip() + ".fasta"
        len_name = 'undefined'
        print("Running for gene " + gene.strip())
        with open(gene_fasta_file) as gene_fastas:
            name = 'undefined'
            for line in gene_fastas:
                # If line found is a specie name, set in var name
                # Else, we are in a sequence of a specie, so we add this sequence
                # to the existing sequence of this specie
                if line.startswith('>'):
                    name = line.split('>')[1].strip() # >Aspcost1\n becomes Aspcost1
                    len_name = name
                else:
                    new_fasta[name] = new_fasta[name] + line.strip()

        # Second, add empty sequences for species where the BUSCO gene is not present
        seq_length = len(new_fasta[len_name]) - total_length
        total_length = len(new_fasta[len_name])
        print(total_length)
        for name, seq in new_fasta.items():
            if len(seq) != total_length:
                empty_seq = "-" * seq_length
                new_fasta[name] = new_fasta[name] + empty_seq
            #new_fasta[name] = new_fasta[name] + 'P'

with open(output_fasta, 'w') as of:
    for name, seq in new_fasta.items():
        of.write(">" + name + "\n")
        of.write(seq + "\n")
