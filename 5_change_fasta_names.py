# The downloaded JGI fasta and GFF files do have additional information
# besides the species tree. This script changes the file names into
# the species name + extesnion
# Input: Fasta folder
# Output: modified file names

# Import libraries
import pathlib
import os

# Change fasta files
def change_fasta_name(file_name):
    if '_AssemblyScaffolds.fasta.gz' in file_name:
        species = file_name.split('_AssemblyScaffolds.fasta.gz')[0]
    elif '_assembly_scaffolds.fasta.gz' in file_name:
        species = file_name.split('_assembly_scaffolds.fasta.gz')[0]
    elif '.AssembledScaffolds.fasta.gz' in file_name:
        species = file_name.split('.AssembledScaffolds.fasta.gz')[0]
    elif '_scaffolds.fasta.gz' in file_name:
        species = file_name.split('_scaffolds.fasta.gz')[0]
    elif '.nuclearAssembly.unmasked.fasta.gz' in file_name:
        species = file_name.split('.nuclearAssembly.unmasked.fasta.gz')[0]
    elif '.AssemblyScaffolds.fasta.gz' in file_name:
        species = file_name.split('.AssemblyScaffolds.fasta.gz')[0]
    elif '.unmasked_nuclear_assembly.fasta.gz' in file_name:
        species = file_name.split('.unmasked_nuclear_assembly.fasta.gz')[0]
    elif '_unmasked_assembly.fasta.gz' in file_name:
        species = file_name.split('_unmasked_assembly.fasta.gz')[0]
    else:
        print(file_name)
        species = 'X'
    old_name = path + file_name
    new_name = path + species + '.fasta.gz'
    os.rename(old_name, new_name)

# Get the right path
path = str(pathlib.Path(__file__).parent.resolve()) + '/'
result = os.listdir(path)

# Save fasta files in list
files = []
for i in result:
    files.append(i)

# Change fasta file names
for file in files:
    if file.endswith('.gz'):
        change_fasta_name(file)
