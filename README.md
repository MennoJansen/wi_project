# wi_project

Common order of bussiness
-------------------------
1. Download annotation and assembly from JGI
2. Unzip annotation to gff folder, only unzipping .gff3 files (use extract_gff_from_zip.py)
3.1 Unzip assembly to fasta folder, unzipping only the files, no dirs (use unzip -j "assembly.zip")
3.2 Rename files to {id}.fasta (remove mito scaffolds, some files need manual renaming, others can be renamed using mmv)
4. Use gffread to produce protein fasta files (use gff_read.py)
5. Use hmmsearch on the protein fasta files with a fungal model
6. Select busco genes you want to use
7. Flatten .fasta files so they become single line sequences (use process_fasta_file.py)
8. Get sequences matching the busco genes for all species (use 7_fasta_per_gene.py)
9. Align and trim all new busco gene files (use mafft_trimal.py)
10. Concat all genes .fasta files (use concat_fasta.py)
11. Get classes for all species (use add_class_name.py)
12. Create annotation file (use create_annotation_file.py)