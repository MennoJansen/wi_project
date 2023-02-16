# wi_project

Common order of bussiness
-------------------------
1. Download annotation and assembly from JGI
2. Unzip annotation to gff folder, only unzipping .gff3 files (use extract_gff_from_zip.py)
3. Unzip assembly to fasta folder, unzipping only the files, no dirs (use unzip -j "assembly.zip")
4. Rename files to {id}.fasta (remove mito scaffolds, some files need manual renaming, others can be renamed using mmv)
5. Use gffread to produce protein fasta files (use gff_read.py)
6. Use hmmsearch on the protein fasta files with a fungal model
7. Select busco genes you want to use
8. Flatten .fasta files so they become single line sequences (use process_fasta_file.py)
9. Get sequences matching the busco genes for all species (use 7_fasta_per_gene.py)
10. Align and trim all new busco gene files (use mafft_trimal.py)
11. Concat all genes .fasta files (use concat_fasta.py)
12. Get classes for all species (use add_class_name.py)
13. Create annotation file (use create_annotation_file.py)