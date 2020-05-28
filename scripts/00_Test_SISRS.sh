#!/bin/bash

# Script takes ~10 minutes with 1 processor

# Set up folder structure using Taxon ID file
python 01_sisrs_folder_setup.py --tid Test_Data/Taxon_IDs

# Move reads to raw read directory
cp -rf Test_Data/* ../Reads/RawReads/

# Trim reads with FastQC analysis with 1 processor
python 02_sisrs_read_trimmer.py --processors 1

# Subset reads for composite genome assembly (assuming a 50Kb genome)
python 03_sisrs_read_subsetter.py --genomesize 50000

# Assemble composite genome with 1 node and 1 processor
python 04_sisrs_composite_genome_assembly.py --nodes 1 --processors 1

# Map and call bases for each species with 1 processor, requiring 2 reads of coverage (minread 2) and all reads supporting a single base (threshold 1)
python 05_sisrs_run_prep.py --processors 1 --minread 2 --threshold 1
python 06_sisrs_run_mapping.py

# Output final alignments, and filter datasets down to those with 0, 1, or 2 species missing per site
python 07_sisrs_output_alignment.py -m 0 1 2
