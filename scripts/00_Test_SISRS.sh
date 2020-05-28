#!/bin/bash

# Script takes ~10 minutes with 1 processor

# Set up folder structure using Taxon ID file
python 01_sisrs_folder_setup.py --tid Test_Data/Taxon_IDs

# Move reads to raw read directory
cp -rf Test_Data/RawReads/* ../Reads/RawReads/

# Trim reads with FastQC analysis with 1 processor
python 02_sisrs_read_trimmer.py --processors 7

# Subset reads for composite genome assembly (assuming a 50Kb genome)
python 03_sisrs_read_subsetter.py --genomesize 50000

# Assemble composite genome with 1 node and 1 processor
python 04_sisrs_composite_genome_assembly.py --nodes 1 --processors 7

# Map and call bases for each species with 1 processor, requiring 2 reads of coverage (minread 2) and all reads supporting a single base (threshold 1)
python 05_sisrs_run_prep.py --processors 7 --minread 2 --threshold 1
python 06_sisrs_run_mapping.py

# Output final alignments, and filter datasets down to those with 0, 1, or 2 species missing per site
python 07_sisrs_output_alignment.py --missing 0 1 2

# Map SISRS contigs from Homo sapiens (HomSap) onto Homo sapiens Chromosome 17 (HomSap_Chr17) using 1 processors
cp /home/ralubuntu/Work/test_SISRS/SISRS_Walkthrough/scripts/Test_Data/Reference_Genome/*bt2 ../Reference_Genome
python 08_sisrs_reference_mapping.py --processors 1 --reference HomSap_Chr17 --taxon HomSap
