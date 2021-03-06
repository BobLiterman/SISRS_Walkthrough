#!/bin/bash

# Script takes ~10 minutes with 1 processor

# Set up folder structure using Taxon ID file
python 01_sisrs_folder_setup.py --tid Test_Data/Taxon_IDs

# Move reads to raw read directory
cp -rf Test_Data/RawReads/* ../Reads/RawReads/

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
python 07_sisrs_output_alignment.py --missing 0 1 2

# Map SISRS contigs from Homo sapiens (HomSap) onto Homo sapiens Chromosome 17 (HomSap_Chr17) using 1 processors
cp -rf Test_Data/Reference_Genome/* ../Reference_Genome
python 08_sisrs_reference_mapping.py --processors 1 --reference HomSap_Chr17 --taxon HomSap

# Extract SISRS sites that map uniquely to the reference genome
python 09_sisrs_alignment_slicer.py --alignment ../SISRS_Run/alignment_bi_m0_nogap.phylip-relaxed --locs ../SISRS_Run/alignment_bi_locs_m0_nogap.txt --subsetlocs ../Reference_Genome/Mapping/HomSap_NonDup_Genome_Sites --name Hg17_Mapped --outputdir ../Reference_Genome/Hg17_Mapped

# Slice SISRS alignment (alignment_bi_m2.phylip-relaxed) down to just sites from CDS or gene annotations
python 10_sisrs_annotation_slicer.py --alignment ../SISRS_Run/alignment_bi_m2_nogap.phylip-relaxed --locs ../SISRS_Run/alignment_bi_locs_m2_nogap.txt --name Biallelic_2Missing
