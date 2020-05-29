#!/usr/bin/env python3

# This script maps SISRS data (or any FASTA file) onto a reference genome and creates a genome-site map
# Arguments: -p/--processors: Number of available processors per node [Default: 1]
# Arguments: -r/--reference: Bowtie2 index name
# Arguments: -t/--taxon: SISRS taxon ID for sample to be mapped
# Arguments: --bash (OPTIONAL): Write command to BASH script (Reference_Genome_Mapping.sh) in 'scripts' folder
# Arguments: --slurm (OPTIONAL): Add SLURM header to mapping scripts
# Output: 

import os
from os import path
import sys
from glob import glob
import pandas as pd
from Bio import SeqIO
import argparse
import subprocess
from subprocess import check_call

# Set script directory
script_dir = sys.path[0]
base_dir = path.dirname(path.abspath(script_dir))

#Create genome mapping folder
reference_dir = base_dir + "/Reference_Genome"

if(not path.isdir(reference_dir+"/Mapping")):
    os.mkdir(reference_dir+"/Mapping")
mapping_dir =reference_dir+"/Mapping"

# Get arguments
my_parser = argparse.ArgumentParser()
my_parser.add_argument('-p','--processors',action='store',default=1)
my_parser.add_argument('-r','--reference',action='store',required=True)
my_parser.add_argument('-t','--taxon',action='store',required=True)
my_parser.add_argument('--bash',action='store_true')
my_parser.add_argument('--slurm',action='store_true')

args = my_parser.parse_args()

processors = int(args.processors)
reference = str(args.reference)
taxon = str(args.taxon)
taxon_contigs = base_dir + "/SISRS_Run/" + taxon + "/contigs.fa"

check_bash = args.bash
check_slurm = args.slurm

# Create script file

f= open(mapping_dir+"/Reference_Genome_Mapping.sh","w+")

if check_slurm: # Generate SLURM header

    slurm_header = """#!/bin/bash
#SBATCH --job-name="SISRS_Ref"
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=PROCESSORS
cd $SLURM_SUBMIT_DIR
"""
    keyList = ['PROCESSORS']
    keyDict = {'PROCESSORS':str(int(processors))}

    for key in keyList:
        slurm_header = slurm_header.replace(key,keyDict[key])
    f.write(slurm_header)
    f.write("\n")

elif check_bash:
    f.write("#!/bin/bash")
    f.write("\n")

ref_map_template="""
bowtie2 -p PROCESSORS -x REFERENCE_DIR/INDEX_NAME -f -U CONTIGS | samtools view -Su -@ PROCESSORS -F 4 - | samtools sort -@ PROCESSORS - -o MAPPING_DIR/TAXON_Temp.bam

samtools view -@ PROCESSORS -H MAPPING_DIR/TAXON_Temp.bam > MAPPING_DIR/TAXON_Header.sam

samtools view -@ PROCESSORS MAPPING_DIR/TAXON_Temp.bam | grep -v "XS:" | cat MAPPING_DIR/TAXON_Header.sam - | samtools view -@ PROCESSORS -b - > MAPPING_DIR/TAXON.bam

rm MAPPING_DIR/*Temp.bam
rm MAPPING_DIR/*Header.sam

samtools view MAPPING_DIR/TAXON.bam | awk 'BEGIN {OFS = "\t"} { print $1, $3, $4, $2, $6}' > MAPPING_DIR/TAXON_MapData.tsv

cut -f1 MAPPING_DIR/TAXON_MapData.tsv | sort > MAPPING_DIR/Uniquely_Mapping_Contigs

python SCRIPT_DIR/Genome_Mapper.py MAPPING_DIR/TAXON_MapData.tsv

sort -k 1,1 -k2,2n MAPPING_DIR/Whole_Genome_Mapping/WholeGenome_TAXON_Mapped_NonDup.bed -o MAPPING_DIR/Whole_Genome_Mapping/WholeGenome_TAXON_Mapped_NonDup.bed

cut -f4 MAPPING_DIR/Whole_Genome_Mapping/WholeGenome_TAXON_Mapped_NonDup.bed | sort > MAPPING_DIR/TAXON_NonDup_Genome_Sites

bedtools merge -i MAPPING_DIR/Whole_Genome_Mapping/WholeGenome_TAXON_Mapped_NonDup.bed > MAPPING_DIR/Whole_Genome_Mapping/WholeGenome_TAXON_Mapped_NonDupMerged.bed
"""

keyList = ['PROCESSORS','REFERENCE_DIR','INDEX_NAME','CONTIGS','MAPPING_DIR','TAXON','SCRIPT_DIR']
keyDict = {'PROCESSORS':str(processors),
           'REFERENCE_DIR':reference_dir,
           'INDEX_NAME':str(reference),
           'CONTIGS':taxon_contigs,
           'MAPPING_DIR':mapping_dir,
           'TAXON':taxon,
           'SCRIPT_DIR':script_dir}

for key in keyList:
    ref_map_template = ref_map_template.replace(key,keyDict[key])

f.write(ref_map_template)
if check_bash or check_slurm:
    f.close()
else:
    with open(mapping_dir+"/out_Reference_Genome_Mapping","w") as file:
        f.close()
        cmd = mapping_dir+'/Reference_Genome_Mapping.sh'
        subprocess.call(['sh',cmd],stdout=file, stderr=subprocess.PIPE)
