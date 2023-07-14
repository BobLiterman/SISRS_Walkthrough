#!/usr/bin/env python3

# This script prepares data for a SISRS run by setting the data up and creating mapping scripts
# Contigs are renamed and moved to the SISRS_Run/Composite_Genome directory
# The composite genome is indexed by Bowtie2 and Samtools
# SISRS scripts are generated from a template and saved to the SISRS_Run/TAXA folder

import os
from os import path
import sys
from glob import glob
import pandas as pd

#Set cwd to script location
script_dir = sys.path[0]
sisrs_dir = path.dirname(path.abspath(script_dir))+"/SISRS_Run"
composite_dir =sisrs_dir+"/Composite_Genome"
sisrs_tax_dirs = sorted(glob(sisrs_dir+"/*/"))

sisrs_template = """#!/bin/bash
#PBS -N TAXA
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=PROCESSORS
#PBS -o out_TAXA_SISRS
#PBS -e err_TAXA_SISRS

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
cd $PBS_O_WORKDIR
cp -as COMPOSITE_DIR/ref .

bbwrap.sh in=READS maxindel=99 strictmaxindel=t sam=1.3 ambiguous=toss out=TAXA.sam append=t
samtools view -Su -@ PROCESSORS -F 4 TAXA.sam | samtools sort -@ PROCESSORS - -o SISRS_DIR/TAXA/TAXA.bam
samtools mpileup -f COMPOSITE_GENOME SISRS_DIR/TAXA/TAXA.bam > SISRS_DIR/TAXA/TAXA.pileups

python SCRIPT_DIR/specific_genome.py SISRS_DIR/TAXA COMPOSITE_GENOME

samtools faidx SISRS_DIR/TAXA/contigs.fa
rm -rf ref TAXA.sam TAXA.bam
bbmap.sh ref=contigs.fa

bbwrap.sh in=READS maxindel=99 strictmaxindel=t sam=1.3 ambiguous=toss out=TAXA.sam append=t
samtools view -Su -@ PROCESSORS -F 4 TAXA.sam | samtools sort -@ PROCESSORS - -o SISRS_DIR/TAXA/TAXA.bam
samtools index SISRS_DIR/TAXA/TAXA.bam

samtools mpileup -f COMPOSITE_GENOME SISRS_DIR/TAXA/TAXA.bam > SISRS_DIR/TAXA/TAXA.pileups

python SCRIPT_DIR/get_pruned_dict.py SISRS_DIR/TAXA COMPOSITE_DIR MINREAD THRESHOLD

"""

#Create links to trimmed read files in SISRS_Run directory
for tax_dir in sisrs_tax_dirs:
    taxa = path.basename(tax_dir[:-1])

    new_sisrs = sisrs_template

    keyList = ['PROCESSORS','BOWTIE2-INDEX','COMPOSITE_GENOME','SCRIPT_DIR','MINREAD','THRESHOLD','TAXA','SISRS_DIR','COMPOSITE_DIR','READS']
    keyDict = {'PROCESSORS':str(int(sys.argv[1])),
               'BOWTIE2-INDEX':composite_dir+"/contigs",
               'COMPOSITE_GENOME':composite_dir+"/contigs.fa",
               'SCRIPT_DIR':script_dir,
               'MINREAD':str(int(sys.argv[2])),
               'THRESHOLD':str(float(sys.argv[3])),
               'TAXA':taxa,
               'SISRS_DIR':sisrs_dir,
               'COMPOSITE_DIR':composite_dir,
               'READS':",".join(glob(tax_dir+"*.fastq.gz"))}
    for key in keyList:
        new_sisrs = new_sisrs.replace(key,keyDict[key])
    with open(sisrs_dir+"/"+taxa+"/"+taxa+".sh", "w") as text_file:
        print(new_sisrs, file=text_file)
    os.system('chmod +x '+sisrs_dir+"/"+taxa+"/"+taxa+".sh")

