#!/usr/bin/env python3

# This script prepares data for a SISRS run by setting the data up and creating mapping scripts
# Contigs are renamed and moved to the SISRS_Run/Composite_Genome directory
# The composite genome is indexed by Bowtie2 and Samtools
# SISRS scripts are generated from a template and saved to the SISRS_Run/TAXA folder
# Arguments: -p/--processors: Number of available processors per node [Default: 1]
# Arguments: -m/--minread: Minimum read coverage to call a base [Default: 3]
# Arguments: -t/--threshold: Percent homozygosity required to call a base [Default: 1 (all reads must support a single base)]
# Arguments: --slurm (OPTIONAL): Add SLURM header to mapping scripts
# Output: Creates composite genome data and mapping scripts into 'SISRS_Run' directory

import os
from os import path
import sys
from glob import glob
import pandas as pd
from Bio import SeqIO
import argparse

# Set script directory
script_dir = sys.path[0]

# Get arguments
my_parser = argparse.ArgumentParser()
my_parser.add_argument('-p','--processors',action='store',default=1,nargs="?")
my_parser.add_argument('-m','--minread',action='store',default=3,nargs="?")
my_parser.add_argument('-t','--threshold',action='store',default=1,nargs="?")
my_parser.add_argument('--slurm',action='store_true')

args = my_parser.parse_args()

processors = args.processors
minread = int(args.minread)
threshold = float(args.threshold)
check_slurm = args.slurm

#Set TrimRead + SISRS directories based off of script folder location
trim_read_dir = path.dirname(path.abspath(script_dir))+"/Reads/TrimReads"
trim_read_tax_dirs = sorted(glob(trim_read_dir+"/*/"))
ray_dir = path.dirname(path.abspath(script_dir))+"/Ray_Composite_Genome"
sisrs_dir = path.dirname(path.abspath(script_dir))+"/SISRS_Run"

trim_read_tax_dirs = [x for x in trim_read_tax_dirs if not x.endswith('trimOutput/')]
trim_read_tax_dirs = [x for x in trim_read_tax_dirs if not x.endswith('fastqcOutput/')]
trim_read_tax_dirs = sorted([x for x in trim_read_tax_dirs if not x.endswith('subsetOutput/')])

#Create composite genome folder
if(not path.isdir(sisrs_dir+"/Composite_Genome")):
    os.mkdir(sisrs_dir+"/Composite_Genome")
composite_dir =sisrs_dir+"/Composite_Genome"

#Create renamed contig file in SISRS_Run/Composite_Genome
#To do: Generalize beyond Ray
rename_command = [
    'rename.sh',
    'in={}'.format(ray_dir+"/Contigs.fasta"),
    'out={}'.format(composite_dir+'/contigs.fa'),
    'prefix=SISRS',
    'addprefix=t',
    'ow=t']
os.system(" ".join(rename_command))

#Index composite genome using Bowtie2 and Samtools
bowtie_command = [
    'bowtie2-build',
    '{}'.format(composite_dir+'/contigs.fa'),
    '{}'.format(composite_dir+'/contigs'),
    '-p',
    '{}'.format(str(processors))]
os.system(" ".join(bowtie_command))

samtools_command = [
    'samtools',
    'faidx',
    '{}'.format(composite_dir+'/contigs.fa')]
os.system(" ".join(bowtie_command))

# Get sequence lengths from assembly
with open(composite_dir+'/contigs_SeqLength.tsv', "w") as lengthfile:
    for seq_record in SeqIO.parse(composite_dir+'/contigs.fa',"fasta"):
        lengthfile.write(str(seq_record.id)+"\t"+str(len(seq_record))+"\n")

#Create file with an entry for each individual site in the alignment
siteCount=0
locListFile = open(composite_dir+'/contigs_LocList','a+')
with open(composite_dir +"/contigs_SeqLength.tsv","r") as filein:
    for line in iter(filein):
        splitline=line.split()
        for x in range(1,(int(splitline[1])+1)):
            locListFile.write((splitline[0] +'/'+str(x)+'\n'))
            siteCount+=1
locListFile.close()

slurm_header = """#!/bin/bash
#SBATCH --job-name="TAXA"
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=PROCESSORS
cd $SLURM_SUBMIT_DIR
"""

sisrs_template = """bowtie2 -p PROCESSORS -N 1 --local -x BOWTIE2-INDEX -U READS | samtools view -Su -@ PROCESSORS -F 4 - | samtools sort -@ PROCESSORS - -o SISRS_DIR/TAXA/TAXA_Temp.bam
samtools view -@ PROCESSORS -H SISRS_DIR/TAXA/TAXA_Temp.bam > SISRS_DIR/TAXA/TAXA_Header.sam
samtools view -@ PROCESSORS SISRS_DIR/TAXA/TAXA_Temp.bam | grep -v "XS:" | cat SISRS_DIR/TAXA/TAXA_Header.sam - | samtools view -@ PROCESSORS -b - > SISRS_DIR/TAXA/TAXA.bam

rm SISRS_DIR/TAXA/TAXA_Temp.bam
rm SISRS_DIR/TAXA/TAXA_Header.sam

samtools mpileup -f COMPOSITE_GENOME SISRS_DIR/TAXA/TAXA.bam > SISRS_DIR/TAXA/TAXA.pileups

python SCRIPT_DIR/specific_genome.py SISRS_DIR/TAXA COMPOSITE_GENOME

samtools faidx SISRS_DIR/TAXA/contigs.fa

bowtie2-build SISRS_DIR/TAXA/contigs.fa SISRS_DIR/TAXA/contigs -p PROCESSORS
bowtie2 -p PROCESSORS -N 1 --local -x SISRS_DIR/TAXA/contigs -U READS | samtools view -Su -@ PROCESSORS -F 4 - | samtools sort -@ PROCESSORS - -o SISRS_DIR/TAXA/TAXA_Temp.bam
samtools view -@ PROCESSORS -H SISRS_DIR/TAXA/TAXA_Temp.bam > SISRS_DIR/TAXA/TAXA_Header.sam
samtools view -@ PROCESSORS SISRS_DIR/TAXA/TAXA_Temp.bam | grep -v "XS:" | cat SISRS_DIR/TAXA/TAXA_Header.sam - | samtools view -@ PROCESSORS -b - > SISRS_DIR/TAXA/TAXA.bam

rm SISRS_DIR/TAXA/TAXA_Temp.bam
rm SISRS_DIR/TAXA/TAXA_Header.sam

samtools index SISRS_DIR/TAXA/TAXA.bam
samtools mpileup -f COMPOSITE_GENOME SISRS_DIR/TAXA/TAXA.bam > SISRS_DIR/TAXA/TAXA.pileups

python SCRIPT_DIR/get_pruned_dict.py SISRS_DIR/TAXA COMPOSITE_DIR MINREAD THRESHOLD
"""

#Create links to trimmed read files in SISRS_Run directory
for tax_dir in trim_read_tax_dirs:
    taxa = path.basename(tax_dir[:-1])
    read_link_command = [
        'cd',
        '{}'.format(sisrs_dir+"/"+taxa),
        ';',
        'cp',
        '-as',
        '{}/*.fastq.gz'.format(tax_dir),
        '.']
    os.system(" ".join(read_link_command))
    
    new_slurm = slurm_header
    new_sisrs = sisrs_template

    keyList = ['PROCESSORS','BOWTIE2-INDEX','COMPOSITE_GENOME','SCRIPT_DIR','MINREAD','THRESHOLD','TAXA','SISRS_DIR','COMPOSITE_DIR','READS']
    keyDict = {'PROCESSORS':str(processors),
               'BOWTIE2-INDEX':composite_dir+"/contigs",
               'COMPOSITE_GENOME':composite_dir+"/contigs.fa",
               'SCRIPT_DIR':script_dir,
               'MINREAD':str(minread),
               'THRESHOLD':str(threshold),
               'TAXA':taxa,
               'SISRS_DIR':sisrs_dir,
               'COMPOSITE_DIR':composite_dir,
               'READS':",".join(glob(tax_dir+"*.fastq.gz"))}
    for key in keyList:
        new_sisrs = new_sisrs.replace(key,keyDict[key])
        new_slurm = new_slurm.replace(key,keyDict[key])
    
    with open(sisrs_dir+"/"+taxa+"/"+taxa+".sh", "w") as text_file:
        if check_slurm:
            print(new_slurm,file=text_file)
            print(new_sisrs, file=text_file)
        else:
            print("#!/bin/bash", file=text_file)
            print(new_sisrs, file=text_file)
            os.system('chmod +x '+sisrs_dir+"/"+taxa+"/"+taxa+".sh")
