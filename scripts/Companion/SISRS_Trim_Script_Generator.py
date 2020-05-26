#!/usr/bin/env python3

# This is a wrapper script for generating trimming scripts for use in SISRS
# All reads for all taxa should be in .fastq.gz format (To change this, find/replace this script, replacing '.fastq.gz' with your chosen extension)
# Paired-end read files must be identically basenamed and end in _1/_2
# Output: (1) Trim script bash file (Trim_Scripts.sh) in 'scripts' directory for all read sets in Reads/RawReads

import os
from os import path
import sys
from glob import glob
import subprocess
from subprocess import check_call

#Set cwd to script location
script_dir = sys.path[0]

#Find BBDuk + Adapter File
cmd = ['which', 'bbduk.sh']
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
o, e = proc.communicate()
bbduk_adapter = path.dirname(o.decode('ascii'))+"/resources/adapters.fa"

#Set RawRead and TrimRead directories based off of script folder location
raw_read_dir = path.dirname(path.abspath(script_dir))+"/Reads/RawReads"
trim_read_dir = path.dirname(path.abspath(script_dir))+"/Reads/TrimReads"

#Find taxa folders within RawRead folder
raw_read_tax_dirs = sorted(glob(raw_read_dir+"/*/"))

#Create folder for BBDuk StdOut
trim_output = raw_read_dir+"/trimOutput/"

#Ensure Trim Log/FastQC output directories not included in taxa list
raw_read_tax_dirs = [x for x in raw_read_tax_dirs if not x.endswith('trimOutput/')]
raw_read_tax_dirs = [x for x in raw_read_tax_dirs if not x.endswith('fastqcOutput/')]

f= open(script_dir+"/Trim_Scripts.sh","w+")

#For each taxa directory...
for tax_dir in raw_read_tax_dirs:

    #List all files and set output dir
    files = sorted(glob(tax_dir+"*.fastq.gz"))
    taxon_ID = path.basename(tax_dir[:-1])
    out_trim_dir = trim_read_dir + "/" + taxon_ID

    left_pairs = list()
    right_pairs = list()
    single_end = list()

    #Find files ending in _1/_2.fastq.gz
    left_files = [s for s in files if "_1.fastq.gz" in s]
    right_files = [s for s in files if "_2.fastq.gz" in s]

    #Strip _1.fastq.gz/_2.fastq.gz and identify pairs based on file name
    left_files = [x.replace('_1.fastq.gz', '') for x in left_files]
    right_files = [x.replace('_2.fastq.gz', '') for x in right_files]
    paired_files = list(set(left_files).intersection(right_files))

    #Reset file names and filter out single-end files
    for pair in paired_files:
        left_pairs.append(pair+"_1.fastq.gz")
        right_pairs.append(pair+"_2.fastq.gz")
    paired_files = sorted(left_pairs + right_pairs)

    single_end = [x for x in files if x not in paired_files]

    #Remove .fastq.gz from lists to make naming easier
    left_pairs = [x.replace('_1.fastq.gz', '') for x in left_pairs]
    right_pairs = [x.replace('_2.fastq.gz', '') for x in right_pairs]
    single_end = [x.replace('.fastq.gz', '') for x in single_end]

    #Trim single-end files if present...
    if len(single_end) > 0:
        for x in single_end:
            se_trim_command = [
                'bbduk.sh',
                'maxns=0',
                'ref={}'.format(bbduk_adapter),
                'qtrim=w',
                'trimq=15',
                'minlength=35',
                'maq=25',
                'in={}'.format(x+'.fastq.gz'),
                'out={}'.format(out_trim_dir+"/"+path.basename(x)+'_Trim.fastq.gz'),
                'k=23',
                'mink=11',
                'hdist=1',
                'hdist2=0',
                'ktrim=r',
                'ow=t',
                '&>',
                '{outDir}out_{fileName}_Trim'.format(outDir=trim_output,fileName=path.basename(x))]
            f.write(' '.join(se_trim_command))
            f.write("\n")

    #Trim paired-end files if present...
    if(len(left_pairs) == len(right_pairs) & len(left_pairs) > 0):
        for x in range(len(left_pairs)):
            file_name = path.basename(left_pairs[x])
            pe_trim_command = [
                'bbduk.sh',
                'maxns=0',
                'ref={}'.format(bbduk_adapter),
                'qtrim=w',
                'trimq=15',
                'minlength=35',
                'maq=25',
                'in={}'.format(left_pairs[x]+'_1.fastq.gz'),
                'in2={}'.format(right_pairs[x]+'_2.fastq.gz'),
                'out={}'.format(out_trim_dir+"/"+path.basename(left_pairs[x])+'_Trim_1.fastq.gz'),
                'out2={}'.format(out_trim_dir+"/"+path.basename(right_pairs[x])+'_Trim_2.fastq.gz'),
                'k=23',
                'mink=11',
                'hdist=1',
                'hdist2=0',
                'ktrim=r',
                'ow=t',
                '&>',
                '{outDir}out_{fileName}_Trim'.format(outDir=trim_output,fileName=file_name)]
            f.write(' '.join(pe_trim_command))
            f.write("\n")

f.close()
