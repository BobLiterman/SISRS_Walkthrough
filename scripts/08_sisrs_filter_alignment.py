#!/usr/bin/env python3

# This script prepares data for a SISRS run by setting the data up and creating mapping scripts
# Contigs are renamed and moved to the SISRS_Run/Composite_Genome directory
# The composite genome is indexed by Bowtie2 and Samtools
# SISRS scripts are generated from a template and saved to the SISRS_Run/TAXA folder

import os
from os import path
import sys
from glob import glob
import subprocess
from itertools import islice

#Set cwd to script location
script_dir = sys.path[0]

#Set TrimRead + SISRS directories based off of script folder location
sisrs_dir = path.dirname(path.abspath(script_dir))+"/SISRS_Run"
composite_dir = sisrs_dir + '/Composite_Genome'
sisrs_tax_dirs = sorted(glob(sisrs_dir+"/*/"))
sisrs_tax_dirs = [x for x in sisrs_tax_dirs if not x.endswith('Composite_Genome/')]
missing_taxa = sys.argv[1]
max_missing = len(sisrs_tax_dirs) - 2

if int(missing_taxa) >= max_missing:
    sys.exit('Specified allowable missing taxa is at or above allowable values based on species counts (Must have 2+ species)')

sisrs_filter_template = """#!/bin/sh
python SCRIPT_DIR/filter_nexus_for_missing_nogap.py SISRS_DIR/alignment_bi.nex MISSING
python SCRIPT_DIR/filter_nexus_for_missing_nogap.py SISRS_DIR/alignment_pi.nex MISSING

grep -oe "SISRS_[^/]*" SISRS_DIR/alignment_bi_locs_mMISSING_nogap.txt | uniq -c | sort -k1 -nr | awk '{print $2}' > SISRS_DIR/alignment_bi_locs_mMISSING_nogap_Clean.txt
grep -oe "SISRS_[^/]*" SISRS_DIR/alignment_pi_locs_mMISSING_nogap.txt | uniq -c | sort -k1 -nr | awk '{print $2}' > SISRS_DIR/alignment_pi_locs_mMISSING_nogap_Clean.txt"""

keyList = ['SISRS_DIR','SCRIPT_DIR','MISSING']
keyDict = {'SISRS_DIR':sisrs_dir,'SCRIPT_DIR':script_dir,'MISSING':str(missing_taxa)}

for key in keyList:
    sisrs_filter_template = sisrs_filter_template.replace(key,keyDict[key])
with open(sisrs_dir+"/Filter_Alignment_m"+str(missing_taxa)+".sh", "w") as text_file:
    print(sisrs_filter_template, file=text_file)

with open(sisrs_dir+"/out_Filter_Alignment_m"+str(missing_taxa),"w") as file:
    cmd = sisrs_dir+"/Filter_Alignment_m"+str(missing_taxa)+".sh"
    subprocess.call(['sh',cmd],stdout=file, stderr=subprocess.PIPE)
