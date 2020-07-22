#!/usr/bin/env python3

# This script outputs the final alignments
# Arguments: -m/--missing (OPTIONAL): Number(s) of taxa allowed to have missing data in final alignment; Must be less than total species count - 2
# Output: Alignments containing all variable and parsimony-informative sites, along with filtered alignments if requested

import os
from os import path
import sys
from glob import glob
import subprocess
from itertools import islice
import argparse

# Set script location
script_dir = os.path.dirname(os.path.abspath(__file__))

#Set TrimRead + SISRS directories based off of script folder location
sisrs_dir = path.dirname(path.abspath(script_dir))+"/SISRS_Run"
composite_dir = sisrs_dir + '/Composite_Genome'
sisrs_tax_dirs = sorted(glob(sisrs_dir+"/*/"))
sisrs_tax_dirs = [x for x in sisrs_tax_dirs if not x.endswith('Composite_Genome/')]
sisrs_tax_count =  len(sisrs_tax_dirs)
max_missing = sisrs_tax_count - 2

# Get missing information
my_parser = argparse.ArgumentParser()
my_parser.add_argument('-m','--missing',action='store',nargs="*",default=False)
args = my_parser.parse_args()

missing=args.missing

filter_template="""
python SCRIPT_DIR/filter_nexus_for_missing.py SISRS_DIR/alignment.nex MISSING
python SCRIPT_DIR/filter_nexus_for_missing.py SISRS_DIR/alignment_bi.nex MISSING
python SCRIPT_DIR/filter_nexus_for_missing.py SISRS_DIR/alignment_pi.nex MISSING
"""

filter_list=[]

# If no missing values are supplied, don't filter alignments
if not missing:
    filter_list=[]

else:
    for miss in missing:
        if int(miss) < int(max_missing):
            new_filter = filter_template
            keyList = ['SISRS_DIR','SCRIPT_DIR','MISSING']
            keyDict = {'SISRS_DIR':sisrs_dir,'SCRIPT_DIR':script_dir,'MISSING':str(int(miss))}

            for key in keyList:
                new_filter = new_filter.replace(key,keyDict[key])
            filter_list.append(new_filter)
        else:
            print("WARNING: Supplied missing taxa value of "+str(miss)+" would result in 2 or fewer taxa left. Cannot filter...")

sisrs_output_template = """#!/bin/sh
python SCRIPT_DIR/get_alignment.py TWOTAXA SISRS_DIR COMPOSITE_DIR
"""

keyList = ['TWOTAXA','SISRS_DIR','SCRIPT_DIR','COMPOSITE_DIR']
keyDict = {'TWOTAXA':str(len(sisrs_tax_dirs) - 2),'SISRS_DIR':sisrs_dir,'SCRIPT_DIR':script_dir,'COMPOSITE_DIR':composite_dir}

for key in keyList:
    sisrs_output_template = sisrs_output_template.replace(key,keyDict[key])

with open(sisrs_dir+"/Output_Alignment.sh", "w+") as text_file:
    print(sisrs_output_template, file=text_file)
    if len(filter_list)>0:
        for missfilter in filter_list:
            print("\n",file=text_file)
            print(missfilter, file=text_file)

with open(sisrs_dir+"/out_SISRS_Alignment","w") as file:
    cmd = sisrs_dir+'/Output_Alignment.sh'
    subprocess.call(['sh',cmd],stdout=file, stderr=subprocess.PIPE)

with open(sisrs_dir+"/out_SISRS_Log","w") as file:
    file.write("\nRead Mapping and SISRS Site Selection:\n")
    for tax_dir in sisrs_tax_dirs:
        taxa = path.basename(tax_dir[:-1])
        bowtie1 = ['grep','-A4',"'of these'",'{}'.format(tax_dir + "err_" + taxa + "_SISRS"),'|','sed','-n','2,6p']
        bowtie2 = ['grep','-A4',"'of these'",'{}'.format(tax_dir + "err_" + taxa + "_SISRS"),'|','sed','-n','9,13p']
        file.write("\n"+ taxa + " Composite Genome Mapping:\n\n")
        file.write((subprocess.check_output(' '.join(bowtie1),shell=True).decode("UTF8")))
        file.write("\n"+taxa+" Specific Genome Mapping:\n\n")
        file.write((subprocess.check_output(' '.join(bowtie2),shell=True).decode("UTF8")))

        with open(tax_dir + "out_" + taxa + "_SISRS") as f:
            file.write("\n"+taxa+" SISRS Site Selection:\n\n")
            for line in f:
                if(str.startswith(line,'Of ')):
                    file.write(line)
    with open(sisrs_dir + "/out_SISRS_Alignment") as f2:
        file.write("\nSISRS Alignment Filtering:\n\n")
        for line in f2:
            file.write(line)
