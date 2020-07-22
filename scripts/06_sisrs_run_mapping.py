#!/usr/bin/env python3

# This script runs the mapping steps of a SISRS run sequentially
# WARNING: These scripts can be run in parallel (e.g. on an HPC-type system),  but this script runs them sequentially
import os
from os import path
import sys
from glob import glob
import subprocess

#Set cwd to script location
script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = path.dirname(path.abspath(script_dir))

sisrs_dir = base_dir+"/SISRS_Run"
sisrs_tax_dirs = sorted(glob(sisrs_dir+"/*/"))
sisrs_tax_dirs = [x for x in sisrs_tax_dirs if not x.endswith('Composite_Genome/')]

#Create links to trimmed read files in SISRS_Run directory
for tax_dir in sisrs_tax_dirs:
    taxa = path.basename(tax_dir[:-1])
    with open(tax_dir+'out_'+taxa+'_SISRS',"w") as file:
        with open(tax_dir+'err_'+taxa+'_SISRS',"w") as file2:
            cmd = tax_dir+taxa+'.sh'
            subprocess.call(['sh',cmd],stdout=file, stderr=file2)
    print("Completed SISRS filtering for "+taxa + "...\n")
