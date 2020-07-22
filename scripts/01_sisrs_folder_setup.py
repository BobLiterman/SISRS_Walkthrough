#!/usr/bin/env python3

# This script preps the folder architecture for a SISRS run. Fill in Taxon_IDs file in base directory to create subfolders.
# Output: Script will create lots of folders, including taxon folders in the RawReads, TrimReads, and SISRS_Run folders

import sys
import os
from os import path
import argparse

# Set script dir and SISRS dir location 
script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = path.dirname(path.abspath(script_dir))

# Get taxon ID file location
my_parser = argparse.ArgumentParser()
my_parser.add_argument('--tid',action='store',default=path.dirname(path.abspath(script_dir))+"/Taxon_IDs",nargs="?")
args = my_parser.parse_args()
taxon_ID_file = args.tid

# Read taxon IDs from file
taxa_list = [line.rstrip('\n') for line in open(taxon_ID_file)]

# Make directories
os.mkdir(base_dir+"/Reads")
os.mkdir(base_dir+"/Reads/RawReads")
os.mkdir(base_dir+"/Reads/TrimReads")
os.mkdir(base_dir+"/Reads/SubsetReads")
os.mkdir(base_dir+"/SISRS_Run")
os.mkdir(base_dir+"/Reference_Genome")
os.mkdir(base_dir+"/Reference_Genome/Annotation")

for x in taxa_list:
    os.mkdir(base_dir+"/Reads/RawReads/"+x)
    os.mkdir(base_dir+"/Reads/TrimReads/"+x)
    os.mkdir(base_dir+"/SISRS_Run/"+x)
