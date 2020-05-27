#!/usr/bin/env python3

# This script preps the folder architecture for a SISRS run.
# Arguments: (1) tid: Path to text file with Taxon IDs on new lines [Default: One level above 'scripts' directory]
# Example: python scripts/folder_setup.py --tid /path/to/TaxonIDs
# Output: Script will create lots of folders, including taxon folders in the RawReads, TrimReads, and SISRS_Run folders

import sys
import os
from os import path
import argparse

# Set script dir and SISRS dir location 
script_dir = sys.path[0]
base_dir = path.dirname(path.abspath(script_dir))

# Get taxon ID file location
my_parser = argparse.ArgumentParser()
my_parser.add_argument('tid',action='store',default=path.dirname(path.abspath(script_dir))+"/Taxon_IDs")
args = my_parser.parse_args()
taxon_ID_file = args.tid

# Read taxon IDs from file
taxa_list = [line.rstrip('\n') for line in open(taxon_ID_file)]

# Make directories
os.mkdir(sisrs_dir+"/Reads")
os.mkdir(sisrs_dir+"/Reads/RawReads")
os.mkdir(sisrs_dir+"/Reads/TrimReads")
os.mkdir(sisrs_dir+"/Reads/SubsetReads")
os.mkdir(sisrs_dir+"/SISRS_Run")
os.mkdir(sisrs_dir+"/Reference_Genome")

for x in taxa_list:
    os.mkdir(sisrs_dir+"/Reads/RawReads/"+x)
    os.mkdir(sisrs_dir+"/Reads/TrimReads/"+x)
    os.mkdir(sisrs_dir+"/SISRS_Run/"+x)
