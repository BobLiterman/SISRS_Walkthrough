#!/usr/bin/env python3

# This script calls a Ray genome assembly on all reads in the SubsetReads directory
# Ray requires working MPI via mpirun, even if running on one node
# This script also calls Bowtie2 and Samtools for genome indexing
# Ensure Ray and Bowtie2 are compiled  with gzip support
#
# Arguments: -n/--nodes: Number of nodes for Ray [Default: 1]
# Arguments: -p/--processors: Number of processors per node [Default: 1]
# Arguments: --bash (Optional): Instead of running Ray, save a BASH script (Ray_Composite.sh) in 'scripts' directory
# Arguments: --slurm (Optional): Instead of running Ray, save a SLURM script (Ray_Composite.sh) in 'scripts' directory [Note: SLURM script has a default time of 48h; Change if necessary]

# Output: Ray assembly will be built in <basedir>/Ray_Composite_Genome (Contigs.fasta)

import os
from os import path
import sys
from glob import glob
import subprocess
from subprocess import check_call
import argparse

#Set cwd to script location
script_dir = os.path.dirname(os.path.abspath(__file__))

# Get arguments
my_parser = argparse.ArgumentParser()
my_parser.add_argument('-p','--processors',action='store',default=1)
my_parser.add_argument('-n','--nodes',action='store',default=1)
my_parser.add_argument('--bash',action='store_true')
my_parser.add_argument('--slurm',action='store_true')

args = my_parser.parse_args()

processors = int(args.processors)
nodes = int(args.nodes)
check_bash = args.bash
check_slurm = args.slurm

# Set total processor count
mpi_processor = processors*nodes

# Create script file
if check_slurm or check_bash:
    f= open(script_dir+"/Ray_Composite.sh","w+")

if check_slurm: # Generate SLURM header

    slurm_header = """#!/bin/bash
#SBATCH --job-name="SISRS_Ray"
#SBATCH --time=48:00:00
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=PROCESSORS
cd $SLURM_SUBMIT_DIR
"""
    keyList = ['NODES','PROCESSORS']
    keyDict = {'NODES':str(int(nodes)),'PROCESSORS':str(int(processors))}

    for key in keyList:
        slurm_header = slurm_header.replace(key,keyDict[key])
    f.write(slurm_header)
    f.write("\n")

elif check_bash:
    f.write("#!/bin/bash")
    f.write("\n")

#Set SubsetRead and Genome directories based off of script folder location
subset_read_dir = path.dirname(path.abspath(script_dir))+"/Reads/SubsetReads"
ray_genome_dir = path.dirname(path.abspath(script_dir))+"/Ray_Composite_Genome"

subset_reads = glob(subset_read_dir+"/*.gz")

ray_command = [
    'mpirun',
    '-n','{}'.format(str(mpi_processor)),
    'Ray',
    '-k',
    '31',
    '{}'.format("-s " + " -s ".join(subset_reads)),
    '-o',
    '{}'.format(ray_genome_dir)]

if check_slurm or check_bash:
    f.write(' '.join(ray_command))
    f.close()
else:
    os.system(" ".join(ray_command))

