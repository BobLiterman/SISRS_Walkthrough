#!/usr/bin/env python3

# This script calls a Ray genome assembly on all reads in the SubsetReads directory
# Ray requires working MPI via mpirun, even if running on  one node
# This script also calls Bowtie2 and Samtools for genome indexing
# Ensure Ray and Bowtie2 are compiled  with gzip support
# Input: (1) Script run-mode (2) Number of nodes for Ray. (3) Number of processors per node. 
# Script run-mode: 'r': Run Ray command directly as a result of running this script
# Script run-mode: 'sh': Write command to 'Assemble_Ray.sh' in scripts directory
# Script run-mode: 'slurm': Write command to 'Assemble_Ray.sh' in scripts directory, with a SLURM header appended
# Output: Ray assembly will be built in <basedir>/Ray_Composite_Genome (Contigs.fasta)

import os
from os import path
import sys
from glob import glob
import subprocess
from subprocess import check_call

if len(sys.argv) is not 4:
    sys.exit('Ray script requires 3 arguments: run-mode, node count, and processors per node')

#Set cwd to script location
script_dir = sys.path[0]

node_count = int(sys.argv[2])
processors_per_node = int(sys.argv[3])

#Get run mode
if sys.argv[1] is 'r':
    run_mode = 'run'

elif sys.argv[1] is 'sh':
    run_mode = 'script'

elif sys.argv[1] is 'slurm':
    run_mode = 'slurm'
    
    slurm_header="""#!/bin/bash
    #SBATCH --job-name="Ray_Assembly"
    #SBATCH --time=96:00:00  # walltime limit (HH:MM:SS)
    #SBATCH --nodes=NODES   # number of nodes
    #SBATCH --ntasks-per-node=PROCESSORS   # processor core(s) per node 
    #SBATCH --exclusive
    # LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
    cd $SLURM_SUBMIT_DIR
    """

    keyList = ['NODES','PROCESSORS']
    keyDict = {'NODES':str(node_count),'PROCESSORS'=str(processors_per_node)}
    for key in keyList:
        slurm_header = slurm_header.replace(key,keyDict[key])
else:
    sys.exit("The first argument must be 'r' (run), 'sh' (save script), or 'slurm' (save script with SLURM header).")

mpi_processor = node_count*processors_per_node

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

if run_mode is 'r':
    os.system(" ".join(ray_command))

elif run_mode is 'sh':
    with open(script_dir+"/Assemble_Ray.sh", "w+") as text_file:
        print" ".join(ray_command), file=text_file)

elif run_mode is 'slurm':
    with open(script_dir+"/Assemble_Ray.sh", "w+") as text_file:
        print(slurm_header, file=text_file)
        print(" ".join(ray_command), file=text_file)
