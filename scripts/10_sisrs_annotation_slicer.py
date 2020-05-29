#!/usr/bin/env python3

# This script finds BED files in BASE_DIR/Reference_Genome/Annotation and slices SISRS sites by annotation
# Arguments: -a/--alignment: Alignment to slice
# Arguments: -l/--locs: LocList for alignment to slice
# Arguments: -n/--name: Prefix for output alignments [Default: Basename of alignment file]
# Arguments: -ow/--overwrite: Process all BED files in Reference_Genome/Annotation, overwriting existing results if present [Default: Process only unprocessed BED files]
# Output: Alignments sliced by BED file
import os
from os import path
import sys
from glob import glob
import subprocess
from itertools import islice
import argparse
import subprocess

# Set script location
script_dir = sys.path[0]
base_dir = path.dirname(path.abspath(script_dir))
reference_dir = base_dir + "/Reference_Genome"
annotation_dir = reference_dir + "/Annotation"
loc_dir = annotation_dir + "/Locs"
slice_dir = annotation_dir + "/Sliced_Alignments"

whole_genome_map = str(glob(reference_dir+'/Mapping/Whole_Genome_Mapping/*_Mapped_NonDup.bed')[0])

if not os.path.isdir(annotation_dir):
    sys.exit("No annotation directory found at "+annotation_dir+" ...exiting...")

if not os.path.isdir(loc_dir):
    os.mkdir(loc_dir)
    
# Get missing information
my_parser = argparse.ArgumentParser()
my_parser.add_argument('-a','--alignment',action='store',required=True)
my_parser.add_argument('-l','--locs',action='store',required=True)
my_parser.add_argument('-n','--name',action='store',nargs="?")
my_parser.add_argument('-ow','--overwrite',action='store_true')
args = my_parser.parse_args()

alignment=os.path.abspath(args.alignment)
locs=os.path.abspath(args.locs)
overwrite=args.overwrite

if args.name is not None:
    dataname = args.name
else:
    dataname = os.path.basename(alignment).rsplit(".",1)[0]
    

bed_file_paths = sorted(glob(annotation_dir+'/*.bed'))

if len(bed_file_paths) < 1:
    sys.exit("No BED files found in "+annotation_dir+" ...exiting...")
    
bed_file_names = [os.path.basename(x).replace('.bed','') for x in bed_file_paths]
existing_locs = sorted(glob(loc_dir+'/*_Locs'))
existing_locs_names = [os.path.basename(x).replace('_Locs','') for x in existing_locs]

if not overwrite:
    bed_file_names = [x for x in bed_file_names if not x in (existing_locs_names)]
    bed_file_paths = [annotation_dir+"/"+x+".bed" for x in bed_file_names]

#Sort annotation BED files
for annotationFile in bed_file_paths:
    sort_annotation = ['sort',
            '-k',
            '1,1',
            '-k2,2n',
            annotationFile,
            '-o',
            annotationFile]
    os.system(' '.join(sort_annotation))

for annotationFile in bed_file_paths:
    annotation=os.path.basename(annotationFile).replace('.bed','')
    output_anno ='{LOCDIR}/{ANNOTATION}_Locs'.format(LOCDIR=loc_dir,ANNOTATION=annotation)
    bed_command=[
        'bedtools',
        'intersect',
        '-sorted',
        '-a',
        '{}'.format(annotationFile),
        '-b',
        '{}'.format(whole_genome_map),
        '-wb',
        '|',
        'awk',
        "'{print",
        '$NF',
        "}'",
        '|',
        'sort',
        '|'
        'uniq'
        '>',
        '{}'.format(output_anno)]
    os.system(' '.join(bed_command))

locsToProcess = [loc_dir+"/"+x+"_Locs" for x in bed_file_names]

for locFile in locsToProcess:
    locname = os.path.basename(locFile).replace('_Locs','')
    if os.stat(locFile).st_size == 0:
        print("No overlapping sites between " + alignment + " and "+ locname)
    else:
        slice_command = ['python',
                '{}/09_sisrs_alignment_slicer.py'.format(script_dir),
                        '--alignment',
                        alignment,
                        '--locs',
                        locs,
                        '--subsetlocs',
                        locFile,
                        '--name',
                        '{DATANAME}_{LOCNAME}'.format(DATANAME=dataname,LOCNAME=locname),
                        '--outputdir',
                        slice_dir]
        output = subprocess.check_output(' '.join(slice_command), shell=True).decode(sys.stdout.encoding).strip()
        with open(slice_dir+"/out_"+dataname+"_"+locname, 'w') as outfile:
            print("\n"+output+"\n")
            print(output,file=outfile)
