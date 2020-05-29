#!/usr/bin/env python3

"""Take an entire nexus file, extracts specific sites based on list of names, outputs a phylip-relaxed alignment with DNA data and another with gap data (if present)

    arguments:
    (1) Complete loc list from original alignment to be sliced
    (2) List of desired locs to pull out
    (3) Original alignment to slice
    (4) Dataset name for output files

    output:
    (1) Phylip-relaxed formatted file ending with <dataset>_NoGaps.phylip-relaxed containing only sites in supplied list (and a binary reformatted <dataset>_RecodedGaps.phylip-relaxed if gaps are present)
    (2) Phylip-relaxed formatted file with DNA and non-recoded gap positions suitable for site-split analysis
    (3) File with list of loci names (One for DNA, one for gaps, one with both)
    (4) Nexus-style partition file is also generated if gaps are present for use in tree-building
    """

from __future__ import division
import sys
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_alphabet,generic_dna
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio import AlignIO, SeqIO
import pandas as pd
import argparse

def translateGaps(seq):
    seqList = list(seq)
    seqList = map(lambda x: (0 if x in ['A','a','T','t','C','c','G','g'] else 1 if x=="-" else "N"),seqList)
    seq = ''.join(str(x) for x in seqList)
    return seq

# Get arguments
my_parser = argparse.ArgumentParser()
my_parser.add_argument('-a','--alignment',action='store',required=True)
my_parser.add_argument('-l','--locs',action='store',required=True)
my_parser.add_argument('-s','--subsetlocs',action='store',required=True)
my_parser.add_argument('-n','--name',action='store',required=True)
my_parser.add_argument('-o','--outputdir',action='store',nargs="?")

args = my_parser.parse_args()

alignment=os.path.abspath(args.alignment)
locs=os.path.abspath(args.locs)
subsetLocs=os.path.abspath(args.subsetlocs)
siteID=str(args.name)

dataPath=os.path.dirname(os.path.abspath(subsetLocs))

if args.outputdir is not None:
    outputdir=args.outputdir
else:
    outputdir=os.path.dirname(os.path.abspath(subsetLocs))

if(not os.path.isdir(outputdir)):
    os.mkdir(outputdir)

#Read in loc list of original alignment
originalLocs=pd.read_csv(locs,header=None)

originalLocs.columns=['Original']

#Read in desired sites to be sliced out
with open(subsetLocs, 'r') as f:
    sitesOfInterest = f.read().splitlines()

#Filter original loci to only include sites of interests, create index-based bi-directional dictionary
newLocs = originalLocs[originalLocs['Original'].isin(sitesOfInterest)]
newLocList = newLocs["Original"].tolist()
newLocIndex = newLocs.index.tolist()
biKey = newLocIndex + newLocList
biVal = newLocList + newLocIndex
newLocDict = dict(zip(biKey,biVal))

#Read in original alignment to be sliced
formats = {'nex':'nexus', 'phy':'phylip-relaxed', 'fa':'fasta', 'phylip-relaxed':'phylip-relaxed'}
fformat = formats[alignment.split('.')[-1]]
originalData = AlignIO.read(alignment,fformat)
originalData_Dict = SeqIO.to_dict(originalData)

#Get species
species = originalData_Dict.keys()

#Fill in alignment data
for k in originalData_Dict:
    originalData_Dict[k] = list(originalData_Dict[k].seq)

#Generate list of locs with and without gaps, as well as their index position relative to the original alignment
dnaLocs = []
gapLocs = []

for loc in newLocList:
    if '-' in [originalData_Dict[sp][newLocDict[loc]] for sp in species]:
        gapLocs.append(loc)
    else:
        dnaLocs.append(loc)

#Create DNA dataset
dnaIndex = map(lambda x: newLocDict[x],dnaLocs)

newDNA = originalData[:,0:0]
dnaColumns = []
for column in dnaIndex:
    dnaColumns.append(originalData[:,column])

dnaList = [''.join(s) for s in zip(*dnaColumns)]

i=0
for record in newDNA:
    record.seq=Seq(dnaList[i],generic_dna)
    i+=1

dnaLocfile = open(outputdir + '/' + siteID + '_NoGaps_LocList.txt', 'w')
dnaLocfile.write("\n".join(dnaLocs))
dnaLocfile.close()

SeqIO.write(newDNA, outputdir + '/' + siteID + '_NoGaps.phylip-relaxed', "phylip-relaxed")
print(' - ' + str(originalLocs.shape[0]) + ' sites from \'' + os.path.basename(alignment) + '\' were reduced to ' + str(newDNA.get_alignment_length()) + ' gappless sites (Dataset: ' + siteID + '). A gap-free alignment file has been generated in ' + outputdir + '/' + siteID + '_NoGaps.phylip-relaxed')

#Create gap dataset
if len(gapLocs)>0:
    gapIndex = map(lambda x: newLocDict[x],gapLocs)

    newGap = originalData[:,0:0]
    newGapRecode = originalData[:,0:0]

    gapColumns = []
    gapColumns_Recoded = []

    for column in gapIndex:
        gapColumns.append(originalData[:,column])

    for i in range(len(gapColumns)):
        gapColumns_Recoded.append(translateGaps(gapColumns[i]))

    newGap = originalData[:,0:0]
    newGap_Recoded = originalData[:,0:0]

    gapList = [''.join(s) for s in zip(*gapColumns)]
    gapRecodeList = [''.join(s) for s in zip(*gapColumns_Recoded)]

    i=0
    for record in newGap:
        record.seq=Seq(gapList[i],generic_dna)
        i+=1

    i=0
    for record in newGapRecode:
        record.seq=Seq(gapRecodeList[i],generic_alphabet)
        i+=1

    gapLocfile = open(outputdir + '/' + siteID + '_Gapped_LocList.txt', 'w')
    gapLocfile.write("\n".join(gapLocs))
    gapLocfile.close()

    comboLocfile = open(outputdir + '/' + siteID + '_AllSites_LocList.txt', 'w')
    allLocs = dnaLocs + gapLocs
    comboLocfile.write("\n".join(allLocs))
    comboLocfile.close()

    SeqIO.write(newGapRecode, outputdir + '/' + siteID + '_RecodedGaps.phylip-relaxed', "phylip-relaxed")
    SeqIO.write(newDNA + newGap, outputdir + '/' + siteID + '_AllSites.phylip-relaxed', "phylip-relaxed")

    partitionFile = open(outputdir + '/' + siteID + '_Partition.nex', 'w')
    partitionFile.write('#nexus\n')
    partitionFile.write('begin sets;\n')
    partitionFile.write('\tcharset part1 = ' + outputdir +  '/' + siteID + '_NoGaps.phylip-relaxed:*;\n')
    partitionFile.write('\tcharset part2 = ' + outputdir +  '/' + siteID + '_RecodedGaps.phylip-relaxed:*;\n')
    partitionFile.write('\tcharpartition mine = GTR+G+ASC:part1, GTR2+ASC:part2;\n')
    partitionFile.write('end;')
    partitionFile.close()

    print(' - ' + str(newGapRecode.get_alignment_length()) + ' selected sites from \'' + alignment + '\' contained potentially informative gap positions. A binary coded gap alignment file and a nexus partition file have been generated in ' + outputdir + '/' + siteID + '_RecodedGaps.phylip-relaxed/' + siteID + '_Partition.nex')
    print(' - To facilitate site-split analysis, a DNA/Gap combined alignment has been generated in ' + dataPath + '/' + siteID + '_AllSites.phylip-relaxed')
else:
    print(' - No informative gap positions were identified in the ' + siteID + ' dataset. The gapless alignment can be used for both tree-building and site-split analyses.')

