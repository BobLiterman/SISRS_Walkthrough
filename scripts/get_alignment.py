#!/usr/bin/env python3

"""
    output alignment of sites useful for phylogenetics (nexus format)
    can specify max number of species missing for each site

    arguments:
        num_missing -- the number of species allowed to be missing at a site so that site will show in output alignment
        sisrs_dir -- Path to SISRS run directory (no traililng /)
        composite_dir -- Path to SISRS composite genome directory (no trailing /)

    output:
        alignment.nex : Nexus formatted alignment (including a header with composite genome position); each site has up to num_missing missing data
        alignment_pi_singletons.nex : Of these, variable sites with possible singletons
        alignment_pi.nex : Of these, variable sites without singletons
        alignment_bi.nex : Of these, variable, biallelic sites without singletons
        alignment_singletons.nex : Sites where the only variation is singleton in nature

"""

import os
import re
import glob
from collections import Counter,defaultdict
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

#########################
class Loc:
    def __init__(self,scaff_loc,flag):
        self.scaff_loc = scaff_loc
        self.flag = flag

class Alignment:
    def __init__(self,locations=[],species_data=dict(),flag=[],singleton=[],othersingle=[],anysingle=[],biallelic=[]):
        self.locations = locations
        self.species_data = species_data
        self.flag = flag
        self.singleton = singleton
        self.othersingle = othersingle
        self.anysingle = anysingle
        self.biallelic = biallelic

    def numsnps(self):
        print(str(len(self.locations))+' total variable sites [alignment.nex]')
        for i in range(len(self.locations)):
            
            # Get valid bases
            bases = [self.species_data[sp][i] for sp in self.species_data if self.species_data[sp][i] in ['A','C','G','T','-']]     #bases for that site
            
            # Count bases and sort by occurence
            c = Counter(bases).most_common(5)
            unique_bases = len(c)
            self.flag.append(unique_bases)
            
            sub_list = c[1:unique_bases]
            
            # Check for singletons
            if all(x[1] == 1 for x in c):
                self.singleton.append(1)
                self.othersingle.append(0)
                self.anysingle.append(1)
                self.biallelic.append(0)

            elif all([x[1] == 1 for x in sub_list]):
                self.singleton.append(1)
                self.othersingle.append(0)
                self.anysingle.append(1)
                self.biallelic.append(0)

            elif any([x[1] == 1 for x in sub_list]):
                self.singleton.append(0)
                self.othersingle.append(1)
                self.anysingle.append(1)
                self.biallelic.append(0)
    
            else:
                if unique_bases == 2:
                    self.biallelic.append(1)
                else:
                    self.biallelic.append(0)
                    
                self.singleton.append(0)
                self.othersingle.append(0)
                self.anysingle.append(0)

        print(str(self.singleton.count(1))+' variable sites are singletons (Invariant except for singletons) [alignment_singletons.nex]')
        print(str(self.othersingle.count(1))+' variable sites contain singletons (Variable, but with 1 or more singletons) [alignment_pi_singletons.nex]')
        print(str(self.anysingle.count(0))+' variable sites contain no singletons (Variable, parsimony-informative [alignment_pi.nex]')
        print(str(self.biallelic.count(1))+' variable sites are biallelic, parsimony-informative sites [alignment_bi.nex]')

        return self.flag.count(2)       # number of biallelic sites

def get_phy_sites(sisrs_dir,composite_dir,num_missing):

    #Fetch contig data
    contigList=glob.glob(composite_dir +'/contigs_LocList')
    assert len(contigList) > 0, 'Total site list not found in assembly folder'

    #Fetch sorted species data
    dataLists = sorted(glob.glob(sisrs_dir+'/*/*_LocList'))
    dataLists = [x for x in dataLists if 'contigs_LocList' not in x]
    splist=[os.path.basename(os.path.dirname(path)) for path in dataLists]
    speciesCount=len(dataLists)
    assert len(dataLists) > 0, 'No species had data from the pileup'

    allLists = contigList+dataLists

    alignment = Alignment()
    alignment.species_data = {species: [] for species in splist}

    files = [open(i, "r") for i in allLists]
    for rows in zip(*files):
        rowList = list(map(lambda foo: foo.replace('\n', ''), list(rows)))
        speciesData = rowList[1:(speciesCount+1)]
        if speciesData.count("N")<=num_missing and len(set(filter(lambda a: a != "N", speciesData)))>1:
            alignment.locations.append(rowList[0])
            for j in range(0,(speciesCount)):
                alignment.species_data[splist[j]].append(speciesData[j])
    return alignment

def write_alignment(fi,alignment,numbi):
    spp = sorted(alignment.species_data.keys())
    ntax = str(len(alignment.species_data))

    #Process alignment.nex
    ALIGNMENT=open(fi,'w')
    ALIGNMENT.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(len(alignment.locations))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
    ALIGNMENT.write('[ '+ " ".join(alignment.locations)+' ]'+"\n")
    for species in spp: #write sequences for each species
        ALIGNMENT.write(species+"\t"+"".join(alignment.species_data[species])+"\n")
    ALIGNMENT.write(';\nend;')
    ALIGNMENT.close()

    #Process alignment_bi.nex
    ALIGNMENTBI=open(fi.replace('.nex','_bi.nex'),'w')
    bi_loc = [alignment.locations[i] for i in range(len(alignment.locations)) if alignment.biallelic[i] == 1]
    bi_sp_data={}
    for species in spp:
        bi_sp_data[species] = [alignment.species_data[species][i] for i in range(len(alignment.locations)) if alignment.biallelic[i] == 1]

    ALIGNMENTBI.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(len(bi_loc))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
    ALIGNMENTBI.write('[ '+ " ".join(bi_loc)+' ]'+"\n")
    for species in spp: #write sequences for each species
        ALIGNMENTBI.write(species+"\t"+("".join(bi_sp_data[species]))+"\n")
    ALIGNMENTBI.write(';\nend;')
    ALIGNMENTBI.close()

    #Process alignment_pi.nex
    ALIGNMENTPI=open(fi.replace('.nex','_pi.nex'),'w')
    pi_loc = [alignment.locations[i] for i in range(len(alignment.locations)) if alignment.anysingle[i] == 0]
    pi_sp_data={}
    for species in spp:
        pi_sp_data[species] = [alignment.species_data[species][i] for i in range(len(alignment.locations)) if alignment.anysingle[i] == 0]

    ALIGNMENTPI.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(len(pi_loc))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
    ALIGNMENTPI.write('[ '+ " ".join(pi_loc)+' ]'+"\n")
    for species in spp: #write sequences for each species
        ALIGNMENTPI.write(species+"\t"+("".join(pi_sp_data[species]))+"\n")
    ALIGNMENTPI.write(';\nend;')
    ALIGNMENTPI.close()

    #Process alignment_pi_singletons.nex
    ALIGNMENTPISINGLETONS=open(fi.replace('.nex','_pi_singletons.nex'),'w')
    pi_singletons_loc = [alignment.locations[i] for i in range(len(alignment.locations)) if (alignment.anysingle[i] == 0 or alignment.othersingle[i] == 1)]
    pi_singletons_sp_data={}
    for species in spp:
        pi_singletons_sp_data[species] = [alignment.species_data[species][i] for i in range(len(alignment.locations)) if (alignment.anysingle[i] == 0 or alignment.othersingle[i] == 1)]

    ALIGNMENTPISINGLETONS.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(len(pi_singletons_loc))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
    ALIGNMENTPISINGLETONS.write('[ '+ " ".join(pi_singletons_loc)+' ]'+"\n")
    for species in spp: #write sequences for each species
        ALIGNMENTPISINGLETONS.write(species+"\t"+("".join(pi_singletons_sp_data[species]))+"\n")
    ALIGNMENTPISINGLETONS.write(';\nend;')
    ALIGNMENTPISINGLETONS.close()

    #Process alignment_singletons.nex
    ALIGNMENTSINGLETONS=open(fi.replace('.nex','_singletons.nex'),'w')
    singletons_loc = [alignment.locations[i] for i in range(len(alignment.locations)) if alignment.singleton[i] == 1]
    singletons_sp_data={}
    for species in spp:
        singletons_sp_data[species] = [alignment.species_data[species][i] for i in range(len(alignment.locations)) if alignment.singleton[i] == 1]

    ALIGNMENTSINGLETONS.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+ntax+' NCHAR='+str(len(singletons_loc))+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
    ALIGNMENTSINGLETONS.write('[ '+ " ".join(singletons_loc)+' ]'+"\n")
    for species in spp: #write sequences for each species
        ALIGNMENTSINGLETONS.write(species+"\t"+("".join(singletons_sp_data[species]))+"\n")
    ALIGNMENTSINGLETONS.write(';\nend;')
    ALIGNMENTSINGLETONS.close()
    
#########################

def main(num_missing, sisrs_dir, composite_dir):

    alignment=get_phy_sites(sisrs_dir,composite_dir,num_missing)
    numbi=alignment.numsnps() #prints numbers of snps, biallelic snps, and singletons
    alignment = write_alignment(sisrs_dir+'/alignment.nex',alignment,numbi)

if __name__ == '__main__':
    num_missing = int(sys.argv[1])
    sisrs_dir = sys.argv[2]
    composite_dir = sys.argv[3]
    main(num_missing, sisrs_dir, composite_dir)
