# **SISRS**  

## Welcome to SISRS! 
SISRS is a bioinformatics pipeline that:  
- Generates genome-scale ortholog datasets  
- Identifies and extracts phylogenetically useful data subsets  
- Requires only whole-genome sequence data  
- Can be run without underlying annotation data or reference genome data (although this always helps).  

Developed in  [2015](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0632-y), this repo contains the most up-to-date scripts, as well as a detailed walkthrough to encourage broad use.  

### How does SISRS work?  
The workhorse of SISRS is the composite genome assembly. Rather than pre-assembling genomes/orthologs for each species and aligning, SISRS pools reads across species and performs a single genome assembly step. The basic steps are:  

1) Basic data preparation (e.g. WGS read trimming, quality assessment)  

2) Based on a user-supplied genome size estimate for the study clade (e.g. 3.5Gb for primates), reads are subset evenly from each species such that the final pooled read depth is ~10X.  

3) Pooled taxa reads are assembled into a single composite genome using Ray.

These contigs represent genomic regions that are conserved enough to be assembled using mixed-read data. All species reads are mapped onto this composite genome, and species-specific bases are called for each position when two key conditions are met:  

- Across all individuals, there is **no intraspecies variation** (i.e. fixed alleles within species, *can be relaxed for exceptionally high-depth datasets*)  
- There are at least **three bases of coverage** for that base (*can be adjusted based on sequencing depth*).  

For each species, any composite bases that don't meet these particular conditions will be denoted as "N". Each position in the composite genome corresponds to putatively othologous base, with deletions denoted by "N", and insertions denoted by "-". In this way, SISRS is able to bypass traditional sequence alignment, directly comparing the bases at each position to enable genome-scale comparisons with reduced computational resources.  

SISRS then performs a few rounds of data partitioning. After all species bases have been called, the following data subsets are automatically extracted:  

1) All variable sites (alignment.nex)  
2) All parsimony-informative sites (alignment_pi.nex; variable sites without singletons)  
3) Parsiomony-informative biallelic sites (**PIBS**; alignment_bi.nex; sites with two detected alleles, and no singletons)  

These sites can be filtered to keep or remove gaps, and also allow different amounts of missing data. If a reference genome and annotation data is available for a focal species (or a very close relative), the data can also be sliced by annotation type.  

PIBS alignments have been shown to perform well for tree building using RAxML and the GTRGAMMA or GTRCAT model, with a correction for ascertainment-bias.  

### Run Requirements and Prerequisite Testing  

The following programs must be installed and in the user path:  
1) Python 3.6+ (with Biopython)  
2) Ray (tested with v2.3.2-devel; Requires working version of mpirun even if using a single node)  
    - Note: In order to use .gz files with Ray, it must be compiled with .gz support.  
```
git clone https://github.com/sebhtml/RayPlatform.git
git clone https://github.com/sebhtml/ray.git
cd ray
make PREFIX=ray-build HAVE_LIBZ=y LDFLAGS="-lpthread -lz"
make install
```
3) Bowtie2 (tested with v2.3.4)  
4) Samtools (tested with v.1.3.1)  
5) BBMap (tested with v.38.73)  
6) BEDTools (if working with annotation data; tested with v.2.26)  
7) FastQC (optional, tested with v.0.11.8)  

If desired, test system-readiness by analyzing the included test data.  
- **Note:** The test script is designed to run using a single processors, and takes ~15m to complete. To increase processors (and decrease run-time), note the sed command below.  

```
git clone https://github.com/BobLiterman/SISRS_Walkthrough.git
cd SISRS_Walkthrough/
# If changing processor counts (i.e use 10 processors)
# sed -i 's/processors 1/processors 10/g' scripts/00_Test_SISRS.sh
bash scripts/00_Test_SISRS.sh
```
This will run a small-scale analysis, including reference genome mapping and annotation slicing.

### Running SISRS With Your Data  

These instructions describe how to run SISRS 'manually', in step-by-step fashion. While a '1-click' version is coming (and will be posted here when completed), many of the steps can actually be run in parallel, and running SISRS in this way allows users to take full use of multinode/HPC-type systems.  

1. **Clone Repo**  
    ```
    git clone https://github.com/BobLiterman/SISRS_Walkthrough.git
    cd SISRS_Walkthrough
      
    # Set folder name if desired
    mv SISRS_Walkthrough <YOUR_ANALYSIS_NAME>
    cd <YOUR_ANALYSIS_NAME>
    ```  
2. **Sample Naming:** For each species/population/OTU in your analysis, add an ID to the TaxonIDs file in the base directory. Avoid spaces or special characters, but for example if analyzing primates, you could use:  
    ```
    Homo_sapiens
    Gorilla_gorilla
    Pan_troglodytes
    Pan_paniscus
    ```
    or  
    ```
    HomSap
    GorGor
    PanTro
    PanPan
    ```

    or  
    ```
    Human
    Gorilla
    Chimp
    Bonobo
    ```  
3. **Folder Setup:** Run the folder setup script, which will set up all the directories needed for a basic SISRS run, based on the names in TaxonIDs.  
    ```
    python scripts/01_sisrs_folder_setup.py
    ```  
4. **Read Trimming:** If you would like to trim your data with our built-in trimming template, put read data for each species into the corresponding subfolder in the raw reads directory (<BASE_DIR>/Reads/RawReads)  
    - **Note:** These files aren't directly altered, so using links (instead of moving or copying raw data) can save space.  
    - **Note:** If your reads are **already trimmed**, those reads can go directly to the trimmed read directory (<BASE_DIR>/Reads/TrimReads)  
    - **Note:** Currently, all reads must be in the **.fastq.gz** format.  
    - **Note:** Single-end and paired-end reads can be added to the folder together, and paired end reads should have a common base name, tailed with **_1.fastq.gz** and **_2.fastq.gz** (e.g HomSap_1.fastq.gz, HomSap_2.fastq.gz)  
    - Interleaved reads will be treated as single-ended, but all SISRS assembly and mapping is single-end based, so this should have limited impact.  
    - **If trimming**, after populating the raw read folders run the trim script:  
    ```
    # Arguments: 
    # -p/--processors: processor count for FastQC (if called) or SLURM script [Default: 1 ]  
    # --skipqc: Don't run FastQC
    # --bash: Don't run trimming, but save a bash script with commands
    # --slurm: Don't run trimming, but save a command script with a SLURM header
      
    # Basic call
    python scripts/02_sisrs_read_trimmer.py
      
    # Output bash script, and skip QC
    python scripts/02_sisrs_read_trimmer.py --bash --skipqc
      
    # Output SLURM script, and request 10 processors
    python scripts/02_sisrs_read_trimmer.py --slurm --processors 10
    ```  
5. **Check data QC:** Putting bad data into SISRS can have a negative impact on site calling, so check the FastQC output paying specific attention to high duplication rates or overrepresented sequences (i.e. non-random sequencing) and very low-quality regions. In many cases (except when sequencing depths are prohibitively low), it would often be better to have less total data that was high-quality by leaving some samples out, as opposed to including more poor-quality data.
6. **Read Subsetting:** The SISRS composite genome is assembled using data sampled evenly among species, and within species, evenly among datasets. The final assembly target depth is ~10X relative to the genome size estimate for the group. To subset your data, run the read subsetting script.  
    ```
    # Arguments:
    # -g/--genomesize (REQUIRED): Genome size estimate for group, in basepairs
      
    # For a 3.5Gb genome (e.g. for primates)...
    python/scripts/03_sisrs_read_subsetter.py -g 3500000000
    ```
7. **Composite Genome Assembly:** After subsetting, SISRS calls Ray to assemble the composite genome.  
    ```
    # Arguments: 
    # -n/--nodes: Number of nodes for Ray [Default: 1]
    # -p/--processors: Number of processors per node [Default: 1]
    # --bash: Don't run Ray, but save a command script
    # --slurm: Don't run Ray, but save a command script with a SLURM header
      
    # Run Ray with 10 nodes, each with 20 processors
    python 04_sisrs_composite_genome_assembly.py --nodes 10 --processors 20
      
    # Make a SLURM script to run Ray, requesting 5 nodes each with 2 processors, 
    python 04_sisrs_composite_genome_assembly.py --nodes 5 --processors 2 --slurm
    ```  
8. **SISRS Run Prep:** After composite genome assembly, the analysis directory (<BASE_DIR>/SISRS_Run) can be populated. This includes indexing the composite genome, and preparing scripts to perform the species-specific read mapping steps which can be run in parallel.  
    ```
    # Arguments:  
    # -p/--processors: Number of available processors per node [Default: 1]
    # -m/--minread: Minimum read coverage to call a base [Default: 3]
    # -t/--threshold: Percent homozygosity required to call a base [Default: 1 (all reads must support a single base)]
    # --slurm (OPTIONAL): Add SLURM header to mapping scripts  
      
    # Prepare scripts for use with 10 processors, with bases only called at 5X coverage, with no intraspecies variation allowed (100% homozygosity threshold)
    python scripts/05_sisrs_run_prep.py --processors 10 --minread 5 --threshold 1     
      
    # Prepare SLURM scripts for use with 20 processors, with bases only called at 50X coverage with 98% homozygosity within species (i.e. allow for rare sequencing errors in high-depth data)
    python scripts/05_sisrs_run_prep.py --processors 20 --minread 50 --threshold 0.98 --slurm 
    ```  
9. **Run Mapping Scripts:** Once the mapping scripts are generated in the analysis directory species folders (<BASE_DIR>/SISRS_Run/<SPECIES>), they can be submitted in a number of ways. These scripts are independent from one another, which is why we leave the option to submit them to the user. If you want to just submit them to run one after another (not recommended for large datasets when parallel processing is available, as run-time is significantly increased), you can use the SISRS run mapping script:  
    ```
    python scripts/06_sisrs_run_mapping.py
    ```
10. **Alignment Output:** Once bases are called for each species, SISRS outputs  alignments. If you want to limit the amount of missing data (i.e. how many species are allowed to have "N" for a given site), you can specify one or more missing data arguments. Filtered alignments will also be split in to data with and without insertions (gaps).  
    ```
    # Arguments:  
    # -m/--missing (OPTIONAL): Number(s) of taxa allowed to have missing data in final alignments; Must be less than total species count - 2  
      
    # Output alignment and filter down to sites with 0, 1, or 2 missing species per site  
    python scripts/07_sisrs_output_alignment.py --missing 0 1 2
    ```  
***
#### The following steps are only relevant when a reference genome is available for one of the focal species (or a close relative)
11. **Reference Genome Mapping:** SISRS ortholog sequences from a user-selected species can be mapped onto a reference genome in order to enable genome-guided site filtering (i.e. only process orthologous that map uniquely in this genome) or generate annotation-specific subsets.  
      - **Note:** If the reference genome FASTA **has already been indexed by Bowtie2**, link or copy the *bt2 Bowtie2 index files to the reference genome directory (<BASE_DIR>/Reference_Genome).
      - **Note:** If the reference genome FASTA **has not been indexed by Bowtie2**, link or copy the FASTA to the reference genome directory (<BASE_DIR>/Reference_Genome) and run:  
      ```
      bowtie2-build <REFERENCE_FASTA> <USER_INDEX_ID> -p <PROCESSORS>
      ``` 
      Then, you can map the orthologs, extract uniquely mapping contigs, and isolate non-overlapping sites.
    ```
    # Arguments:
    # -p/--processors: Number of available processors per node [Default: 1]
    # -r/--reference: Bowtie2 index name
    # -t/--taxon: SISRS taxon ID for sample to be mapped
    # --bash: Don't run, but write command to BASH script (Reference_Genome_Mapping.sh) in 'scripts' folder
    # --slurm: Don't run, and add SLURM header to mapping script
      
    # Map orthologs from human (HomSap) onto human genome (Hg19) using 10 processors 
    python scripts/08_sisrs_reference_mapping.py --processors 10 --reference Hg19 --taxon HomSap
    ```  
12. **Alignment Slicing**: Alignments generated by SISRS can be sliced down to a specified subset of sites using the following script:  
    ```
    # Arguments:
    # -a/--alignment: Path to alignment to be sliced
    # -l/--locs: Path to LocList for original alignment (typically generated when filtering SISRS alignments for missing data)
    # -s/--subsetlocs:  Path to  LocList of desired sites to extract 
    # -n/--name: An ID (file prefix) for the slice subset 
    # -o/--outputdir: Where to save sliced sites [Default: Directory where slice locs are located]
      
    # Slice biallelic alignment with 0 missing taxa down to just those sites that uniquely map to the human genome
    python scripts/09_sisrs_alignment_slicer.py --alignment <BASE_DIR>/SISRS_Run/alignment_bi_m0_nogap.phylip-relaxed --locs <BASE_DIR>/SISRS_Run/alignment_bi_locs_m0_nogap.txt --subsetlocs <BASE_DIR>/Reference_Genome/Mapping/HomSap_NonDup_Genome_Sites --name Hg19_Mapped --outputdir ../Reference_Genome/Hg19_Mapped
    ```
13. **Annotation Slicing:** If you want to subset alignments by annotation types from the reference genome, add/link BED files for each desired annotation subset into the annotation directory (<BASE_DIR>/Reference_Genome/Annotation).  
    - **Note:** This script will check to see if an annotation has already been processed for a given annotation/dataset. If you want to reprocess a BED file/alignment combination (e.g. if you've updated either file), enable overwriting with --overwrite. If you want to ignore certain BED files, use the --ignore flag
    ```
    # Arguments: 
    # -a/--alignment: Path to alignment to be sliced
    # -l/--locs: Path to LocList for alignment to be sliced
    # -i/--ignore (OPTIONAL): The basenames of one or more annotation files to ignore (for example -i CDS Gene)
    # Arguments: -n/--name: Prefix for output alignments [Default: Basename of alignment file]
    # Arguments: -ow/--overwrite: Process all BED files in Reference_Genome/Annotation, overwriting existing results if present [Default: Process only unprocessed BED files, don't overwrite existing data]
      
    # Slice SISRS alignment (alignment_bi_m2.phylip-relaxed) down to just sites from CDS or Gene annotations (if CDS.bed and Gene.bed are in annotation directory)
    python scripts/10_sisrs_annotation_slicer.py --alignment <BASE_DIR>/SISRS_Run/alignment_bi_m2_nogap.phylip-relaxed --locs <BASE_DIR>/SISRS_Run/alignment_bi_locs_m2_nogap.txt --name Biallelic_2Missing
      
    # Slice SISRS alignment (alignment_bi_m2.phylip-relaxed) down to just sites from CDS (even though Gene.bed is in annotation directory)
    python scripts/10_sisrs_annotation_slicer.py --alignment <BASE_DIR>/SISRS_Run/alignment_bi_m2_nogap.phylip-relaxed --locs <BASE_DIR>/SISRS_Run/alignment_bi_locs_m2_nogap.txt --name Biallelic_2Missing --ignore Gene

    ``` 