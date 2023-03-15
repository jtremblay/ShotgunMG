# ShotgunMG

This repository contains an implementation of the ShotgunMG pipeline (https://doi.org/10.1093/bib/bbac443) for Nextflow. The original pipeline implemented with the GenPipes workflow management system is available here: https://bitbucket.org/jtremblay514/nrc_pipeline_public. 

All the modules defined in the ```shotgunmg.config``` file should be installed and functional. The nrc_tools utilities can be found here: https://bitbucket.org/jtremblay514/nrc_tools_public . Briefly, this pipeline takes a set of raw reads (i.e. short Illumina reads), performs quality control and co-assemble the QC-controlled reads. These reads are then mapped against the co-assembly to generate contig and gene abundance matrices. The co-assembly is also processed through a gene caller (i.e. Prodigal). Resulting genes are functionally annotated using hmmsearch vs pfam; hmmsearch vs kofam and rpsblast vs COG. Taxonomic annotations are assigned using the CAT package. Finally, MAGs are generated using MetaBAT2. End results are included in the ```output/``` directory. Ultimately, this pipeline processes raw fastqs into gene and contig abundance matrices and functional and taxonomic annotation files.

A replicated simple mock community dataset is available here https://doi.org/10.5281/zenodo.7140751 and is a good dataset to test this pipeline. A fully functional implementation of the pipeline is available as a Docker image: https://cloud.docker.com/u/julio514/repository/docker/julio514/centos

This project is in development - more coming soon. 

```
###############################################################################
              _____ _           _                    __  __  _____ 
             / ____| |         | |                  |  \/  |/ ____|  
            | (___ | |__   ___ | |_ __ _ _   _ _ __ | \  / | |  __ 
             \___ \| '_ \ / _ \| __/ _` | | | | '_ \| |\/| | | |_ |
             ____) | | | | (_) | || (_| | |_| | | | | |  | | |__| |
            |_____/|_| |_|\___/ \__\__, |\__,_|_| |_|_|  |_|\_____|
                                    __/ |                          
                                   |___/    for N E X T F L O W 

                Github: https://github.com/jtremblay/ShotgunMG
             Home page: jtremblay.github.io/pipelines.html
               Version: 1.4.0-beta                                         
###############################################################################
```

![Figure 1](./images/overview_shotgunmg.png)
*Figure 1. Overview of ShotgunMG.
1) Reads of each library are controlled for quality. 2) Quality controlled reads are co-assembled into one single de novo assembly. Gene coordinates are computed on each contig. Quality controlled reads are mapped on the co-assembly to estimate contig and gene abundance. 3) Contig and gene abundance are summarized into abundance matrices where columns = samples/liraries and rows = contig or gene identifiers. 4) Genes are annotated for taxonomy and functions and compiled on one single database (5). 6) These end results can then be used for downstream analyses.*

## Usage:
Once Nextflow (and an appropriate version of Java) is installed, you can run the pipeline like this:
```
nextflow run -c ./shotgunmg.config ./test5.nf -resume
```

In the ```shotgun.config``` file are all the parameters used for every steps of the pipeline. There you can customize the amounts of resources of each step, depending on this size and complexity of the dataset to analyze. The pipeline relies on environment modules (https://modules.readthedocs.io/en/latest/) which means that each software required by the pipeline have to be available through a module. For instance, for the co-assembly step, the MEGAHIT should be made available by first loading the module : (i.e. ```module load nrc/megahit/1.2.9```) and then running the software (i.e. ```megahit -h```).

## Databases
The pipeline relies on many databases in order to run the various annotations.

### CAT
Go here - https://github.com/dutilh/CAT - and follow the instructions under the preconstructed databases section,

### CheckM
Follow the instructions : https://github.com/Ecogenomics/CheckM

### COG
https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Cog_LE.tar.gz

### KOG 
https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Kog_LE.tar.gz

### NCBI nr
ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz (Do not forget to run ```diamond makedb --in nr -p 4 --db nr.dmnd``` once downloaded).

### PFAM-A
Available here: http://ftp.ebi.ac.uk/pub/databases/Pfam. Use the latest version. May have to run ```hmmpress``` once downloaded.

### Contaminants
Available here: http://jtremblay.github.io/files/contaminants.tar.gz. Contains known Illumina contaminants, other sequencing artefacts and adapter sequences from  various kits (NexteraXT, TruSeq, etc.). In fasta format. 

### KEGG
KEGG orthologs (KO) assignment is done using kofamscan available here: https://www.genome.jp/ftp/tools/kofam_scan/kofam_scan-1.3.0.tar.gz. In order to run kofamscan, you will need to have the HMM profiles of each KO - available here: https://www.genome.jp/ftp/db/kofam/profiles.tar.gz - and also the KO link file - available here: https://www.genome.jp/ftp/db/kofam/ko_list.gz. Another file needed to link each KO to their associated pathway and/or module is available here: https://https://jtremblay.github.io/files/kegg_ref_pathways_modules_combined.tsv.gz.
Once downloaded, uncompress these files and move them into their location which should be $INSTALL_HOME/databases/kofam/<date>/. Double check their path in the .ini file under the [kofamscan] and [parse_kofam] sections.
KOfamscan generates lots of intermediate files which can actually make it impossible to use for large metagenomes. To circumvent this issue, you can concatenate all kofam individual profiles (```cat ./profiles/K* > kofam.hmm``` and run hmmpress ```hmmpress kofam.hmm```). Once done this kofam.hmm file can be used with the ```hmmsearch``` software - see KEGG section downstream. Our internal benchmarks showed that hmmsearch and kofamscan gave identical results (i.e. for both methods, each gene pointed to the same KO).


## Setting up files needed by the pipeline
The first step in running ShotgunMG is to setup the files required by the pipeline to run. Fastq libraries usually come in the form of demultiplexed paired end sequencing libraries - one library per sample. These ```.fastq.gz``` files should be stored in a directory labeled ```raw_reads/```.


## Diagram of the pipeline
Here is the Mermaid diagram of the pipeline.

```mermaid
flowchart TD
    p0((Channel.fromFilePairs))
    p1[TRIMMOMATIC]
    p2(( ))
    p3(( ))
    p4[BBDUK]
    p5(( ))
    p6(( ))
    p7((Channel.empty))
    p8(( ))
    p9[BBMAP_SUBTRACT]
    p10(( ))
    p11([collect])
    p12([collect])
    p13([collect])
    p14((Channel.empty))
    p15[MEGAHIT]
    p16(( ))
    p17(( ))
    p18([collect])
    p19([ifEmpty])
    p20([mix])
    p21([collect])
    p22[PRODIGAL]
    p23(( ))
    p24([collect])
    p25[BEDFILE_CONTIGS]
    p26[BEDFILE_GENES]
    p27([collect])
    p28[MAKE_BWA_INDEX]
    p29([collect])
    p30[BWAMEM_PE]
    p31(( ))
    p32[BEDTOOLS_COV_CONTIGS]
    p33[BEDTOOLS_COV_GENES]
    p34([map])
    p35([collect])
    p36[MERGE_COV_CONTIGS]
    p37([map])
    p38([collect])
    p39[MERGE_COV_GENES]
    p40([collect])
    p41[EXONERATE_CONTIGS]
    p42(( ))
    p43(( ))
    p44(( ))
    p45[EXONERATE_GENES]
    p46(( ))
    p47(( ))
    p48(( ))
    p49(( ))
    p50([map])
    p51([transpose])
    p52[DIAMOND_BLASTP_NR]
    p53([collect])
    p54[MERGE_DIAMOND_BLASTP_NR]
    p55[HMMSEARCH_PFAM]
    p56([collect])
    p57([collect])
    p58[MERGE_PFAM]
    p59(( ))
    p60(( ))
    p61[RPSBLAST_COG]
    p62([collect])
    p63[MERGE_COG]
    p64[COG_OVERREP]
    p65(( ))
    p66[HMMSEARCH_KEGG]
    p67([collect])
    p68([collect])
    p69[MERGE_KEGG]
    p70(( ))
    p71[PARSE_KEGG]
    p72[KO_OVERREP]
    p73(( ))
    p74[CONVERT_IDS_FOR_CAT]
    p75[CAT]
    p76[GENERATE_FEATURE_TABLES]
    p77(( ))
    p78(( ))
    p79[SUMMARIZE_TAXONOMY]
    p80(( ))
    p81(( ))
    p82[BETA_DIVERSITY_BACTARCH]
    p83(( ))
    p84(( ))
    p85(( ))
    p86[BETA_DIVERSITY_ALL]
    p87(( ))
    p88(( ))
    p89(( ))
    p90[ALPHA_DIVERSITY_CONTIGS]
    p91(( ))
    p92[ALPHA_DIVERSITY_GENES]
    p93(( ))
    p94[COG_MATRIX_RPOB]
    p95[COG_MATRIX_RECA]
    p96[ALPHA_DIVERSITY_RPOB]
    p97(( ))
    p98[ALPHA_DIVERSITY_RECA]
    p99(( ))
    p0 -->|raw_reads_channel| p1
    p1 --> p4
    p1 --> p3
    p1 --> p2
    p4 --> p9
    p4 --> p6
    p4 --> p5
    p7 -->|ch_qced_reads| p8
    p9 --> p11
    p9 --> p10
    p11 -->|sample_id| p15
    p9 --> p12
    p12 -->|R1| p15
    p9 --> p13
    p13 -->|R2| p15
    p14 -->|ch_assembly| p20
    p15 --> p18
    p15 --> p17
    p15 --> p16
    p18 --> p19
    p19 --> p20
    p20 -->|ch_assembly| p21
    p21 --> p22
    p22 --> p23
    p22 --> p45
    p22 --> p26
    p22 --> p75
    p20 -->|ch_assembly| p24
    p24 --> p25
    p25 --> p32
    p26 --> p33
    p20 -->|ch_assembly| p27
    p27 --> p28
    p28 --> p30
    p28 --> p30
    p28 --> p30
    p28 --> p30
    p28 --> p30
    p28 --> p30
    p20 -->|ch_assembly| p29
    p29 --> p30
    p9 --> p30
    p30 --> p32
    p30 --> p31
    p32 --> p34
    p30 --> p33
    p33 --> p37
    p34 --> p35
    p35 --> p36
    p36 --> p76
    p37 --> p38
    p38 --> p39
    p39 --> p64
    p20 -->|ch_assembly| p40
    p40 --> p41
    p41 --> p44
    p41 --> p43
    p41 --> p42
    p45 --> p49
    p45 --> p48
    p45 --> p47
    p45 -->|NUMCHUNKS| p46
    p45 --> p50
    p50 --> p51
    p51 -->|ch_gene_chunks| p52
    p52 --> p53
    p53 --> p54
    p54 --> p74
    p51 -->|ch_gene_chunks| p55
    p55 --> p57
    p55 --> p56
    p56 --> p58
    p57 --> p58
    p58 --> p60
    p58 --> p59
    p51 -->|ch_gene_chunks| p61
    p61 --> p62
    p62 --> p63
    p63 --> p64
    p64 --> p65
    p51 -->|ch_gene_chunks| p66
    p66 --> p68
    p66 --> p67
    p67 --> p69
    p68 --> p69
    p69 --> p71
    p69 --> p70
    p71 --> p72
    p39 -->|gene_abundance| p72
    p72 --> p73
    p22 -->|infile_gff| p74
    p74 --> p75
    p15 -->|contigs_fna| p75
    p75 --> p76
    p76 --> p78
    p76 --> p79
    p76 --> p82
    p76 --> p77
    p79 --> p81
    p79 --> p80
    p82 --> p85
    p82 --> p84
    p82 --> p83
    p76 -->|feature_table| p86
    p86 --> p89
    p86 --> p88
    p86 --> p87
    p36 -->|contig_abundance| p90
    p90 --> p91
    p39 -->|gene_abundance| p92
    p92 --> p93
    p39 -->|gene_abundance| p94
    p63 -->|rpsblast_cog| p94
    p94 --> p96
    p39 -->|gene_abundance| p95
    p63 -->|rpsblast_cog| p95
    p95 --> p98
    p96 --> p97
    p98 --> p99
```
