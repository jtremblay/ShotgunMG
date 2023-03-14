# ShotgunMG

This repository contains an implementation of the ShotgunMG pipeline (https://doi.org/10.1093/bib/bbac443) for Nextflow. The original pipeline implemented with the GenPipes workflow management system is available here: https://bitbucket.org/jtremblay514/nrc_pipeline_public. 

All the modules defined in the ```shotgunmg.config``` file should be installed and functional. The nrc_tools utilities can be found here: https://bitbucket.org/jtremblay514/nrc_tools_public . Briefly, this pipeline takes a set of raw reads (i.e. short Illumina reads), performs quality control and co-assemble the QC-controlled reads. These reads are then mapped against the assembly to generate contigs and genes abundance matrices. The assembly is also processed through a gene caller (i.e. Prodigal). Resulting genes are functionally annotated using hmmsearch vs pfam; hmmsearch vs kofam and rpsblast vs COG. Taxonomic annotations are assigned using the CAT package. 

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
