//params.raw_reads = "$projectDir/raw_reads/*_R{1,2}.fastq.gz"
//params.outdir = "$projectDir/output/"

//raw_reads = $params.DEFAULT.raw_reads_dir/$params.DEFAULT.raw_reads_pattern
//outdir = ${params.DEFAULT.outdir}

/*
    parameters are defined in the shotgunmg.config file
*/


log.info """
###############################################################################
              _____ _           _                    __  __  _____ 
             / ____| |         | |                  |  \\/  |/ ____|
            | (___ | |__   ___ | |_ __ _ _   _ _ __ | \\  / | |  __ 
             \\___ \\| '_ \\ / _ \\| __/ _` | | | | '_ \\| |\\/| | | |_ |
             ____) | | | | (_) | || (_| | |_| | | | | |  | | |__| |
            |_____/|_| |_|\\___/ \\__\\__, |\\__,_|_| |_|_|  |_|\\_____|
                                    __/ |                          
                                   |___/     for  N E X T F L O W 

               Support: jtremblay514@gmail.com
             Home page: jtremblay.github.io/pipelines.html
               Version: 1.4.0-beta
###############################################################################

 reads   : ${params.DEFAULT.raw_reads}
 outdir  : ${params.DEFAULT.outdir}
    """.stripIndent()
/*
  Import processes from external files
  It is common to name processes with UPPERCASE strings, to make
  the program more readable
*/
include { TRIMMOMATIC; BBDUK; BBMAP_SUBTRACT; MEGAHIT; METASPADES; PRODIGAL; EXONERATE_CONTIGS; EXONERATE_GENES; BEDFILE_CONTIGS; BEDFILE_GENES; MAKE_BWA_INDEX; BWAMEM_PE; BEDTOOLS_COV_CONTIGS; BEDTOOLS_COV_GENES; MERGE_COV_GENES; MERGE_COV_CONTIGS; DIAMOND_BLASTP_NR; MERGE_DIAMOND_BLASTP_NR; HMMSEARCH_KEGG; MERGE_KEGG; PARSE_KEGG; KO_OVERREP; HMMSEARCH_PFAM; MERGE_PFAM; RPSBLAST_COG; MERGE_COG; COG_OVERREP; CONVERT_IDS_FOR_CAT; CAT; GENERATE_FEATURE_TABLES; SUMMARIZE_TAXONOMY; BETA_DIVERSITY_BACTARCH; BETA_DIVERSITY_ALL; ALPHA_DIVERSITY_CONTIGS; ALPHA_DIVERSITY_GENES; COG_MATRIX_RPOB; COG_MATRIX_RECA; ALPHA_DIVERSITY_RPOB; ALPHA_DIVERSITY_RECA } from './modules/shotgun_metagenomics'
raw_reads_channel = Channel.fromFilePairs(params.DEFAULT.raw_reads, checkIfExists:true)

workflow {
    /*
        QC
    */
    TRIMMOMATIC(raw_reads_channel)
    BBDUK(TRIMMOMATIC.out.reads)

    // If host/contaminant DNA is to be removed
    ch_qced_reads = Channel.empty()
    qced_reads = []
    if(params.DEFAULT.ref_genome_to_subtract != "null" && params.DEFAULT.ref_genome_to_subtract != ""){
        if(file(params.DEFAULT.ref_genome_to_subtract, checkIfExists: true)){
            BBMAP_SUBTRACT(BBDUK.out.reads)
            //Re-group fastqs for assembly.
            sample_id = BBMAP_SUBTRACT.out.reads.collect{ it[0] }
            R1 = BBMAP_SUBTRACT.out.reads.collect{ it[1] }
            R2 = BBMAP_SUBTRACT.out.reads.collect{ it[2] }
        }else{
            //Re-group fastqs for assembly.
            ch_qced_reads = ch_qced_reads.mix(BBDDUK.out.reads.collect().ifEmpty([]))
            sample_id = BBDUK.out.reads.collect{ it[0] }
            R1 = BBDUK.out.reads.collect{ it[1] }
            R2 = BBDUK.out.reads.collect{ it[2] }
        }
    }else{
        exit 1, 'Invalid path'
    }
    /*
        Co-assembly, gene prediction and split faa and fna in small chunks for annotations,.
    */
    ch_assembly = Channel.empty()
    if(params.DEFAULT.assembler == "megahit"){
        MEGAHIT(sample_id, R1, R2)
        ch_assembly = ch_assembly.mix(MEGAHIT.out.assembly.collect().ifEmpty([]))

    }else if(params.DEFAULT.assembler == "metaspades"){
        METASPADES(sample_id, R1, R2)
        ch_assembly = ch_assembly.mix(SPADES.out.assembly.collect().ifEmpty([]))

    }else{
        exit 1, 'Invalid assembler'
    }

    /*
        Predict CDS of genes
    */
    PRODIGAL(ch_assembly.collect())

    /*
        Abundance (aligned QCed reads against co-assembly).
    */
    BEDFILE_CONTIGS(ch_assembly.collect())
    BEDFILE_GENES(PRODIGAL.out.predicted_genes_gff)
    MAKE_BWA_INDEX(ch_assembly.collect())
    BWAMEM_PE(
        MAKE_BWA_INDEX.out.amb, MAKE_BWA_INDEX.out.ann, MAKE_BWA_INDEX.out.bwt, MAKE_BWA_INDEX.out.fai, MAKE_BWA_INDEX.out.pac, MAKE_BWA_INDEX.out.sa, 
        ch_assembly.collect(), 
        BBMAP_SUBTRACT.out.reads
    )
    //GENERATE_ABUNDANCE_MATRICES()
    BEDTOOLS_COV_CONTIGS(BWAMEM_PE.out.bam, BEDFILE_CONTIGS.out.bed_contigs)
    BEDTOOLS_COV_GENES(BWAMEM_PE.out.bam, BEDFILE_GENES.out.bed_genes)
    
    MERGE_COV_CONTIGS(
        BEDTOOLS_COV_CONTIGS.out.cov
            .map{sample_id, file -> file}.collect()
    )
    MERGE_COV_GENES(
        BEDTOOLS_COV_GENES.out.cov
            .map{sample_id, file -> file}.collect()
    )

    /*
        Annotations
        see documentation https://jtremblay.github.io/shotgunmg_guide_v1.3.2.html
        Contig and gene files are first split in smaller chunks, which are
        then annotated using diamond_blastp, hmmsearch. The generated  output files
        are then merged together and filtered.
    */
    EXONERATE_CONTIGS(ch_assembly.collect())
    EXONERATE_GENES(PRODIGAL.out.predicted_genes_faa)
    /* ugh this one gave me a hard time. Briefly I had to figure a way to give each
       output file name a unique identifier, because downstream, we'll have to merge the 
       blast output files that were split by exonerate. But this can't be done if all 
       files have the same name, becuase Nextflow will try to create many symlinks of 
       the same name in the same directory which will return an error... transpose() was
       the key here.
    */
    ch_gene_chunks = EXONERATE_GENES.out.outfiles
        .map { it -> [it.name, it] }
        .transpose()
    //ch_gene_chunks.view()
    //prefixes = ch_tmp.collect{ it[0] }
    //files = ch_tmp.collect{ it[1] }
    //prefixes.view()
    //files.view()

    /* NCBI NR. Will be used later in CAT (taxonomy assignation) */
    DIAMOND_BLASTP_NR(ch_gene_chunks)
    MERGE_DIAMOND_BLASTP_NR(DIAMOND_BLASTP_NR.out.diamond_blastp_outfiles.collect())
    
    /* PFAM */
    HMMSEARCH_PFAM(ch_gene_chunks)
    MERGE_PFAM(
        HMMSEARCH_PFAM.out.tblout_outfiles.collect(),
        HMMSEARCH_PFAM.out.domtblout_outfiles.collect()
    )
    //TODO PFAM_OVERREP()

    /* COG */
    RPSBLAST_COG(ch_gene_chunks)
    MERGE_COG(RPSBLAST_COG.out.outfiles.collect())
    COG_OVERREP(
        MERGE_COG.out.tsv,
        MERGE_COV_GENES.out.gene_abundance_matrix
    )

    /* KEGG KO */
    HMMSEARCH_KEGG(ch_gene_chunks)
    MERGE_KEGG(
        HMMSEARCH_KEGG.out.tblout_outfiles.collect(), 
        HMMSEARCH_KEGG.out.domtblout_outfiles.collect()
    )
    PARSE_KEGG(MERGE_KEGG.out.tblout)
    KO_OVERREP(
        PARSE_KEGG.out.KOs_parsed, 
        MERGE_COV_GENES.out.gene_abundance_matrix
    )

    /*
        Taxonomy annotation
    */
    CONVERT_IDS_FOR_CAT(PRODIGAL.out.predicted_genes_gff, MERGE_DIAMOND_BLASTP_NR.out.tsv)
    CAT(
        MEGAHIT.out.assembly, 
        PRODIGAL.out.predicted_genes_raw_faa, 
        CONVERT_IDS_FOR_CAT.out.orf_ids
    )
    GENERATE_FEATURE_TABLES(
        CAT.out.classification_with_names, 
        MERGE_COV_CONTIGS.out.contig_abundance_matrix
    )
    /*
        Microbial ecology metrics
    */
    SUMMARIZE_TAXONOMY(GENERATE_FEATURE_TABLES.out.feature_table)
    BETA_DIVERSITY_BACTARCH(GENERATE_FEATURE_TABLES.out.feature_table_ba)
    BETA_DIVERSITY_ALL(GENERATE_FEATURE_TABLES.out.feature_table)
    ALPHA_DIVERSITY_CONTIGS(MERGE_COV_CONTIGS.out.contig_abundance_matrix)
    ALPHA_DIVERSITY_GENES(MERGE_COV_GENES.out.gene_abundance_matrix)
    COG_MATRIX_RPOB(MERGE_COV_GENES.out.gene_abundance_matrix, MERGE_COG.out.tsv)
    COG_MATRIX_RECA(MERGE_COV_GENES.out.gene_abundance_matrix, MERGE_COG.out.tsv)
    ALPHA_DIVERSITY_RPOB(COG_MATRIX_RPOB.out.rpob_abundance)
    ALPHA_DIVERSITY_RECA(COG_MATRIX_RECA.out.reca_abundance)

    /*
        MAGs
    */
    //GENERATE_MAGS()


}

