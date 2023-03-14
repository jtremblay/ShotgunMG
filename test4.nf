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
include { TRIMMOMATIC; BBDUK; BBMAP_SUBTRACT; MEGAHIT; METASPADES; PRODIGAL; EXONERATE_CONTIGS; EXONERATE_GENES; BEDFILE_CONTIGS; BEDFILE_GENES; MAKE_BWA_INDEX; BWAMEM_PE } from './modules/shotgun_metagenomics'
raw_reads_channel = Channel.fromFilePairs(params.DEFAULT.raw_reads, checkIfExists:true)

workflow {
    /*
        QC
    */
    TRIMMOMATIC(raw_reads_channel)
    BBDUK(TRIMMOMATIC.out.reads)

    // If host/contaminant DNA is to be removed
    ch_qced_reads = Channel.empty()
    if(params.DEFAULT.ref_genome_to_subtract != "null" && params.DEFAULT.ref_genome_to_subtract != ""){
        if(file(params.DEFAULT.ref_genome_to_subtract, checkIfExists: true)){
            //println "valid path or file"
            //println params.DEFAULT.ref_genome_to_subtract
            BBMAP_SUBTRACT(BBDUK.out.reads)
            tmp_channel = BBMAP_SUBTRACT.out.reads.collect().map { sample_id, reads -> [ it[1], it[2] ] }
            tmp_channel.view()
            ch_qced_reads = ch_qced_reads.mix(tmp_channel)
            //ch_qced_reads = ch_qced_reads.mix(BBMAP_SUBTRACT.out.reads.collect().ifEmpty([]))
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
    //ch_assembly.collect().view()
    //ch_qced_reads.view()
    //raw_reads_channel.view()
    //BWAMEM_PE(ch_assembly.collect(), ch_qced_reads)
    //GENERATE_ABUNDANCE_MATRICES()

    /*
        Annotations
        see documentation https://jtremblay.github.io/shotgunmg_guide_v1.3.2.html
    */
    //EXONERATE_CONTIGS(ch_assembly.collect())
    //EXONERATE_GENES(PRODIGAL.out.predicted_genes_faa)
    //DIAMOND_BLASTNR()
    //HMMSEARCH_KEGG()
    //HMMSEARCH_PFAM()
    //CAT()
    
    /*
        Microbial ecology metrics
    */
    //ALPHA_DIVERSITY()
    //BETA_DIVERSITY()
    //TAXONOMIC_PROFILES()

    /*
        MAGs
    */
    //GENERATE_MAGS()


}

