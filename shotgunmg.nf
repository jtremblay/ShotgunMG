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
                                   |___/                           

                           for  N E X T F L O W

               Support: jtremblay514@gmail.com
             Home page: jtremblay.github.io/pipelines.html
               Version: 1.4.0-beta
###############################################################################

 reads        : ${params.DEFAULT.raw_reads}
 outdir       : ${params.DEFAULT.outdir}
    """.stripIndent()
/*
  Import processes from external files
  It is common to name processes with UPPERCASE strings, to make
  the program more readable
*/
include { TRIMMOMATIC; BBDUK; BBMAP_SUBTRACT; MEGAHIT; METASPADES; PRODIGAL; EXONERATE } from './modules/shotgun_metagenomics'
raw_reads_channel = Channel.fromFilePairs(params.DEFAULT.raw_reads, checkIfExists:true)


workflow {
    /*
        QC
    */
    TRIMMOMATIC(raw_reads_channel)
    BBDUK(TRIMMOMATIC.out.reads)

    // If host/contaminant DNA is to be removed
    if(params.DEFAULT.ref_genome_to_subtract != "null" && params.DEFAULT.ref_genome_to_subtract != ""){
        if(file(params.DEFAULT.ref_genome_to_subtract, checkIfExists: true)){
            //println "valid path or file"
            //println params.DEFAULT.ref_genome_to_subtract
            BBMAP_SUBTRACT(BBDUK.out.reads)
            sample_id = BBMAP_SUBTRACT.out.reads.collect{ it[0] }
            R1 = BBMAP_SUBTRACT.out.reads.collect{ it[1] }
            R2 = BBMAP_SUBTRACT.out.reads.collect{ it[2] }
        }else{
            //Re-group fastqs for assembly.
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
    if(params.DEFAULT.assembler == "megahit"){
        MEGAHIT(sample_id, R1, R2)
        PRODIGAL(MEGAHIT.out.assembly)
        ch = Channel.of("contigs")
        //EXONERATE(ch, MEGAHIT.out.assembly)
    }else if(params.DEFAULT.assembler == "metaspades"){
        METASPADES(sample_id, R1, R2)
        //PRODIGAL(METASPADES.out.assembly)
        //EXONERATE(METASPADES.out.assembly)
    }else{
        exit 1, 'Invalid assembler'
    }
    /*
        Abundance (aligned QCed reads against co-assembly).
    */
    //BEDFILES(PRODIGAL.out.predicted_genes_faa)
    //MAKE_INDEX()
    //ALIGN()
    //GENERATE_ABUNDANCE_MATRICES()

    /*
        Annotations
        see documentation https://jtremblay.github.io/shotgunmg_guide_v1.3.2.html
    */
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

