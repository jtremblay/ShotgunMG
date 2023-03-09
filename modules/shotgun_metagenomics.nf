/*
    Shotgun metagenomics functions 
*/

process TRIMMOMATIC {
    publishDir "$params.DEFAULT.outdir/qced_reads/", mode: 'symlink'
    debug true
    cpus params.fake_trimmomatic.cluster_cpus
    memory params.fake_trimmomatic.cluster_memory
    time params.fake_trimmomatic.cluster_time


    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_paired_R1.fastq.gz"), path("${sample_id}_paired_R2.fastq.gz"), emit: reads 
        path("${sample_id}_trimmomatic_log.txt"), emit: log
        path("${sample_id}_trimmomatic_stats.tsv"), emit: stats

    script:
        command = """
        module load ${params.modules.trimmomatic} $params.modules.java && \\
        java -XX:ParallelGCThreads=$params.trimmomatic.threads -Xmx2G -jar \$TRIMMOMATIC_JAR PE \\
            -threads ${params.trimmomatic.threads} \\
            -phred${params.trimmomatic.quality_offset} \\
            ${reads[0]} ${reads[1]} \\
            ${sample_id}_paired_R1.fastq.gz ${sample_id}_single_R1.fastq.gz ${sample_id}_paired_R2.fastq.gz ${sample_id}_single_R2.fastq.gz \\
            ILLUMINACLIP:${params.trimmomatic.adapter_fasta}${params.trimmomatic.illumina_clip_settings} \\
            TRAILING:${params.trimmomatic.trailing_min_quality} \\
            SLIDINGWINDOW:${params.trimmomatic.sliding_window1}:${params.trimmomatic.sliding_window2} \\
            MINLEN:${params.trimmomatic.min_length} \\
            HEADCROP:${params.trimmomatic.headcrop}"""
        if(params.trimmomatic.crop){
            command += """ CROP:${params.trimmomatic.crop}"""
        }
        command += """ 2> ${sample_id}_trimmomatic_log.txt && \\
        grep ^Input ${sample_id}_trimmomatic_log.txt | \\
        perl -pe 's/^Input Read Pairs: (\\d+).*Both Surviving: (\\d+).*Forward Only Surviving: (\\d+).*\$/Raw Fragments,\\1\\nFragment Surviving,\\2\\nSingle Surviving,\\3/' >  ${sample_id}_trimmomatic_stats.tsv
        """
        return command

}

process BBDUK {
    publishDir "$params.DEFAULT.outdir/qced_reads/", mode: 'symlink'
    debug true
    cpus params.fake_bbduk.cluster_cpus
    memory params.fake_bbduk.cluster_memory
    time params.fake_bbduk.cluster_time

    //carefull with tuple. If only one variable, dont include the damn keywork tuple!!!
    input:
        tuple val(sample_id), path(infile_R1), path(infile_R2)

    output:
        tuple val(sample_id), path("${sample_id}_paired_ncontam_R1.fastq.gz"), path("${sample_id}_paired_ncontam_R2.fastq.gz"), emit: reads
        tuple val(sample_id), path("${sample_id}_paired_contam_R1.fastq.gz"), path("${sample_id}_paired_contam_R2.fastq.gz"), emit: bad_reads
        path("${sample_id}_bbduk_log.txt"), emit: log

    script:
        """
        module load ${params.modules.bbmap} ${params.modules.java} && \\
        bbduk.sh \\
            in=${infile_R1} \\
            in2=${infile_R2} \\
            out=${sample_id}_paired_ncontam_R1.fastq.gz \\
            out2=${sample_id}_paired_ncontam_R2.fastq.gz \\
            outm=${sample_id}_paired_contam_R1.fastq.gz \\
            outm2=${sample_id}_paired_contam_R2.fastq.gz \\
            stats=${sample_id}_bbduk_log.txt \\
            k=${params.bbduk.k} \\
            minkmerhits=${params.bbduk.c} \\
            ref=${params.DEFAULT.contaminants} \\
            overwrite=true \\
            threads=1
        """
}

process BBMAP_SUBTRACT {
    publishDir "$params.DEFAULT.outdir/qced_reads/", mode: 'symlink'
    debug true
    cpus params.bbmap_sub.cluster_cpus
    memory params.bbmap_sub.cluster_memory
    time params.bbmap_sub.cluster_time

    //carefull with tuple. If only one variable, dont include the damn keywork tuple!!!
    input:
        tuple val(sample_id), path(infile_R1), path(infile_R2)

    output:
        tuple val(sample_id), path("${sample_id}_paired_ncontam_sub_R1.fastq.gz"), path("${sample_id}_paired_ncontam_sub_R2.fastq.gz"), emit: reads
        path("${sample_id}_bbduksub_log.txt"), emit: log

    script:
        """
        module load ${params.modules.bbmap} ${params.modules.java} && \\
        bbmap.sh \\
            -threads=${params.bbmap_sub.num_threads} \\
            -Xmx${params.bbmap_sub.ram} \\
            in=${infile_R1} in2=${infile_R2} \\
            outu=${sample_id}_paired_ncontam_sub_R1.fastq.gz outu2=${sample_id}_paired_ncontam_sub_R2.fastq.gz \\
            minid=${params.bbmap_sub.min_id} \\
            sam=1.3 nhtag=t mdtag=t xmtag=t amtag=t nmtag=t xstag=us \\
            maxindel=3 bwr=0.16 bw=12 fast=t overwrite=t minhits=2 \\
            path=${params.DEFAULT.ref_genome_to_subtract} \\
            statsfile=${sample_id}_bbduksub_log.txt
        """
}

process MEGAHIT {
    publishDir "$params.DEFAULT.outdir/assembly/", mode: 'symlink'
    debug true
    cpus params.megahit.cluster_cpus
    memory params.megahit.cluster_memory
    time params.megahit.cluster_time
    
    input:
        val(sample_id)
        path(reads1)
        path(reads2)

    output:
        path("megahit/final.contigs.fa"), emit: assembly
        path("megahit/log"), emit: log
        path("megahit/assembly_stats.txt"), emit: stats

    script:
        def input_R1 = "-1 " + reads1.join(",")
        def input_R2 = "-2 " + reads2.join(",")
        def requested_memory_in_bytes = params.megahit.cluster_memory.toBytes()

        """
        module load ${params.modules.megahit} ${params.modules.tools} && \\
        rm -rf megahit -rf && \\
        megahit -t ${params.megahit.num_threads} --k-min ${params.megahit.kmin} \\
            --k-max ${params.megahit.kmax} --k-step ${params.megahit.kstep} \\
            --min-contig-len ${params.megahit.min_contig_length} \\
            ${input_R1} \\
            ${input_R2} \\
            --memory ${requested_memory_in_bytes} \\
            --out-dir megahit && \\
            compileAssemblyResultsSingle.pl --infile megahit/final.contigs.fa > megahit/assembly_stats.txt
        """
}

process METASPADES {
    publishDir "$params.DEFAULT.outdir/assembly/", mode: 'symlink'
    debug true
    
    input:
        val(sample_id)
        path(reads1)
        path(reads2)

    output:
        path("spades/scaffolds_gt${params.spades.length}.fasta"), emit: assembly
        path("spades/assembly_stats.txt"), emit: stats
        //path("log.txt"), emit: log

    script:
        def input_R1 = " -1 " + reads1.join(" -1 ")
        def input_R2 = " -2 " + reads2.join(" -2 ")
        """
        module load ${params.modules.spades} ${params.modules.tools} && \\
        rm spades_outdir -rf && \\
        spades.py --meta -t {params.spades.num_threads} -k {params.spades.kmers} \\
            ${input_R1} \\
            ${input_R2} \\
            -o spades_outdir -m ${params.spades.memory} && \\
        filterFastaByLength.pl \\
            --infile spades_outdir/scaffolds.fasta \\
            --length ${params.spades.length} > \\
            spades/scaffolds_gt${params.spades.length}.fasta && \\
        compileAssemblyResultsSingle.pl --infile spades/scaffolds_gt${params.spades.length}.fasta > spades/assembly_stats.txt
        """
}

process PRODIGAL {
    publishDir "$params.DEFAULT.outdir/gene_prediction/", mode: 'symlink'
    debug true
    
    input:
        path(infile)

    output:
        path("contigs_renamed.fna"), emit: predicted_genes_fna
        path("contigs_renamed.faa"), emit: predicted_genes_faa
        path("contigs_renamed.gff"), emit: predicted_genes_gff

    script:
        """
        module load ${params.modules.prodigal} ${params.modules.tools} && \\
        prodigal -i ${infile} -f gff -p meta \\
            -o contigs.gff \\
            -a contigs.faa \\
            -d contigs.fna && \\
        convertProdigalNames.pl \\
            --gff contigs.gff \\
            --fna contigs.fna \\
            --faa contigs.faa \\
            --renamed_gff contigs_renamed.gff \\
            --renamed_faa contigs_renamed.faa \\
            > contigs_renamed.fna
        """
}

process EXONERATE {
    debug true
    publishDir "$params.DEFAULT.outdir/assembly/", mode: 'symlink'
    
    input:
        val(sequence_type)
        path(infile)

    output:
        path("estimated_number_of_chunks_${sequence_type}.txt"), emit: chunks
        path("${sequence_type}_chunks/exonerate.done"), emit: done

    script:
        def to_delete = ""
        println ${infile.baseName}
        
        """
        module load ${params.modules.exonerate} ${params.modules.tools} && \\
        echo estimated_number_of_chunks_${sequence_type}.txt && \\
        echo estimated_number_of_chunks_${sequence_type}.txt
        """
}

/*
        if(${infile} =~ /.faa/){
            to_delete = ${sequence_type} + "/*.faa"
        }else if( ${infile} =~ /.fna/){
            to_delete = ${sequence_type} + "/*.fna"
        }else if(file_extension =~ /.fasta/){
            to_delete = ${sequence_type} + "/*.fasta"
        }
    if(${sequence_type} == "contigs"){
        publishDir "$params.DEFAULT.outdir/assembly/", mode: 'symlink'
    }else if(${sequence_type} == "genes"){
        publishDir "$params.DEFAULT.outdir/gene_prediction/", mode: 'symlink'
    }else{
        exit 1, "Invalid sequence type for EXONERATE"
    }
        estimateChunkFileSize.pl \\
            --infile {infile} \\
            --targeted_chunk_size ${params.exonerate.targeted_chunk_file_size_${sequence_type}} > \\
            estimated_number_of_chunks_${sequence_type}.txt && \\
        rm -rf {to_delete} && mkdir -p ${sequence_type} && \\
        fastasplit -f {infile} \\
          -o ${sequence_type}_chunks/ \\
          -c \`cat estimated_number_of_chunks_${sequence_type}.txt\` &&
        doublecheckExonerate.pl --indir {outdir} --prefix ${infile.baseName} --no_chunks \`cat estimated_number_of_chunks_${sequence_type}.txt\` && \\
        touch exonerate.done && \\
        rm -rf {to_delete}


*/
