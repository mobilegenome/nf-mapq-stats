params.input_bam = "$projectDir/test/data/sample{A,B}.bam"
params.genome_file = "$projectDir/test/data/genome.fa"
params.outdir = "$projectDir/test/results/"
params.window_size = 500

sampleNameSeparator = "."


// get channel of BAM files
// in the format of tuples consisting of (filename, bampath)
bam_files_ch = Channel.fromPath(params.input_bam, checkIfExists: true).map {
    tuple( it.name.split("\\.")[0], it)
    }

bam_files_ch.view()


log.info """\
    MAPQ STATISTICS
    ===============
    bams: ${params.input_bam}
    genome_fasta: ${params.genome_file}
    outdir: ${params.outdir}
    window_size: ${params.window_size}
    """
    .stripIndent(true)

process createGenomeSizeFile {
    memory = 2.GB
    cpus = 2
    time = 1.h

    conda 'environment.yml'


    input:
    path genome_fasta
    

    output: 
    path 'genomesize_file'

    script:
    """
    samtools version && 
    samtools faidx --fai-idx /dev/stdout ${genome_fasta} | cut -f1-2 > genomesize_file
    """

}

process createGenomicWindows {
    memory = 8.GB
    cpus = 2
    time = 8.h
    conda 'environment.yml'


    publishDir params.outdir, mode: 'copy'

    input:
    path genomesize_file
    val window_size
    

    output: 
    path "genomic_windows.${window_size}.bed"

    script:
    """
    bedtools makewindows -g ${genomesize_file} -w ${window_size} > 'genomic_windows.${window_size}.bed'
    """

}


process getMAPQinWindows {
    memory = 8.GB
    cpus = 2
    time = 8.h
    
    conda 'environment.yml'


    tag "get MAPQs for ${bamfile} (Sample: ${sample_id})"
    publishDir params.outdir, mode: 'copy'

    input: 
    path window_file
    tuple val(sample_id), path(bamfile)

    output:
    path "mapq_${sample_id}.bed.gz"
    
    script: 
    """
    bedtools intersect -a ${window_file} -wa -wb -bed -b ${bamfile} -wa | bedtools merge -d -1 -c 8 -o collapse,count | gzip -c > mapq_${sample_id}.bed.gz
    """
}

workflow {
    genomesize_ch = createGenomeSizeFile(params.genome_file)
    genomic_windows_ch = createGenomicWindows(genomesize_ch, params.window_size)
    getMAPQinWindows(genomic_windows_ch, bam_files_ch)
}


