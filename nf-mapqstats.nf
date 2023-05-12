params.input_bam = "$projectDir/test/data/sample{A,B}.bam"
params.genome_file = "$projectDir/test/data/genome.fa"
params.output_bed_files = "$projectDir/test/data/{sample}_mapqstats.bed"
params.outdir = "$projectDir/test/results/"
params.window_size = 500

sampleNameSeparator = "."

bam_files_ch = Channel.fromPath(params.input_bam, checkIfExists: true).map {
    tuple( it.name.split("\\.")[0], it)
    }
.subscribe onNext: { println it }, onComplete: { println 'Done' }


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
    input:
    path genome_fasta
    

    output: 
    path 'genomesize_file'

    script:
    """
    samtools faidx --fai-idx /dev/stdout ${genome_fasta} | cut -f1-2 > genomesize_file
    """

}

process createGenomicWindows {
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
    publishDir params.outdir, mode: 'copy'

    input: 
    tuple val(sample_id), path(bamfile)
    path window_file

    output:
    path "mapq_${sample_id}.bed.gz"
    
    script: 
    """
    bedtools intersect -a ${window_file} -wa -wb -bed -b ${bamfile} -wa | \
    bedtools merge -d -1 -c 8 -o collapse,count | gzip -c > mapq_${sample_id}.bed.gz"
    """
}

workflow {
    genomesize_ch = createGenomeSizeFile(params.genome_file)
    genomic_windows_ch = createGenomicWindows(genomesize_ch, params.window_size)
    bam_files_ch.view()
    mapq_in_windows_ch = getMAPQinWindows(bam_files_ch, genomic_windows_ch)
    genomic_windows_ch.view() 
}


