params.output_prefix = ""

process getMAPQforRegions {
    // set system requirements
    memory { 2.GB * task.attempt }
    cpus = 2
    time = 8.h
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
    maxRetries = 8

    // define dependencies for conda
    conda (params.enable_conda ? "bioconda::bedtools=2.25" : null)


    tag "get MAPQs for ${bamfile} (Sample: ${sample_id})"
    publishDir params.outdir, mode: 'copy'

    input: 
    path window_file
    path genome_file
    tuple val(sample_id), path(bamfile)
    

    output:
    path "${params.output_prefix}_${sample_id}.bed.gz"
    
    script: 
    """
    set -euf -o pipefail
    bedtools intersect -sorted -g ${genome_file} -a ${window_file} -wa -wb -bed -b ${bamfile} -wa | bedtools merge -d -1 -c 8 -o collapse,count | gzip -c > ${params.output_prefix}_${sample_id}.bed.gz
    """
}
