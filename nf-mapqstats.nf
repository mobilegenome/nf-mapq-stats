params.input_bam = "$projectDir/test/data/sample{A,B}.bam"
params.input_vcf = "$projectDir/test/data/*.vcf.gz"
params.genome_file = "$projectDir/test/data/genome.fa"
params.outdir = "$projectDir/test/results/"
params.window_size = 500
params.sv_qual_min_threshold = 30
params.sv_len_max_threshold = 20000
params.mapqinwindows = false

sampleNameSeparator = "."


// get channel of BAM files
// in the format of tuples consisting of (filename, bampath)
def bam_files_ch = Channel.fromPath(params.input_bam, checkIfExists: true)
    .map {
    tuple( it.name.split("\\.")[0], it)
    }

bam_files_ch.view()

// get channel of VCF files
// in the format of tuples consisting of (filename, bampath)
def vcf_files_ch = Channel.fromPath(params.input_vcf, checkIfExists: true)
    .map {
    tuple( it.name.split("\\.")[0], it)
    }

vcf_files_ch.view()


log.info """\
    MAPQ STATISTICS
    ===============
    bams: ${params.input_bam}
    vcfs: ${params.input_vcf}
    genome_fasta: ${params.genome_file}
    outdir: ${params.outdir}
    window_size: ${params.window_size}
    """
    .stripIndent(true)

process createGenomeSizeFile {
    memory = 2.GB
    cpus = 2
    time = 1.h

    conda (params.enable_conda ? "bioconda::samtools=1.17" : null)


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
    conda (params.enable_conda ? "bioconda::samtools=1.17 bioconda::bedtools=2.25" : null)


    publishDir params.outdir, mode: 'copy'

    input:
    path genomesize_file
    val window_size
    

    output: 
    path "genomic_windows.${window_size}.bed"

    script:
    """
    bedtools makewindows -g ${genomesize_file} -w ${window_size} > 'genomic_windows.${window_size}.bed'
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch 'genomic_windows.${window_size}.bed'
    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}

process convertVCFtoBED {
    // set system requirements
    memory { 2.GB * task.attempt }
    cpus = 2
    time = 8.h
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
    maxRetries = 8

    // define dependencies for conda
    conda (params.enable_conda ? "bioconda::cyvcf python=3.6.8" : null)
    
    tag "convert VCF to BED for  ${vcf_file}"
    publishDir params.outdir, mode: 'copy'


    input: 
    tuple val(sample_id), path(vcf_file)
    val qual_min_threshold
    val len_max_threshold

    output:
    tuple val(sample_id), path("sv_${sample_id}.bed")

    script:
    """
    #!/usr/bin/env python
    from cyvcf2 import VCF

    output_fh = open("sv_${sample_id}.bed", "w")

    for variant in VCF("${vcf_file}"): # or VCF('some.bcf')
        if (variant.INFO.get("SVTYPE") == "DEL") and \
            (variant.QUAL >= ${qual_min_threshold}):

            chrom = variant.CHROM
            start = variant.POS + variant.INFO.get("CIPOS")[1]  # get position from inner confidence interval
            qual = round(variant.QUAL,3)
            end = variant.INFO.get("END") - variant.INFO.get("CIEND")[0]  # get position from inner confidence interval
            strand = "." # igore strand
            id = variant.ID

            # check end larger start
            assert end > start
            
            if abs(variant.INFO.get("SVLEN")) >= ${len_max_threshold}:
                continue

            bed_line = [
                chrom, 
                start,
                end,
                id,
                qual,
                strand
            ]
            

            print("\t".join([str(elem) for elem in bed_line]), file=output_fh)
    
    output_fh.close()
    
    """
}

process SortBed {
    // set system requirements
    memory { 2.GB * task.attempt }
    cpus = 2
    time = 8.h
    errorStrategy { task.exitStatus == 130 ? 'retry' : 'terminate' }
    maxRetries = 8

    conda (params.enable_conda ? "bioconda::bedtools=2.25" : null)

    input:
    tuple val(sample_id), path(bed_unsorted)
    path genomesize_file

    output:
    path("sv_${sample_id}.sorted.bed")
    //tuple val(sample_id), path("sv_${sample_id}.sorted.bed")


	script:
    """
    # use `-faidx genimesize_file to maintain correct chromosome/contig order 
	bedtools sort -faidx ${genomesize_file} -i ${bed_unsorted} | cut -f1-3 >  sv_${sample_id}.sorted.bed

    """
}

include { getMAPQforRegions as getMAPQforGenotypes} from './modules/getMAPQforRegions.nf' addParams(output_prefix: 'svMAPQ')
include { getMAPQforRegions as getMAPQforWindows} from './modules/getMAPQforRegions.nf'  addParams(output_prefix: 'windowMAPQ')

workflow {


    genomesize_ch = createGenomeSizeFile(params.genome_file)
    if ( params.mapqinwindows ) {
        genomic_windows_ch = createGenomicWindows(genomesize_ch, params.window_size)
        getMAPQforWindows(genomic_windows_ch, genomesize_ch, bam_files_ch)
    }
    
    svGenotypeBED_ch = convertVCFtoBED(vcf_files_ch, params.sv_qual_min_threshold, params.sv_len_max_threshold)
    svGenotypeBED_ch.view()
    svGenotypeSortedBed_ch = SortBed(svGenotypeBED_ch, genomesize_ch)
    svGenotypeSortedBed_ch.view()
    getMAPQforGenotypes(svGenotypeSortedBed_ch, genomesize_ch, bam_files_ch)

}


