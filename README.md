# nf-mapstats

A nextflow pipeline to extract MAPQ across genomic windows from a BAM file using samtools and bedtools. 

## How to run

```
nextflow run nf-mapqstats.nf \
	 --genome_file <FASTA_FILE> \ 
	 --input_bam <BAM_FILE> 
