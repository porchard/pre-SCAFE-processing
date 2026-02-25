#!/usr/bin/env nextflow

nextflow.enable.dsl=2

BAM_GLOB = params.bam_glob

process trim {

    publishDir "${params.results}/trimmed"
    container 'docker://porchard/general:20220406125608'
    memory '5 GB'
    time '10h'

    input:
    tuple val(library), path("in.bam")

    output:
    tuple val(library), path("${library}.trimmed.bam"), path("${library}.trimmed.bam.bai")

    """
    trim-for-scafe.py --input-bam in.bam --output-bam ${library}.trimmed.unsorted.bam --trim-to 100 --ggg-mismatches-allowed 1 --max-softclipping 4
    samtools sort -m 3G -o ${library}.trimmed.bam ${library}.trimmed.unsorted.bam
    samtools index ${library}.trimmed.bam
    rm ${library}.trimmed.unsorted.bam
    """

}


process filter_to_most_supported_5prime_ends {

    publishDir "${params.results}/filtered"
    container 'docker://porchard/general:20220406125608'
    memory '75 GB'
    time '10h'
    //label 'largemem'

    input:
    tuple val(library), path(bam), path(bam_index)

    output:
    tuple val(library), path("${library}.filtered.bam"), path("${library}.filtered.bam.bai")
    path("*.png")

    """
    filter-bam-to-most-supported-5prime-ends.py --bam-in $bam --bam-out ${library}.filtered.bam --prefix ${library}.
    samtools index ${library}.filtered.bam
    """

}


workflow {

    bams = Channel.fromPath(BAM_GLOB).map({it -> [it.getName().tokenize('.')[0], it]})

    trim(bams) | filter_to_most_supported_5prime_ends

}
