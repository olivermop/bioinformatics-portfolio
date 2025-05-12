#!/usr/bin/env nextflow

process ANALYSIS {
    input:
    path input_file

    output:
    path "results/*.txt"

    script:
    """
    python3 /Users/Olivermop/Documents/bioinformatics_portfolio/exercises/rna_seq_analysis.py GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt
    """
}

workflow {
    params.input_file = file('data/GEO/GSE176078_Wu_etal_2021_bulkRNAseq_raw_counts.txt')

    ANALYSIS(params.input_file)
}
