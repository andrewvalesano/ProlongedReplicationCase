# Analysis of sequences from patient with prolonged SARS-CoV-2 infection

This repository holds code and data for the analysis of longitudinal specimens from a patient with prolonged infection with SARS-CoV-2. The code and data is provided here for reproducibility. Below is a map of the repository.

# Overview
--------

    project
    |- README          # The top level description of content. You are here.
    |
    |- data  
    |  |- metadata/   # Sample metadata for the specimens sequenced in this study.
    |  |- reference/  # Reference fasta files used in pipelines.
    |  |- raw/        # Raw data.
    |- scripts/       # Code for analyses presented in manuscript.
    |- pipelines/     # Snakemake files for alignment, consensus calling, and variant calling.
    
  --------

# Notes

Raw sequence reads are available at the Sequence Read Archive, BioProject accession PRJNA662589.

# Contact

If you have questions, please contact the [Lauring Lab](https://lauringlab.wordpress.com/contacts/).
