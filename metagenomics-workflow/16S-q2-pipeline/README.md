# 16S Paired-End Microbiome Analysis with QIIME 2 (via Docker)

This project contains a complete and reproducible pipeline for 16S rRNA paired-end sequence analysis using QIIME 2 (2024.2) in Docker.

The goal is to convert raw sequencing data into meaningful biological insights: denoising, taxonomic classification, visualizations, and exportable results.

---

## Project Structure

16S-q2-pipeline/
├── data/ # Input files (FASTQ and classifier)
│ └── emp-paired-end-sequences/
│ ├── forward.fastq.gz
│ ├── reverse.fastq.gz
│ └── barcodes.fastq.gz
│ └── silva-138-99-nb-classifier.qza
├── metadata/
│ └── sample-metadata.tsv
├── scripts/ # Python scripts for each QIIME 2 step
├── results/ # Output artifacts (.qza, .qzv)
├── exports/ # Final exportable data (TSV, FASTA)
└── README.md


---

##  Requirements

- [Docker Desktop](https://www.docker.com/products/docker-desktop/)
- QIIME 2 Docker image (amplicon flavor)

---

## Running the pipeline with Docker


