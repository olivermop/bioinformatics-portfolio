#bioinformatics-portfolio/scripts/02_demultiplexar.py
# This script demultiplexes the imported sequences, separating them by sample
# using the metadata file and barcodes.
# Output: Sample-separated reads (.qza)
import subprocess

comando = [
    "qiime", "demux", "emp-paired",
    "--m-barcodes-file", "../metadata/sample-metadata.tsv",
    "--m-barcodes-column", "barcode-sequence",
    "--p-rev-comp-mapping-barcodes",
    "--i-seqs", "../results/emp-paired-end-sequences.qza",
    "--o-per-sample-sequences", "../results/demux-full.qza",
    "--o-error-correction-details", "../results/demux-details.qza"
]

subprocess.run(comando, check=True)
print("Demultiplexing completed successfully.")
