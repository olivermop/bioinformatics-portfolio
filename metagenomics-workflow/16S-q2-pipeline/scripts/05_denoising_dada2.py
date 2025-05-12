#bioinformatics-portfolio/scripts/05_denoising_dada2.py
# This script performs cleanup and error correction on reads
# using QIIME 2's DADA2 algorithm for paired-end sequences.
# It uses the trimming parameters determined in the quality visualization.
# Outputs:
# - Feature table (table.qza)
# - Representative sequences (rep-seqs.qza)
# - Denoising statistics (denoising-stats.qza)

import subprocess

comando = [
    "qiime", "dada2", "denoise-paired",
    "--i-demultiplexed-seqs", "../results/demux-subsample.qza",  # Note: We will use filtered demux in step 6
    "--p-trim-left-f", "13",
    "--p-trim-left-r", "13",
    "--p-trunc-len-f", "150",
    "--p-trunc-len-r", "150",
    "--o-table", "../results/table.qza",
    "--o-representative-sequences", "../results/rep-seqs.qza",
    "--o-denoising-stats", "../results/denoising-stats.qza"
]

subprocess.run(comando, check=True)
print("Denoising completed with DADA2.")
