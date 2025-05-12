#bioinformatics-portfolio/scripts/03_subsampling.py
# This script takes a fraction (30%) of the demultiplexed sequences.
# Useful for quick tests or when you want to reduce computation time.
import subprocess

comando = [
    "qiime", "demux", "subsample-paired",
    "--i-sequences", "../results/demux-full.qza",
    "--p-fraction", "0.3",  # value 30%
    "--o-subsampled-sequences", "../results/demux-subsample.qza"
]

subprocess.run(comando, check=True)
print("Subsampling completed.")
