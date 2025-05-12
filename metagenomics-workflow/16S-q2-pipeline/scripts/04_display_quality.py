#bioinformatics-portfolio/scripts/04_display_quality.py
# This script generates a visual summary (.qzv) of the quality of the subsampled sequences.
# The result allows you to view the quality of forward and reverse reads by position.
# This is essential for deciding the trimming parameters (trim/trunc) in the denoising step with DADA2.

import subprocess

comando = [
    "qiime", "demux", "summarize",
    "--i-data", "../results/demux-subsample.qza",
    "--o-visualization", "../results/demux-subsample.qzv"
]

subprocess.run(comando, check=True)
print("Quality visualization generated.")
