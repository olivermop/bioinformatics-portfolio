#bioinformatics-portfolio/06_visualize_denoising_stats.py
# This script generates visualizations for the output of the DADA2 denoising step.
# It includes:
# - A summary of the feature table showing the number of sequences per sample
# - A table of representative sequences
# - A summary of denoising statistics
# These visualizations can be viewed at https://view.qiime2.org


import subprocess

comandos = [
    [
        "qiime", "feature-table", "summarize",
        "--i-table", "../results/table.qza",
        "--o-visualization", "../results/table.qzv",
        "--m-sample-metadata-file", "../metadata/sample-metadata.tsv"
    ],
    [
        "qiime", "feature-table", "tabulate-seqs",
        "--i-data", "../results/rep-seqs.qza",
        "--o-visualization", "../results/rep-seqs.qzv"
    ],
    [
        "qiime", "metadata", "tabulate",
        "--m-input-file", "../results/denoising-stats.qza",
        "--o-visualization", "../results/denoising-stats.qzv"
    ]
]

for cmd in comandos:
    subprocess.run(cmd, check=True)

print("Views generated: tabla.qzv, secuencias.qzv, denoising-stats.qzv")
