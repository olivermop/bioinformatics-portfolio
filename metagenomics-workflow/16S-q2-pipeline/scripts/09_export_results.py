#bioinformatics-portfolio/09_export_results.py
# This script exports QIIME 2 artifacts to plain-text formats for external analysis.
# Outputs include:
# - Feature table in BIOM format
# - Taxonomic assignments in TSV
# - Representative sequences in FASTA
# These files can be used in R, Excel, Python, or for further downstream tools.


import subprocess
import os

# Create export directory if it doesn't exist
os.makedirs("../exports", exist_ok=True)

comandos = [
    # Export feature table as BIOM
    [
        "qiime", "tools", "export",
        "--input-path", "../results/table.qza",
        "--output-path", "../exports/feature-table"
    ],
    # Export taxonomy
    [
        "qiime", "tools", "export",
        "--input-path", "../results/taxonomy.qza",
        "--output-path", "../exports/taxonomy"
    ],
    # Export representative sequences
    [
        "qiime", "tools", "export",
        "--input-path", "../results/rep-seqs.qza",
        "--output-path", "../exports/rep-seqs"
    ]
]

for cmd in comandos:
    subprocess.run(cmd, check=True)

print("✅ Export completed. Files available in the 'exports/' folder.")
