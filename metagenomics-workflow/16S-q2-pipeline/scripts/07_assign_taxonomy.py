#bioinformatics-portfolio/07_assign_taxonomy.py
# This script performs taxonomic classification on the representative sequences
# using a pre-trained Naive Bayes classifier (e.g., SILVA or Greengenes).
# The classifier assigns taxonomy from domain to species level, depending on confidence.
# Output: taxonomy.qza — a QIIME 2 artifact containing taxonomic annotations.

import subprocess

comando = [
    "qiime", "feature-classifier", "classify-sklearn",
    "--i-classifier", "../data/gg-classifier.qza",
    "--i-reads", "../results/rep-seqs.qza",
    "--o-classification", "../results/taxonomy.qza"
]

subprocess.run(comando, check=True)
print("Taxonomy classification completed.")
