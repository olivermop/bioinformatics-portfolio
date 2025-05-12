#bioinformatics-portfolio/08_visualize_taxonomy_barplot.py
# This script generates an interactive taxonomy bar plot visualization.
# It maps the assigned taxonomy to the feature table, allowing you to see
# the relative abundance of taxonomic groups across samples.
# Output: taxonomy-barplot.qzv — view at https://view.qiime2.org


import subprocess

comando = [
    "qiime", "taxa", "barplot",
    "--i-table", "../results/table.qza",
    "--i-taxonomy", "../results/taxonomy.qza",
    "--m-metadata-file", "../metadata/sample-metadata.tsv",
    "--o-visualization", "../results/taxonomy-barplot.qzv"
]

subprocess.run(comando, check=True)
print("✅ Taxonomy bar plot created successfully.")
