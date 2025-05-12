#bioinformatics-portfolio/scripts/01_import_data.py
# This script imports paired .fastq.gz files (forward, reverse, and barcodes)
# as a QIIME 2 artifact (.qza) of type EMPPairedEndSequences.
# It is the first step in the analysis to prepare the raw data for demultiplexing.
import subprocess

comando = [
    "qiime", "tools", "import",
    "--type", "EMPPairedEndSequences",
    "--input-path", "../data/emp-paired-end-sequences",
    "--output-path", "../results/emp-paired-end-sequences.qza"
]

subprocess.run(comando, check=True)
print("Sequences imported successfully.")
