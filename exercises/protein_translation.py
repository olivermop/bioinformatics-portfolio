# /bioinformatics-portfolio/exercises/protein_translation.py
# This script reads DNA sequences from the FASTA files for BRCA1 and BRCA2,
# translates them into protein sequences, and saves the results in FASTA format.
# It uses Biopython’s SeqIO module to parse the FASTA and the translate() method
# to convert nucleotide sequences to amino acids.

from Bio.Seq import Seq
from Bio import SeqIO
import os

# Paths to the FASTA files for BRCA1 and BRCA2 DNA sequences
brca1_file = "/Users/Olivermop/Documents/bioinformatics_portfolio/data/BRCA1_datasets/ncbi_dataset/data/gene.fna"
brca2_file = "/Users/Olivermop/Documents/bioinformatics_portfolio/data/BRCA2_datasets/ncbi_dataset/data/gene.fna"

# Create the "results" directory if it does not exist
results_dir = "results"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Output file paths for translated proteins
brca1_output_file_path = os.path.join(results_dir, "protein_translation_brca1.txt")
brca2_output_file_path = os.path.join(results_dir, "protein_translation_brca2.txt")

# Function to translate DNA to protein and write results in FASTA format
def translate_and_write(fasta_file, gene_name, output_path):
    with open(output_path, "w") as output_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = record.seq

            # Pad the sequence with 'N's if its length is not a multiple of three
            if len(sequence) % 3 != 0:
                sequence += 'N' * (3 - len(sequence) % 3)

            # Translate DNA sequence to protein sequence
            protein_sequence = sequence.translate()

            # Write the translated protein in FASTA format
            output_file.write(f">{gene_name}_translated_protein\n")
            output_file.write(f"{protein_sequence}\n\n")

# Translate and write results for BRCA1
translate_and_write(brca1_file, "BRCA1", brca1_output_file_path)

# Translate and write results for BRCA2
translate_and_write(brca2_file, "BRCA2", brca2_output_file_path)

print(f"Translation complete. Results saved to {brca1_output_file_path} and {brca2_output_file_path}")
