# /bioinformatics-portfolio/exercises/sequence_analysis.py
# This script analyzes the frequency of nucleotide bases (A, T, C, G) in the BRCA1 and BRCA2 DNA sequences.
# It uses Biopython’s SeqIO module to read sequences from FASTA files and computes the frequency of each base.

from Bio import SeqIO
import os

# Paths to the FASTA files for BRCA1 and BRCA2
brca1_file = "/Users/Olivermop/Documents/bioinformatics_portfolio/data/BRCA1_datasets/ncbi_dataset/data/gene.fna"
brca2_file = "/Users/Olivermop/Documents/bioinformatics_portfolio/data/BRCA2_datasets/ncbi_dataset/data/gene.fna"

# Create the "results" directory if it does not exist
results_dir = "results"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Create the output file for base frequency results
output_file_path = os.path.join(results_dir, "base_frequency_brca1_brca2.txt")
with open(output_file_path, "w") as output_file:
    output_file.write("Base Frequency Analysis for BRCA1 and BRCA2 Genes:\n")
    output_file.write("This analysis calculates the frequency of nucleotide bases (A, T, C, G) in the DNA sequences.\n\n")

    # Function to calculate and write base frequencies for a given FASTA file
    def calculate_base_frequency(fasta_file, gene_name):
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = record.seq
            base_counts = {base: sequence.count(base) for base in ['A', 'T', 'C', 'G']}
            total_bases = len(sequence)

            # Write results to the output file
            output_file.write(f"Results for {gene_name}:\n")
            output_file.write(f"Total bases: {total_bases}\n")
            for base, count in base_counts.items():
                frequency = (count / total_bases) * 100
                output_file.write(f"{base}: {count} ({frequency:.2f}%)\n")
            output_file.write("\n")

    # Calculate and write base frequencies for BRCA1
    calculate_base_frequency(brca1_file, "BRCA1")

    # Calculate and write base frequencies for BRCA2
    calculate_base_frequency(brca2_file, "BRCA2")

print(f"Analysis complete. Results saved to {output_file_path}")
