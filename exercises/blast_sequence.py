# /bioinformatics-portfolio/exercises/blast_sequence.py
# This script takes the translated protein sequences of BRCA1 and BRCA2 and performs a BLAST search on NCBI.
# It splits the sequences into smaller fragments (<10,000 residues) to avoid CPU limit errors
# when performing multiple BLAST searches. The results are stored in XML format for each fragment.
# XML is used because it is the standard format for NCBI BLAST search results.

from Bio.Blast import NCBIWWW
import os

# Path to translated protein files
brca1_protein_file = "results/protein_translation_brca1.txt"
brca2_protein_file = "results/protein_translation_brca2.txt"

# Function to read protein sequences
def read_protein_sequence(file_path):
    # Read the protein sequence from a FASTA-like file (ignoring the header)
    with open(file_path) as f:
        return f.read().strip().split('\n')[1]

# Read the protein sequences
brca1_sequence = read_protein_sequence(brca1_protein_file)
brca2_sequence = read_protein_sequence(brca2_protein_file)

# Function to split sequences into smaller fragments
def split_sequence(sequence, max_length=10000):
    # Divide the sequence into chunks of up to max_length
    return [sequence[i:i + max_length] for i in range(0, len(sequence), max_length)]

# Split sequences into fragments
brca1_fragments = split_sequence(brca1_sequence)
brca2_fragments = split_sequence(brca2_sequence)

# Results directory
results_dir = "results"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Function to perform BLAST search
def perform_blast(sequence_fragments, gene_name):
    # Iterate over each fragment and submit a BLAST search
    for i, fragment in enumerate(sequence_fragments, 1):
        print(f"Starting BLAST search on NCBI for {gene_name}, fragment {i}...")
        try:
            blast_handle = NCBIWWW.qblast("blastp", "nr", fragment)
            output_file = os.path.join(results_dir, f"blast_result_{gene_name}_fragment_{i}.xml")
            with open(output_file, "w") as result_file:
                result_file.write(blast_handle.read())
            print(f"Completed BLAST search for {gene_name}, fragment {i}. Results saved in {output_file}")
        except Exception as e:
            print(f"Error during BLAST search for {gene_name}, fragment {i}: {e}")

# Perform BLAST searches for BRCA1 and BRCA2
perform_blast(brca1_fragments, "BRCA1")
perform_blast(brca2_fragments, "BRCA2")
