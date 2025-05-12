# /bioinformatics-portfolio/exercises/clump_finding2.py
# This code reads the genome from the file and searches for all 9-mers that
# appear at least 3 times within any 500-nucleotide window.
# The result is the number of distinct 9-mers that meet the criterion.

def read_genome(file_path):
    # Read the file content and remove newline characters
    with open(file_path, 'r') as file:
        return file.read().replace('\n', '')


def find_clumps(genome, k, L, t):
    # Create a set to store k-mers that form clumps
    patterns = set()
    n = len(genome)

    # Iterate over the genome to analyze each window of length L
    for i in range(n - L + 1):
        window = genome[i:i + L]
        freq_map = {}

        # Build a frequency map for k-mers in this window
        for j in range(L - k + 1):
            kmer = window[j:j + k]
            if kmer in freq_map:
                freq_map[kmer] += 1
            else:
                freq_map[kmer] = 1

        # Check if any k-mer appears at least t times
        for kmer, count in freq_map.items():
            if count >= t:
                patterns.add(kmer)

    return patterns


# Path to the E. coli genome file
genome_file_path = '/Users/Olivermop/Documents/bioinformatics_portfolio/data/Genomes/e_coli.txt'
genome = read_genome(genome_file_path)

# Parameters: k=9, L=500, t=3
k = 9
L = 500
t = 3

# Find clumps and print the number of distinct 9-mers forming clumps
clumps = find_clumps(genome, k, L, t)
print(f"Number of distinct 9-mers forming (500,3)-clumps: {len(clumps)}")
