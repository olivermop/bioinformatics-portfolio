# /bioinformatics-portfolio/exercises/frequent_words_with_mismatches_and_reverse_complements.py
# This script finds the most frequent k-mers with up to d mismatches and their reverse complements in a given DNA string.

from itertools import product

# Function to generate all possible neighbors within d mismatches
def neighbors(pattern, d):
    nucleotides = ['A', 'C', 'G', 'T']
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return set(nucleotides)

    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)
    for text in suffix_neighbors:
        # If the suffix has fewer than d differences, substitute the first character
        if hamming_distance(pattern[1:], text) < d:
            for nucleotide in nucleotides:
                neighborhood.add(nucleotide + text)
        else:
            neighborhood.add(pattern[0] + text)

    return neighborhood

# Function to calculate the Hamming distance between two sequences
def hamming_distance(p, q):
    # Count the positions where the sequences differ
    return sum(1 for i in range(len(p)) if p[i] != q[i])

# Function to compute the reverse complement of a DNA string
def reverse_complement(pattern):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # Build the reverse complement sequence
    return ''.join(complement[base] for base in reversed(pattern))

# Function to find the most frequent k-mers with mismatches and reverse complements
def frequent_words_with_mismatches_and_reverse_complements(Text, k, d):
    freq_map = {}  # Dictionary to store frequency of each k-mer and its reverse complement
    n = len(Text)

    # Iterate over all k-mers in the text
    for i in range(n - k + 1):
        pattern = Text[i:i + k]
        rev_pattern = reverse_complement(pattern)
        neighborhood = neighbors(pattern, d) | neighbors(rev_pattern, d)

        # Count each k-mer and its reverse complement in the neighborhood
        for neighbor in neighborhood:
            freq_map[neighbor] = freq_map.get(neighbor, 0) + 1

    # Determine the maximum frequency
    max_count = max(freq_map.values())

    # Collect all k-mers with the maximum frequency
    most_frequent_kmers = [pattern for pattern, count in freq_map.items() if count == max_count]

    return most_frequent_kmers

# Example input
Text = "AACGTATACTTGTACAACGCGTATATATACCGTACTACTACCGAAAACGCGTACAATACGAACGCGCGTACGTACGAAAATACGTATACTACTTGTAAACGCGTTGTACTTGTACTACTACCGTATACTACGAATTGTACTACGTATTGTTGTATTGAATATTGTTGCGTTGCGAATTGCGTTGTTGTACTTGTAAATACTTGTACTACTACAATAAA"
k = 6
d = 3

# Compute the most frequent k-mers with mismatches and reverse complements
result = frequent_words_with_mismatches_and_reverse_complements(Text, k, d)

print(" ".join(result))
