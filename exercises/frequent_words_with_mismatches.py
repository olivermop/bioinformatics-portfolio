# /bioinformatics-portfolio/exercises/frequent_words_with_mismatches.py
# This script finds the most frequent k-mers with up to d mismatches in a given DNA string.

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

# Function to find the most frequent k-mers with mismatches in a DNA string
def frequent_words_with_mismatches(Text, k, d):
    freq_map = {}  # Dictionary to store frequency of each k-mer
    n = len(Text)

    # Iterate over all k-mers in the text
    for i in range(n - k + 1):
        pattern = Text[i:i + k]
        neighborhood = neighbors(pattern, d)  # Get the d-neighborhood of the pattern

        # Count each k-mer in the neighborhood
        for neighbor in neighborhood:
            if neighbor not in freq_map:
                freq_map[neighbor] = 1
            else:
                freq_map[neighbor] += 1

    # Determine the maximum frequency
    max_count = max(freq_map.values())

    # Collect all k-mers with the maximum frequency
    most_frequent_kmers = [pattern for pattern, count in freq_map.items() if count == max_count]

    return most_frequent_kmers

# Example input
Text = "GATGAGCTACTACGATGTGTTACAGCTGTTCCTCCGATGTCCGATGTGTTCCAGCAGCTACAGCAGCTGTTGTTGTTCCTGTTGTTACTCCGATGGATGTACTCCTACTACTCCAGCTGTTGTGATGTCCGATGTACGATGTCCTGTGATGTGTTACTGTTCCTACTCCAGCTACTCCAGCGATGTGTAGCAGCGATGTACAGCAGCTGTTGTGATGGATGTCCTACGATGTGTTGTTACGATGTACAGCTACTCCTCCAGCTGTGATGTGTGATGAGCAGCGATGTCCTGTTACTACTGTTGTAGCTCCTCCTCCTCCTCCAGCTCCAGCTACTGTTGTAGCTCCAGCTACTCCTAC"
k = 7
d = 2

# Compute the most frequent k-mers with mismatches
result = frequent_words_with_mismatches(Text, k, d)

# Print the result
print(" ".join(result))
