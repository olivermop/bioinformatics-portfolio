# bioinformatics-portfolio/exercises/approximate_pattern_count.py
# This script implements the ApproximatePatternCount algorithm.
# Counts how many times a pattern appears in the text with up to d differences (mismatches).

# Function to calculate the Hamming distance between two sequences
def hamming_distance(p, q):
    # Count the differences between the two sequences
    return sum(1 for i in range(len(p)) if p[i] != q[i])


# Function to count approximate pattern occurrences in the text with up to d mismatches
def approximate_pattern_count(Text, Pattern, d):
    count = 0  # Initialize the occurrence counter
    pattern_length = len(Pattern)

    # Iterate over the text and compare each substring of pattern length
    for i in range(len(Text) - pattern_length + 1):
        substring = Text[i:i + pattern_length]
        # If the Hamming distance is less than or equal to d, increment the counter
        if hamming_distance(Pattern, substring) <= d:
            count += 1

    return count  # Return the total number of occurrences


# Example input
Pattern = "ACCGC"
Text = "CATTCACGCCGCAGCCTAACAATGAATTAGTACGACTCTGAGGTCATGGATTAGGGATTTCGTTACAAAATGAACGGTCGGTTCAGAATGGACCTGGTGCAGTCACCTGTTTCGCCAGCCGGAGGTCATGACAGTTGGACAATCACGGGTTGGGCGCGCTGATCGGTCTTAAATGCTCGCGGGGGCTCAGGACCTTGGAGGTGAGTTCTGGGGTATTACACGATCCATAAATACTCATAACCGATTAGCTTTGATTCGATTAATGCCAAACCCTTAATTTCACTGCTGGACTTACCCCTAATTACCGCGCCGGTTTGTATCTCACCTGATCGTAATACTTTTCTATCGGCTTGAACCATGACC"
d = 2

# Compute the result
result = approximate_pattern_count(Text, Pattern, d)

print(result)
