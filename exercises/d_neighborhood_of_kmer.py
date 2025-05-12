# /bioinformatics-portfolio/exercises/d_neighborhood_of_kmer.py
# Python code to find the d-neighborhood of a k-mer (in this case TGCAT)

# Function to calculate the Hamming distance between two sequences
def hamming_distance(seq1, seq2):
    # Check that the sequences have the same length
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length.")

    # Initialize counter for the Hamming distance
    distance = 0

    # Compare each position of the sequences
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            distance += 1  # Increment the counter if the positions do not match

    return distance


# Function to generate all possible k-mers of a given length
def generate_kmers(k):
    from itertools import product
    # Generate all possible combinations of A, C, G, T
    return [''.join(p) for p in product('ACGT', repeat=k)]


# Function to find the d-neighborhood of a k-mer
def d_neighborhood(pattern, d):
    k = len(pattern)
    all_kmers = generate_kmers(k)  # Generate all k-mers of length k
    neighborhood = []

    # Add all k-mers that have a Hamming distance <= d
    for kmer in all_kmers:
        if hamming_distance(pattern, kmer) <= d:
            neighborhood.append(kmer)

    return neighborhood


# Example pattern and d value
pattern = "CCAGTCAATG"
d = 1

# Find the d-neighborhood of the pattern
neighborhood = d_neighborhood(pattern, d)

# Print the neighborhood and its size (number of k-mers)
print(f"The d-neighborhood of {pattern} is: {neighborhood}")
print(f"Number of k-mers in the neighborhood: {len(neighborhood)}")
