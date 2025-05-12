# /bioinformatics-portfolio/exercises/questionnaire1.py
# Code to solve five DNA sequence–based questions

# Function 1: Count how many times a pattern appears in a sequence
def count_pattern_occurrences(sequence, pattern):
    count = 0
    pattern_length = len(pattern)
    for i in range(len(sequence) - pattern_length + 1):
        if sequence[i:i + pattern_length] == pattern:
            count += 1
    return count

# Function 2: Find the most frequent k-mer in a given sequence
def most_frequent_kmer(sequence, k):
    kmer_frequencies = frequency_table(sequence, k)
    most_frequent = max(kmer_frequencies, key=kmer_frequencies.get)
    return most_frequent, kmer_frequencies[most_frequent]

# Function 3: Get the reverse complement of a DNA sequence
def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement[base] for base in reversed(dna_sequence))

# Function 4: Find all starting positions of a pattern in a sequence (0-based indexing)
def pattern_matching_positions(text, pattern):
    positions = []
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i + len(pattern)] == pattern:
            positions.append(i)
    return positions

# Function 5: Build a frequency table of k-mers in a text window
def frequency_table(text, k):
    frequency_map = {}
    for i in range(len(text) - k + 1):
        k_mer = text[i:i + k]
        frequency_map[k_mer] = frequency_map.get(k_mer, 0) + 1
    return frequency_map

# Question 1: Count occurrences of "CGCG" in a specific sequence
sequence1 = "CGCGATACGTTACATACATGATAGACCGCGCGATCATATCGCGATTATC"
pattern1 = "CGCG"
count1 = count_pattern_occurrences(sequence1, pattern1)
print(f"Question 1 – Count of '{pattern1}': {count1}")

# Question 2: Find the most frequent 3-mer in a sequence
sequence2 = "TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT"
most_frequent_3mer, freq = most_frequent_kmer(sequence2, 3)
print(f"Question 2 – Most frequent 3-mer is '{most_frequent_3mer}' with {freq} occurrences.")

# Question 3: Get the reverse complement of a given sequence
sequence3 = "GCTAGCT"
reverse_complement3 = reverse_complement(sequence3)
print(f"Question 3 – Reverse complement of '{sequence3}' is '{reverse_complement3}'.")

# Question 4: Find all starting positions of "CGC" in a sequence
sequence4 = "ATGACTTCGCTGTTACGCGC"
pattern4 = "CGC"
positions4 = pattern_matching_positions(sequence4, pattern4)
print(f"Question 4 – Positions of '{pattern4}': {' '.join(map(str, positions4))}")

# Question 5: Find all starting positions of "CGC" in a sequence
# (Note: same pattern as Question 4—replace pattern if intended)
sequence5 = "ATGACTTCGCTGTTACGCGC"
pattern5 = "CGC"
positions5 = pattern_matching_positions(sequence5, pattern5)
print(f"Question 5 – Positions of '{pattern5}': {' '.join(map(str, positions5))}")
