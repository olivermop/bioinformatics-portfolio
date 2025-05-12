# /bioinformatics-portfolio/exercises/pattern_occurrence_counter.py
# This script counts the occurrences of a specific pattern within a given DNA sequence.

def count_pattern_occurrences(sequence, pattern):
    # Initialize a counter for occurrences
    count = 0
    pattern_length = len(pattern)

    # Traverse the sequence to search for the pattern
    for i in range(len(sequence) - pattern_length + 1):
        if sequence[i:i + pattern_length] == pattern:
            count += 1  # Increment the counter when the pattern is found

    return count

# Example sequence and pattern from the exercise
sequence = "CGTGACAGTGTATGGGCATCTTT"
pattern = "TGT"

# Compute the number of pattern occurrences in the sequence
count = count_pattern_occurrences(sequence, pattern)
print(count)
