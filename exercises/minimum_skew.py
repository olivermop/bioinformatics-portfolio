# /bioinformatics-portfolio/exercises/minimum_skew.py
# This script solves the Minimum Skew Problem for a DNA string read from a file.

# Function to calculate the skew of a DNA sequence.
# Skew increases by 1 for each 'G' and decreases by 1 for each 'C'.
def calculate_skew(sequence):
    skew = [0]  # Initialize skew with 0 at position 0
    for i in range(len(sequence)):
        if sequence[i] == 'G':
            skew.append(skew[-1] + 1)  # Increment skew for each 'G'
        elif sequence[i] == 'C':
            skew.append(skew[-1] - 1)  # Decrement skew for each 'C'
        else:
            skew.append(skew[-1])  # No change in skew for 'A' or 'T'
    return skew

# Function to find all positions where the skew reaches its minimum value.
def find_min_skew_positions(sequence):
    skew = calculate_skew(sequence)
    min_skew_value = min(skew)  # Find the minimum skew value
    # Collect positions where the skew equals the minimum value
    min_positions = [i for i, value in enumerate(skew) if value == min_skew_value]
    return min_positions

# Path to the input file containing the genome sequence
file_path = "/Users/Olivermop/Documents/bioinformatics_portfolio/data/Genomes/dataset_30277_10 (1).txt"

# Read the genome sequence from the file
with open(file_path, 'r') as file:
    genome = file.read().strip()  # Read content and remove whitespace

# Find positions where the skew is minimal
min_skew_positions = find_min_skew_positions(genome)

print(f"Minimum skew positions: {min_skew_positions}")
