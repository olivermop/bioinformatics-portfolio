# /bioinformatics-portfolio/exercises/clump_finding
# This code implements the "clump" search algorithm,
# which is useful in detecting regions of the genome that have a high concentration of certain patterns (k-mers).
# The idea behind the origin of replication (ori) is that certain sequences, such as DNAA boxes,
# appear multiple times in a relatively short region of the genome.

def frequency_table(text, k):
    frequency_map = {}
    # Iterate through the text and generate k-mers
    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        # If the pattern is already in the frequency map, increment its count; otherwise initialize it to 1
        if pattern in frequency_map:
            frequency_map[pattern] += 1
        else:
            frequency_map[pattern] = 1
    return frequency_map

# Function to find patterns that form "clumps" in a genome
def find_clumps(genome, k, L, t):
    patterns = set()  # Use a set to store distinct k-mers that form clumps
    n = len(genome)   # Length of the genome

    # Slide a window of length L across the genome
    for i in range(n - L + 1):
        window = genome[i:i + L]
        # Build a frequency table of k-mers in the current window
        frequency_map = frequency_table(window, k)

        # Check if any k-mer occurs at least t times in the current window
        for kmer, count in frequency_map.items():
            if count >= t:
                patterns.add(kmer)  # Add the k-mer to the set (duplicates are ignored)

    return patterns

# Example dataset
genome = ""  # Insert the genome sequence here
k = 10       # Length of the k-mer
L = 100      # Window size
t = 4        # Threshold for k-mer repetitions within the window

# Find the clumps in the genome
clumps = find_clumps(genome, k, L, t)

# Print the result as a space-separated list of k-mers
print(" ".join(clumps))

