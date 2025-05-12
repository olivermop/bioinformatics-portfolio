#/bioinformatics-portfolio/exercises/estimate_ktl.py
# Heuristic-based code to automatically obtain optimal values of k, L, and t by analyzing a genome.
# Use this before clump finding to detect regions with high concentrations of repeated patterns (clumps).
import collections


# Function to calculate k-mer diversity in the genome
def estimate_k(genome):
    k_values = range(4, 13)  # Range of k to test (from 4 to 12)
    diversity = {}

    for k in k_values:
        k_mer_set = set()
        for i in range(len(genome) - k + 1):
            k_mer_set.add(genome[i:i + k])
        diversity[k] = len(k_mer_set)  # Store the number of unique k-mers for each k value

    # Select the k that maximizes diversity without being too high
    optimal_k = max(diversity, key=diversity.get)
    return optimal_k, diversity


# Function to calculate the average distance between k-mer repeats and estimate the window size L
def estimate_L(genome, k):
    distances = []
    k_mer_positions = collections.defaultdict(list)

    # Find all positions of each k-mer
    for i in range(len(genome) - k + 1):
        k_mer = genome[i:i + k]
        k_mer_positions[k_mer].append(i)

    # Calculate distances between repeat positions of each k-mer
    for positions in k_mer_positions.values():
        if len(positions) > 1:
            for j in range(1, len(positions)):
                distances.append(positions[j] - positions[j - 1])

    # Use the average distance as window size L
    if distances:
        avg_distance = sum(distances) // len(distances)
        return avg_distance
    else:
        return None  # Not enough repeats found to estimate L


# Function to estimate t based on the number of repeats in each window
def estimate_t(genome, k, L):
    repetition_counts = []

    # Slide windows of length L
    for i in range(len(genome) - L + 1):
        window = genome[i:i + L]
        frequency_map = collections.defaultdict(int)

        # Count k-mer repeats within the window
        for j in range(len(window) - k + 1):
            k_mer = window[j:j + k]
            frequency_map[k_mer] += 1

        # Get the maximum number of repeats of any k-mer in the window
        max_repetition = max(frequency_map.values(), default=0)
        repetition_counts.append(max_repetition)

    # Use the average of the maximum repeats as the estimate for t
    if repetition_counts:
        return sum(repetition_counts) // len(repetition_counts)
    else:
        return None  # No windows to estimate t


# Main function to estimate k, L, and t
def estimate_clump_parameters(genome):
    # Estimate k
    k, diversity = estimate_k(genome)
    print(f"Estimated optimal k: {k}")
    print(f"Diversity of k-mers for each k: {diversity}")

    # Estimate L
    L = estimate_L(genome, k)
    if L:
        print(f"Estimated window size L: {L}")
    else:
        print("Could not estimate L due to insufficient repeats.")

    # Estimate t
    if L:
        t = estimate_t(genome, k, L)
        if t:
            print(f"Estimated repetition threshold t: {t}")
        else:
            print("Could not estimate t due to insufficient repetitions.")
    else:
        t = None

    return k, L, t


# Example usage
genome = ""  # Insert the genome sequence here
k, L, t = estimate_clump_parameters(genome)

print(f"Estimated parameters: k = {k}, L = {L}, t = {t}")
