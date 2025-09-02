# exercises/randomized_motif_search.py
# Run:
#   python exercises/randomized_motif_search.py
#
# Implements RandomizedMotifSearch with pseudocounts.
# Runs the algorithm 1000 times and prints the best motifs found.

import argparse
import random
from typing import List, Tuple

DNA_ALPHABET = "ACGT"

# Default dataset path (update this to your repo dataset path)
DEFAULT_INPUT_PATH = "data/raw/Randomized_Motif_Search/dataset_30307_5.txt"

def hamming_distance(s1: str, s2: str) -> int:
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def profile_with_pseudocounts(motifs: List[str]) -> dict:
    """Build a profile matrix with pseudocounts (Laplace smoothing)."""
    k = len(motifs[0])
    profile = {base: [1] * k for base in DNA_ALPHABET}  # pseudocounts start at 1
    t = len(motifs)

    for motif in motifs:
        for j, base in enumerate(motif):
            profile[base][j] += 1

    # Normalize to probabilities
    for j in range(k):
        col_sum = sum(profile[base][j] for base in DNA_ALPHABET)
        for base in DNA_ALPHABET:
            profile[base][j] /= col_sum

    return profile

def pr_kmer_profile(kmer: str, profile: dict) -> float:
    """Probability of a k-mer given a profile matrix."""
    prob = 1.0
    for j, base in enumerate(kmer):
        prob *= profile[base][j]
    return prob

def profile_most_probable_kmer(text: str, k: int, profile: dict) -> str:
    """Return the profile-most probable k-mer in text."""
    best_prob = -1
    best_kmer = text[:k]
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        prob = pr_kmer_profile(kmer, profile)
        if prob > best_prob:
            best_prob = prob
            best_kmer = kmer
    return best_kmer

def score(motifs: List[str]) -> int:
    """Score motifs: total number of mismatches to the consensus."""
    k = len(motifs[0])
    t = len(motifs)
    score = 0
    for j in range(k):
        col = [motif[j] for motif in motifs]
        max_count = max(col.count(base) for base in DNA_ALPHABET)
        score += t - max_count
    return score

def randomized_motif_search(dna: List[str], k: int, t: int) -> List[str]:
    """One run of Randomized Motif Search."""
    motifs = [dna[i][random.randint(0, len(dna[i]) - k):][:k] for i in range(t)]
    best_motifs = motifs[:]

    while True:
        profile = profile_with_pseudocounts(motifs)
        motifs = [profile_most_probable_kmer(seq, k, profile) for seq in dna]
        if score(motifs) < score(best_motifs):
            best_motifs = motifs[:]
        else:
            return best_motifs

def read_input(file_path: str = None) -> Tuple[int, int, List[str]]:
    """Read dataset from file (or DEFAULT_INPUT_PATH if none provided)."""
    path = file_path or DEFAULT_INPUT_PATH
    with open(path) as f:
        tokens = f.read().strip().split()
        k, t = map(int, tokens[:2])
        dna = tokens[2:]
    return k, t, dna

def main():
    ap = argparse.ArgumentParser(description="Randomized Motif Search with pseudocounts.")
    ap.add_argument("--file", "-f", type=str, default=None,
                    help=f"Path to dataset file (default: {DEFAULT_INPUT_PATH})")
    ap.add_argument("--iterations", type=int, default=1000,
                    help="Number of times to repeat RandomizedMotifSearch (default=1000).")
    args = ap.parse_args()

    k, t, dna = read_input(args.file)

    best_motifs = randomized_motif_search(dna, k, t)
    best_score = score(best_motifs)

    for _ in range(args.iterations - 1):
        motifs = randomized_motif_search(dna, k, t)
        s = score(motifs)
        if s < best_score:
            best_score = s
            best_motifs = motifs

    print(" ".join(best_motifs))

if __name__ == "__main__":
    main()
