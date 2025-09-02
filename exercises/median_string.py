# exercises/median_string.py
# Run:
#   python exercises/median_string.py

"""
Median String (brute force).

Given an integer k and a collection of DNA strings Dna,
find a k-mer Pattern that minimizes d(Pattern, Dna),
where
    d(Pattern, Dna) = sum over each string s in Dna of
                      min_{k-mer x in s} HammingDistance(Pattern, x)

Implementation notes:
- By default, reads input from DEFAULT_INPUT_PATH below so you can run the script
  without any flags.
- You can still override with --file to point to another dataset.
- Enumerates all 4^k DNA k-mers (A,C,G,T). Exact but exponential—best for small k.
- Prints one optimal pattern by default (any one); with --all prints all ties.
- With --show-score also prints the minimal total distance.

Input file format:
    Line 1: k
    Next tokens (space-separated, across one or more lines): DNA strings
"""

import argparse
import itertools
from typing import List, Tuple

# Default dataset:
DEFAULT_INPUT_PATH = "data/raw/Median_String/dataset_30304_9.txt"
DNA_ALPHABET = "ACGT"

def hamming_distance(a: str, b: str) -> int:
    """Return the Hamming distance between equal-length strings a and b."""
    return sum(x != y for x, y in zip(a, b))

def all_kmers_in_text(text: str, k: int) -> List[str]:
    """Return all overlapping k-mers in 'text'."""
    n = len(text)
    if k > n:
        return []
    return [text[i:i+k] for i in range(n - k + 1)]

def distance_pattern_to_text(pattern: str, text_kmers: List[str]) -> int:
    """
    For a fixed 'pattern' and a precomputed list of k-mers from a text,
    return min Hamming distance between 'pattern' and any k-mer in that text.
    """
    k = len(pattern)
    best = k + 1  # upper bound
    for x in text_kmers:
        d = hamming_distance(pattern, x)
        if d < best:
            best = d
            if best == 0:  # can't get better
                return 0
    return best

def distance_pattern_to_dna(pattern: str, dna_kmers: List[List[str]]) -> int:
    """
    Sum of distances between 'pattern' and each string in Dna,
    where each distance is min over that string's k-mers.
    """
    return sum(distance_pattern_to_text(pattern, kmers) for kmers in dna_kmers)

def median_string(dna: List[str], k: int, alphabet: str = DNA_ALPHABET,
                  return_all: bool = False) -> Tuple[List[str], int]:
    """
    Brute-force Median String solver.
    Returns:
      (patterns, best_score)
    - patterns: one optimal pattern by default; if return_all=True, all ties.
    - best_score: minimal total distance d(Pattern, Dna).
    """
    # Precompute all k-mers for each DNA string to speed up scoring.
    dna_kmers = [all_kmers_in_text(s, k) for s in dna]

    best_score = float("inf")
    best_patterns: List[str] = []

    # Enumerate all k-mers over the alphabet (4^k for DNA).
    for tup in itertools.product(alphabet, repeat=k):
        pat = "".join(tup)
        score = distance_pattern_to_dna(pat, dna_kmers)

        if score < best_score:
            best_score = score
            best_patterns = [pat]
        elif score == best_score and return_all:
            best_patterns.append(pat)

    return best_patterns, best_score

def parse_input_tokens(tokens: List[str]) -> Tuple[int, List[str]]:
    """
    Expected tokens:
      tokens[0] = k (int)
      tokens[1:] = DNA strings (space-separated; may span multiple lines)
    """
    if not tokens:
        raise ValueError("Empty input.")
    k = int(tokens[0])
    dna = [t.strip().upper() for t in tokens[1:] if t.strip()]
    if not dna:
        raise ValueError("No DNA strings provided.")
    return k, dna

def read_input(file_path: str = None) -> Tuple[int, List[str]]:
    """Read input from the given file path (or DEFAULT_INPUT_PATH if None)."""
    path = file_path or DEFAULT_INPUT_PATH
    try:
        with open(path, "r") as f:
            content = f.read().strip().split()
    except FileNotFoundError as e:
        raise SystemExit(
            f"[ERROR] Could not open input file: {path}\n"
            f"Hint: ensure the dataset exists or pass --file PATH"
        ) from e
    return parse_input_tokens(content)

def main():
    ap = argparse.ArgumentParser(description="Median String (brute force) for DNA.")
    ap.add_argument("--file", "-f", type=str, default=None,
                    help=f"Path to input file (default: {DEFAULT_INPUT_PATH}).")
    ap.add_argument("--all", action="store_true",
                    help="Print all tie-optimal patterns (space-separated).")
    ap.add_argument("--show-score", action="store_true",
                    help="Also print the minimal total distance.")
    ap.add_argument("--alphabet", type=str, default=DNA_ALPHABET,
                    help="Alphabet to use (default: ACGT).")
    args = ap.parse_args()

    k, dna = read_input(args.file)
    patterns, best_score = median_string(dna, k, alphabet=args.alphabet, return_all=args.all)

    if args.all:
        print(" ".join(sorted(patterns)))
    else:
        print(patterns[0])

    if args.show_score:
        print(f"# minimal total distance: {best_score}")

if __name__ == "__main__":
    main()
