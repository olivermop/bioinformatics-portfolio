# exercises/motif_enumeration.py
# Run:
#   python exercises/motif_enumeration.py

"""
Motif Enumeration (brute force / exhaustive search).

Given: integers k and d, followed by a collection of DNA strings (Dna).
Finds all (k, d)-motifs that appear in every string in Dna with at most d mismatches.

Input formats supported:
1) From a file (--file path/to/input.txt), where the first line has "k d"
   and the next line (or the same line) contains space-separated DNA strings.
2) From stdin (same format).

Output:
A single line with space-separated motifs in lexicographic order (duplicates removed).
"""

import argparse
from typing import List, Set

def hamming_distance(a: str, b: str) -> int:
    """Return the Hamming distance between equal-length strings a and b."""
    return sum(x != y for x, y in zip(a, b))

def neighbors(pattern: str, d: int) -> Set[str]:
    """
    Generate the d-neighborhood of 'pattern': all strings whose Hamming distance
    from pattern is <= d (over alphabet {A,C,G,T}).
    """
    if d == 0:
        return {pattern}
    if len(pattern) == 0:
        return {""}

    alphabet = ["A", "C", "G", "T"]
    suffix_neighbors = neighbors(pattern[1:], d)
    neighborhood = set()

    for sn in suffix_neighbors:
        # Count mismatches in suffix
        dist = hamming_distance(pattern[1:], sn)
        if dist < d:
            # We can change the first character to any base
            for base in alphabet:
                neighborhood.add(base + sn)
        else:
            # Must keep the first character
            neighborhood.add(pattern[0] + sn)

    return neighborhood

def appears_with_mismatches(pattern: str, text: str, d: int) -> bool:
    """Check if pattern appears in text with at most d mismatches (overlapping allowed)."""
    k = len(pattern)
    for i in range(len(text) - k + 1):
        if hamming_distance(pattern, text[i:i+k]) <= d:
            return True
    return False

def motif_enumeration(dna: List[str], k: int, d: int) -> Set[str]:
    """
    Brute-force MotifEnumeration:
    - Enumerate all k-mers from the FIRST string only (key optimization).
    - For each k-mer, enumerate its d-neighborhood.
    - Keep those neighbors that appear (<= d mismatches) in every string in dna.
    """
    first = dna[0]
    candidates = set()

    # 1) Collect all neighbors of every k-mer in the first string
    for i in range(len(first) - k + 1):
        kmer = first[i:i+k]
        candidates |= neighbors(kmer, d)

    # 2) Filter candidates by checking presence (<= d mismatches) in all strings
    motifs = set()
    for pat in candidates:
        if all(appears_with_mismatches(pat, s, d) for s in dna):
            motifs.add(pat)

    return motifs

def parse_input_tokens(tokens: List[str]):
    """Parse tokens: first two are k and d, the rest are the DNA strings."""
    if len(tokens) < 3:
        raise ValueError("Input must provide: k d and at least one DNA string.")
    k = int(tokens[0])
    d = int(tokens[1])
    dna = tokens[2:]
    return k, d, dna

def read_input(file_path: str = None):
    """Read and tokenize input from a file or stdin."""
    if file_path:
        with open(file_path, "r") as f:
            content = f.read().strip().split()
    else:
        import sys
        content = sys.stdin.read().strip().split()
    return parse_input_tokens(content)

def main():
    ap = argparse.ArgumentParser(description="Brute-force Motif Enumeration (k, d)-motifs.")
    ap.add_argument("--file", "-f", type=str, default=None,
                    help="Path to input file. If omitted, reads from stdin.")
    args = ap.parse_args()

    k, d, dna = read_input(args.file)
    motifs = motif_enumeration(dna, k, d)
    print(" ".join(sorted(motifs)))

if __name__ == "__main__":
    main()
