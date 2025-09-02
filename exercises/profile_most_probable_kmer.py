# exercises/profile_most_probable_kmer.py
# Run:
#   python exercises/profile_most_probable_kmer.py


"""
Profile-most Probable k-mer

Reads:
  Line 1: Text (DNA string)
  Line 2: k (integer)
  Next 4 lines (or 4*k numbers across lines): A, C, G, T rows of a 4×k profile matrix
             (unlabeled rows assumed to be in A, C, G, T order; labels like 'A:' allowed)

Prints:
  The Profile-most probable k-mer in Text (first occurrence in case of ties).

This script is intentionally minimal to match Coursera/Stepik graders: it prints only the answer.
"""

import argparse
from typing import Dict, List

DEFAULT_DATASET = "data/raw/Profile_Most_Probable/dataset_30305_3.txt"
ALPHABET = "ACGT"

def _chunk(lst: List[float], size: int) -> List[List[float]]:
    return [lst[i:i+size] for i in range(0, len(lst), size)]

def parse_dataset(path: str):
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    if len(lines) < 6:
        # Could still be valid if profile numbers are spread across fewer/more lines;
        # we'll parse flexibly below.
        pass

    text = lines[0].strip().upper()
    k = int(lines[1].strip())

    # Try labeled rows first (A:, C:, G:, T:)
    label_map: Dict[str, List[float]] = {}
    remaining_lines = lines[2:]

    labeled = True
    for ln in remaining_lines:
        parts = ln.replace("\t", " ").split()
        if parts[0].rstrip(":").upper() in {"A","C","G","T"} and len(parts) > 1:
            base = parts[0].rstrip(":").upper()
            vals = [float(x) for x in parts[1:]]
            label_map[base] = label_map.get(base, []) + vals
        else:
            labeled = False
            break

    if labeled and all(b in label_map for b in "ACGT"):
        # Ensure each row has exactly k entries; if more, truncate; if less, error
        prof = {}
        for b in "ACGT":
            if len(label_map[b]) != k:
                raise ValueError(f"Row {b} must have exactly {k} values (got {len(label_map[b])}).")
            prof[b] = label_map[b]
        return text, k, prof

    # Fallback: unlabeled numeric rows in A,C, G, T order (possibly split across lines)
    nums: List[float] = []
    for ln in remaining_lines:
        for tok in ln.split():
            # Skip any accidental labels; accept raw numbers otherwise
            if tok.rstrip(":").upper() in {"A","C","G","T"}:
                continue
            nums.append(float(tok))
    if len(nums) != 4 * k:
        raise ValueError(f"Expected {4*k} profile numbers, got {len(nums)}.")
    rows = _chunk(nums, k)
    prof = {b: row for b, row in zip("ACGT", rows)}
    return text, k, prof

def kmer_prob(kmer: str, profile: Dict[str, List[float]]) -> float:
    p = 1.0
    for j, ch in enumerate(kmer):
        p *= profile[ch][j]
        if p == 0.0:
            return 0.0
    return p

def profile_most_probable_kmer(text: str, k: int, profile: Dict[str, List[float]]) -> str:
    if len(text) < k:
        raise ValueError("Text shorter than k.")
    best_kmer = text[0:k]
    best_p = -1.0
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        p = kmer_prob(kmer, profile)
        if p > best_p:
            best_p = p
            best_kmer = kmer
            # ties automatically keep the first occurrence
    return best_kmer

def main():
    ap = argparse.ArgumentParser(description="Profile-most Probable k-mer (grader-friendly).")
    ap.add_argument("--file", "-f", type=str, default=DEFAULT_DATASET,
                    help=f"Input file path (default: {DEFAULT_DATASET})")
    args = ap.parse_args()

    text, k, profile = parse_dataset(args.file)
    ans = profile_most_probable_kmer(text, k, profile)
    print(ans)

if __name__ == "__main__":
    main()
