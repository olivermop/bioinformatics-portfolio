# exercises/greedy_motif_search.py
# Run:
#   python exercises/greedy_motif_search.py

"""
Greedy Motif Search (with deterministic tie-breaking).

Input format (grader-friendly):
- First two integers: k t
- Then a space-separated collection of t DNA strings (may span lines)

Output:
- BestMotifs as space-separated strings (or newline-separated with --newline)

Notes:
- By default, builds a profile from current motifs without pseudocounts (classic version).
- You can enable Laplace pseudocounts (+1) with --pseudocounts for robustness.
- Tie-breaking: when multiple Profile-most probable k-mers have the same probability,
  we return the first occurring k-mer in the text (leftmost index).
"""

import argparse
from typing import List, Dict

DEFAULT_DATASET = "data/raw/Greedy_Motif_Search/dataset_30305_5.txt"
ALPHABET = "ACGT"

def hamming_distance(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))

def score(motifs: List[str]) -> int:
    """Score = sum over columns of (t - max column count). Lower is better."""
    if not motifs:
        return 0
    t = len(motifs)
    k = len(motifs[0])
    sc = 0
    for j in range(k):
        counts = {b: 0 for b in ALPHABET}
        for m in motifs:
            counts[m[j]] += 1
        sc += t - max(counts.values())
    return sc

def build_profile(motifs: List[str], pseudocounts: bool = False) -> Dict[str, List[float]]:
    """Return profile[b][j] as probabilities per column (A,C,G,T)."""
    t = len(motifs)
    k = len(motifs[0])
    add = 1 if pseudocounts else 0
    denom = t + 4 * add
    profile = {b: [0.0] * k for b in ALPHABET}
    for j in range(k):
        counts = {b: add for b in ALPHABET}
        for m in motifs:
            counts[m[j]] += 1
        for b in ALPHABET:
            profile[b][j] = counts[b] / denom
    return profile

def pr_kmer_given_profile(kmer: str, profile: Dict[str, List[float]]) -> float:
    p = 1.0
    for j, ch in enumerate(kmer):
        p *= profile[ch][j]
        if p == 0.0:
            return 0.0
    return p

def profile_most_probable_kmer(text: str, k: int, profile: Dict[str, List[float]]) -> str:
    """Return the highest-probability k-mer in text under profile; ties -> leftmost."""
    best_kmer = text[0:k]
    best_p = -1.0
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        p = pr_kmer_given_profile(kmer, profile)
        if p > best_p:
            best_p = p
            best_kmer = kmer
            # ties keep the first occurrence automatically (strict '>' here)
    return best_kmer

def greedy_motif_search(dna: List[str], k: int, t: int, pseudocounts: bool = False) -> List[str]:
    """Classic Greedy Motif Search across t strings, seeding with each k-mer of dna[0]."""
    best_motifs = [s[:k] for s in dna]  # initial naive motifs
    best_score = score(best_motifs)

    first = dna[0]
    for i in range(len(first) - k + 1):
        motifs = [first[i:i+k]]
        for r in range(1, t):
            prof = build_profile(motifs, pseudocounts=pseudocounts)
            next_kmer = profile_most_probable_kmer(dna[r], k, prof)
            motifs.append(next_kmer)
        sc = score(motifs)
        if sc < best_score:
            best_score = sc
            best_motifs = motifs
    return best_motifs

def parse_input_tokens(tokens: List[str]):
    """
    tokens[0]=k, tokens[1]=t, tokens[2:]=Dna strings (space-separated; may span multiple lines).
    """
    if len(tokens) < 3:
        raise ValueError("Input must provide: k t and at least one DNA string.")
    k = int(tokens[0])
    t = int(tokens[1])
    dna = [tok.strip().upper() for tok in tokens[2:] if tok.strip()]
    if len(dna) != t:
        # Some graders split lines oddly—do a soft check but proceed if >= t
        if len(dna) < t:
            raise ValueError(f"Expected {t} DNA strings, got {len(dna)}.")
        dna = dna[:t]
    return k, t, dna

def read_input(path: str):
    with open(path, "r") as f:
        tokens = f.read().strip().split()
    return parse_input_tokens(tokens)

def main():
    ap = argparse.ArgumentParser(description="Greedy Motif Search (grader-friendly).")
    ap.add_argument("--file", "-f", type=str, default=DEFAULT_DATASET,
                    help=f"Input file path (default: {DEFAULT_DATASET})")
    ap.add_argument("--pseudocounts", action="store_true",
                    help="Use Laplace pseudocounts (+1) when building profiles.")
    ap.add_argument("--newline", action="store_true",
                    help="Print motifs one per line instead of space-separated.")
    args = ap.parse_args()

    k, t, dna = read_input(args.file)
    best = greedy_motif_search(dna, k, t, pseudocounts=args.pseudocounts)

    if args.newline:
        print("\n".join(best))
    else:
        print(" ".join(best))

if __name__ == "__main__":
    main()
