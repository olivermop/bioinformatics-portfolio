# exercises/greedy_motif_search_pseudocounts.py
# Run:
#   python exercises/greedy_motif_search_pseudocounts.py

"""
Greedy Motif Search with pseudocounts (+1 Laplace smoothing).

Problem:
  Input: integers k and t, followed by a space-separated collection of t DNA strings (Dna).
  Output: BestMotifs from applying GreedyMotifSearch(Dna, k, t) with pseudocounts.
  Tie-breaking: if multiple profile-most probable k-mers have equal probability,
                choose the leftmost occurrence in the string.

Notes:
  - Pseudocounts are applied when building the profile from the CURRENT motifs
    (i.e., if we have r motifs so far, denom = r + 4).
  - Score(Motifs) is the usual column score (no pseudocounts in scoring).
  - Prints motifs space-separated by default, or one per line with --newline.

Sample (from the spec):
  Input:
    3 5
    GGCGTTCAGGCA AAGAATCAGTCA CAAGGAGTTCGC CACGTCAATCAC CAATAATATTCG
  Output:
    TTC ATC TTC ATC TTC
"""

import argparse
from typing import Dict, List

DEFAULT_DATASET = "data/raw/Greedy_Motif_Search/dataset_30306_9.txt"
ALPHABET = "ACGT"

# ---------- Utilities ----------

def score(motifs: List[str]) -> int:
    """Score = sum over columns of (t - max column count). Lower is better."""
    if not motifs:
        return 0
    t = len(motifs)
    k = len(motifs[0])
    total = 0
    for j in range(k):
        counts = {b: 0 for b in ALPHABET}
        for m in motifs:
            counts[m[j]] += 1
        total += t - max(counts.values())
    return total

def build_profile_with_pseudocounts(motifs: List[str]) -> Dict[str, List[float]]:
    """
    Build a profile with Laplace pseudocounts (+1) from the CURRENT motifs.
    If r = len(motifs) and motif length = k, then for each column j:
        count[b,j] = (# of base b at j among motifs) + 1
        profile[b,j] = count[b,j] / (r + 4)
    """
    r = len(motifs)
    k = len(motifs[0])
    denom = r + 4  # +1 per base
    profile = {b: [0.0] * k for b in ALPHABET}
    for j in range(k):
        counts = {b: 1 for b in ALPHABET}  # start at 1 (pseudocount)
        for m in motifs:
            counts[m[j]] += 1
        for b in ALPHABET:
            profile[b][j] = counts[b] / denom
    return profile

def pr_kmer_given_profile(kmer: str, profile: Dict[str, List[float]]) -> float:
    """Product of per-column probabilities; independent columns assumption."""
    p = 1.0
    for j, ch in enumerate(kmer):
        p *= profile[ch][j]
        if p == 0.0:
            return 0.0
    return p

def profile_most_probable_kmer(text: str, k: int, profile: Dict[str, List[float]]) -> str:
    """
    Return the highest-probability k-mer in 'text' under 'profile'.
    Tie-breaking: first (leftmost) k-mer wins because we update only on strict '>'.
    """
    best_idx = 0
    best_p = -1.0
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        p = pr_kmer_given_profile(kmer, profile)
        if p > best_p:
            best_p = p
            best_idx = i
    return text[best_idx:best_idx+k]

# ---------- Greedy with pseudocounts ----------

def greedy_motif_search_with_pseudocounts(dna: List[str], k: int, t: int) -> List[str]:
    """
    Greedy Motif Search:
      - Seed with each k-mer of the first string.
      - Build profile from current motifs using pseudocounts.
      - Pick the profile-most probable k-mer in each subsequent string (leftmost on ties).
      - Keep the motif set with minimal score.
    """
    # Start with the naive set: first k-mer from each string
    best_motifs = [s[:k] for s in dna]
    best_score = score(best_motifs)

    first = dna[0]
    for i in range(len(first) - k + 1):
        motifs = [first[i:i+k]]
        for r in range(1, t):
            profile = build_profile_with_pseudocounts(motifs)
            next_kmer = profile_most_probable_kmer(dna[r], k, profile)
            motifs.append(next_kmer)
        sc = score(motifs)
        if sc < best_score:
            best_score = sc
            best_motifs = motifs
    return best_motifs

# ---------- I/O ----------

def parse_input_tokens(tokens: List[str]):
    """
    tokens[0]=k, tokens[1]=t, tokens[2:]=t DNA strings (space-separated; may span lines).
    """
    if len(tokens) < 3:
        raise ValueError("Input must provide: k t and at least one DNA string.")
    k = int(tokens[0])
    t = int(tokens[1])
    dna = [tok.strip().upper() for tok in tokens[2:] if tok.strip()]
    if len(dna) < t:
        raise ValueError(f"Expected {t} DNA strings, got {len(dna)}.")
    if len(dna) > t:
        dna = dna[:t]
    # Basic sanity
    if any(len(s) < k for s in dna):
        raise ValueError("All DNA strings must have length >= k.")
    return k, t, dna

def read_input(path: str):
    with open(path, "r") as f:
        tokens = f.read().strip().split()
    return parse_input_tokens(tokens)

def main():
    ap = argparse.ArgumentParser(description="Greedy Motif Search with pseudocounts (+1).")
    ap.add_argument("--file", "-f", type=str, default=DEFAULT_DATASET,
                    help=f"Input file path (default: {DEFAULT_DATASET})")
    ap.add_argument("--newline", action="store_true",
                    help="Print motifs one per line instead of space-separated.")
    args = ap.parse_args()

    k, t, dna = read_input(args.file)
    best = greedy_motif_search_with_pseudocounts(dna, k, t)

    if args.newline:
        print("\n".join(best))
    else:
        print(" ".join(best))

if __name__ == "__main__":
    main()
