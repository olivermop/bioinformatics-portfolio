# exercises/gibbs_sampler.py
# Run:
#   python exercises/gibbs_sampler.py

#!/usr/bin/env python3
"""
Gibbs Sampler for Motif Finding (with pseudocounts and weighted sampling).

Assignment alignment:
  Input:  integers k, t, N, followed by t DNA strings (space-separated; may span lines).
  Output: BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts.
  Notes:  Use +1 Laplace pseudocounts; select k-mers via probability proportional to Pr(k-mer|Profile).
          We allow --iters to override N if you want to experiment, but by default we use N from file.

Algorithm (one run):
  - Randomly initialize one k-mer motif per string.
  - Repeat N iterations:
      * Pick a random index i in [0..t-1].
      * Build a profile from all motifs except Motifs[i] (with pseudocounts).
      * Replace Motifs[i] by a k-mer sampled from Dna[i] with probability ~ Pr(k-mer|Profile).
      * Track the best motif set by the standard Score (lower is better).
  - Return the best motif set seen in this run.
Then do 20 independent random restarts and return the best across runs.
"""

import argparse
import random
from typing import Dict, List, Tuple

DEFAULT_INPUT_PATH = "data/raw/Gibbs_Sampler/dataset_30309_11.txt"
ALPHABET = "ACGT"

# ---------------- Utilities ----------------

def hamming_distance(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))

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
    Build a profile with Laplace pseudocounts (+1).
    If r = len(motifs), denom per column is r + 4.
    """
    r = len(motifs)
    k = len(motifs[0])
    denom = r + 4
    profile = {b: [0.0] * k for b in ALPHABET}
    for j in range(k):
        counts = {b: 1 for b in ALPHABET}  # +1 pseudocounts
        for m in motifs:
            counts[m[j]] += 1
        for b in ALPHABET:
            profile[b][j] = counts[b] / denom
    return profile

def pr_kmer_given_profile(kmer: str, profile: Dict[str, List[float]]) -> float:
    """Product of per-column probabilities for kmer under profile."""
    p = 1.0
    for j, ch in enumerate(kmer):
        p *= profile[ch][j]
        if p == 0.0:
            return 0.0
    return p

def all_kmer_probs_in_text(text: str, k: int, profile: Dict[str, List[float]]) -> List[float]:
    """Return Pr for each overlapping k-mer in text under profile (not normalized)."""
    n = len(text)
    probs = []
    for i in range(n - k + 1):
        probs.append(pr_kmer_given_profile(text[i:i+k], profile))
    return probs

def weighted_index(weights: List[float]) -> int:
    """
    Return an index i with probability proportional to weights[i].
    If all weights are zero, fall back to uniform random choice.
    """
    total = sum(weights)
    if total <= 0.0:
        return random.randrange(len(weights))
    r = random.random() * total  # uniform in [0, total)
    acc = 0.0
    for i, w in enumerate(weights):
        acc += w
        if r <= acc:
            return i
    return len(weights) - 1  # numerical safeguard

def sample_kmer_from_profile(text: str, k: int, profile: Dict[str, List[float]]) -> str:
    """Sample a k-mer from text proportional to Pr(k-mer | profile)."""
    weights = all_kmer_probs_in_text(text, k, profile)
    idx = weighted_index(weights)
    return text[idx:idx + k]

def random_initial_motifs(dna: List[str], k: int) -> List[str]:
    """Pick one random k-mer (uniform position) from each string."""
    motifs = []
    for s in dna:
        if len(s) < k:
            raise ValueError("All DNA strings must have length >= k.")
        start = random.randint(0, len(s) - k)
        motifs.append(s[start:start + k])
    return motifs

# ---------------- Gibbs Sampler ----------------

def gibbs_sampler_once(dna: List[str], k: int, t: int, iters: int) -> List[str]:
    """One GibbsSampler run with random initialization for exactly 'iters' steps."""
    motifs = random_initial_motifs(dna, k)
    best = motifs[:]
    best_sc = score(best)

    for _ in range(iters):
        i = random.randrange(t)
        # Build profile from all motifs except i
        excluded = [m for j, m in enumerate(motifs) if j != i]
        profile = build_profile_with_pseudocounts(excluded)
        # Resample the i-th motif according to the current profile
        motifs[i] = sample_kmer_from_profile(dna[i], k, profile)
        sc = score(motifs)
        if sc < best_sc:
            best = motifs[:]
            best_sc = sc

    return best

def gibbs_sampler(dna: List[str], k: int, t: int, iters: int, restarts: int = 20) -> List[str]:
    """
    Multiple restarts (default 20) to improve robustness.
    Returns the best motif set across all runs by Score().
    """
    best_overall = gibbs_sampler_once(dna, k, t, iters)
    best_score = score(best_overall)

    for _ in range(restarts - 1):
        candidate = gibbs_sampler_once(dna, k, t, iters)
        sc = score(candidate)
        if sc < best_score:
            best_overall = candidate
            best_score = sc
    return best_overall

# ---------------- I/O ----------------

def parse_tokens(tokens: List[str]) -> Tuple[int, int, int, List[str]]:
    """
    tokens[0]=k, tokens[1]=t, tokens[2]=N, tokens[3:]=t DNA strings (space-separated; may span lines).
    """
    if len(tokens) < 4:
        raise ValueError("Input must provide: k t N and at least one DNA string.")
    k = int(tokens[0]); t = int(tokens[1]); N = int(tokens[2])
    dna = [tok.strip().upper() for tok in tokens[3:] if tok.strip()]
    if len(dna) < t:
        raise ValueError(f"Expected {t} DNA strings, got {len(dna)}.")
    if len(dna) > t:
        dna = dna[:t]
    if any(len(s) < k for s in dna):
        raise ValueError("All DNA strings must have length >= k.")
    return k, t, N, dna

def read_input(path: str = None) -> Tuple[int, int, int, List[str]]:
    p = path or DEFAULT_INPUT_PATH
    with open(p, "r") as f:
        tokens = f.read().strip().split()
    return parse_tokens(tokens)

def main():
    ap = argparse.ArgumentParser(description="Gibbs Sampler for Motif Finding (with pseudocounts).")
    ap.add_argument("--file", "-f", type=str, default=None,
                    help=f"Input file path (default: {DEFAULT_INPUT_PATH})")
    ap.add_argument("--iters", type=int, default=None,
                    help="Override number of iterations per run (otherwise uses N from dataset)")
    ap.add_argument("--restarts", type=int, default=20,
                    help="Number of independent random starts (default: 20)")
    ap.add_argument("--seed", type=int, default=None,
                    help="Random seed for reproducibility")
    ap.add_argument("--newline", action="store_true",
                    help="Print motifs one per line instead of space-separated")
    args = ap.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    k, t, N, dna = read_input(args.file)
    iters = args.iters if args.iters is not None else N

    best = gibbs_sampler(dna, k, t, iters=iters, restarts=args.restarts)

    if args.newline:
        print("\n".join(best))
    else:
        print(" ".join(best))

if __name__ == "__main__":
    main()
