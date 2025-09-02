# exercises/profile_probability_tools.py
# Run:
#   python exercises/profile_probability_tools.py

"""
Profile tools for computing Pr(pattern | profile), finding the most-probable k-mer,
and sampling k-mers by "rolling the profile dice".

Behavior:
- If run with no flags, uses the embedded example profile and prints probabilities
  for two patterns: TCGGGGATTTCC and TCGTGGATTTCC.

Profile file format (rows A,C,G,T; k columns with probabilities summing to ~1 per column):
A: 0.2 0.1 0.0 0.7
C: 0.3 0.3 0.2 0.2
G: 0.1 0.5 0.7 0.0
T: 0.4 0.1 0.1 0.1
"""

import argparse
import math
import random
from typing import Dict, List, Tuple

ALPHABET = "ACGT"

# ---------------------------
# Embedded example profile
# ---------------------------
# Consensus: T C G G G G A T T T C C
# Consensus per-column probs:
# [0.7, 0.6, 1.0, 1.0, 0.9, 0.9, 0.9, 0.5, 0.8, 0.7, 0.4, 0.6]
# Remaining mass split equally across the other three bases (0 if consensus prob = 1.0).
def _build_embedded_example_profile() -> Dict[str, List[float]]:
    consensus = "TCGGGGATTTCC"
    cons_p = [0.7, 0.6, 1.0, 1.0, 0.9, 0.9, 0.9, 0.5, 0.8, 0.7, 0.4, 0.6]
    k = len(consensus)
    profile = {b: [0.0] * k for b in ALPHABET}
    for j, (b_star, p_star) in enumerate(zip(consensus, cons_p)):
        remaining = max(0.0, 1.0 - p_star)
        for b in ALPHABET:
            if b == b_star:
                profile[b][j] = p_star
            else:
                profile[b][j] = remaining / 3.0 if remaining > 0 else 0.0
    return profile

# ---------------------------
# Profile parsing
# ---------------------------
def _parse_profile_from_lines(lines: List[str]) -> Dict[str, List[float]]:
    rows: Dict[str, List[float]] = {}
    for raw in lines:
        s = raw.strip()
        if not s:
            continue
        parts = s.replace("\t", " ").split()
        head = parts[0].rstrip(":").upper()
        if head in {"A", "C", "G", "T"} and (len(parts) > 1):
            key = head
            vals = [float(x) for x in parts[1:]]
        else:
            # unlabeled row; assign in A,C,G,T order
            key = None
            vals = [float(x) for x in parts]
            for b in ALPHABET:
                if b not in rows:
                    key = b
                    break
            if key is None:
                raise ValueError("Too many unlabeled rows; expected exactly 4.")
        if key in rows:
            raise ValueError(f"Duplicate row for base {key}.")
        rows[key] = vals

    if set(rows.keys()) != set(ALPHABET):
        raise ValueError(f"Profile must include rows for A,C,G,T. Found: {sorted(rows.keys())}")
    lens = {len(v) for v in rows.values()}
    if len(lens) != 1:
        raise ValueError("All rows must have the same number of columns (k).")
    k = lens.pop()

    for j in range(k):
        col_sum = sum(rows[b][j] for b in ALPHABET)
        if not (abs(col_sum - 1.0) <= 1e-6 or (0.999 <= col_sum <= 1.001)):
            raise ValueError(f"Column {j} probs sum to {col_sum:.6f}, not ~1.0")

    return rows

def load_profile(path: str) -> Dict[str, List[float]]:
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    return _parse_profile_from_lines(lines)

# ---------------------------
# Core computations
# ---------------------------
def prob_of_pattern_given_profile(pattern: str, profile: Dict[str, List[float]]) -> float:
    """Compute Pr(pattern | profile) = product over j of profile[pattern[j]][j]."""
    pattern = pattern.upper()
    k = len(next(iter(profile.values())))
    if len(pattern) != k:
        raise ValueError(f"Pattern length ({len(pattern)}) must equal profile width ({k}).")
    p = 1.0
    for j, ch in enumerate(pattern):
        if ch not in ALPHABET:
            raise ValueError(f"Invalid base '{ch}' at pos {j}. Allowed: A,C,G,T.")
        p *= profile[ch][j]
        if p == 0.0:
            return 0.0
    return p

def log2_prob_of_pattern(pattern: str, profile: Dict[str, List[float]]) -> float:
    """Return log2 Pr(pattern | profile). Uses -inf if any column prob is 0."""
    pattern = pattern.upper()
    k = len(next(iter(profile.values())))
    if len(pattern) != k:
        raise ValueError(f"Pattern length ({len(pattern)}) must equal profile width ({k}).")
    logp = 0.0
    for j, ch in enumerate(pattern):
        p = profile[ch][j]
        if p <= 0.0:
            return float("-inf")
        logp += math.log2(p)
    return logp

def consensus_from_profile(profile: Dict[str, List[float]]) -> str:
    """Return the consensus string (argmax per column)."""
    k = len(next(iter(profile.values())))
    out = []
    for j in range(k):
        best_b = max(ALPHABET, key=lambda b: profile[b][j])
        out.append(best_b)
    return "".join(out)

def most_probable_kmer_in_text(text: str, k: int, profile: Dict[str, List[float]]) -> Tuple[str, float]:
    """Return the k-mer in 'text' with the highest Pr(kmer | profile), and its probability."""
    text = text.upper()
    if len(text) < k:
        raise ValueError("Text shorter than k.")
    best_kmer, best_p = text[0:k], -1.0
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        p = 1.0
        for j, ch in enumerate(kmer):
            p *= profile[ch][j]
            if p == 0.0:
                break
        if p > best_p:
            best_p = p
            best_kmer = kmer
    return best_kmer, best_p

def sample_kmer_from_profile(profile: Dict[str, List[float]]) -> str:
    """Sample a k-mer by rolling the profile dice (independent columns)."""
    k = len(next(iter(profile.values())))
    out = []
    for j in range(k):
        r = random.random()
        acc = 0.0
        chosen = None
        for b in ALPHABET:
            acc += profile[b][j]
            if r <= acc:
                chosen = b
                break
        if chosen is None:  # numerical safeguard
            chosen = max(ALPHABET, key=lambda b: profile[b][j])
        out.append(chosen)
    return "".join(out)

# ---------------------------
# CLI
# ---------------------------
def main():
    ap = argparse.ArgumentParser(description="Compute Pr(pattern | profile), most-probable k-mer, and sampling.")
    ap.add_argument("--profile", type=str, default=None, help="Path to profile file (4 rows: A,C,G,T).")
    ap.add_argument("--example-profile", action="store_true",
                    help="Use the embedded example profile for consensus TCGGGGATTTCC.")
    ap.add_argument("--pattern", type=str, default=None, help="Pattern to score under the profile.")
    ap.add_argument("--show-log2", action="store_true", help="Also print log2 probability.")
    ap.add_argument("--most-probable-in", type=str, default=None,
                    help="Text in which to find the most-probable k-mer.")
    ap.add_argument("--k", type=int, default=None, help="k for --most-probable-in (must match profile width).")
    ap.add_argument("--sample", type=int, default=0, help="Number of k-mers to sample from the profile.")
    ap.add_argument("--seed", type=int, default=None, help="Random seed for reproducible sampling.")
    args = ap.parse_args()

    # Default: if no profile is provided, use the embedded example profile.
    if args.profile and args.example_profile:
        raise SystemExit("Provide either --profile or --example-profile, not both.")
    if args.profile:
        profile = load_profile(args.profile)
    else:
        # Either --example-profile or no profile flag at all → use embedded example
        profile = _build_embedded_example_profile()

    # If the user provides no actions, show a useful default demo
    if not args.pattern and not args.most_probable_in and args.sample == 0:
        demo_patterns = ["TCGGGGATTTCC", "TCGTGGATTTCC"]
        cons = consensus_from_profile(profile)
        print(f"[Default demo] Using embedded example profile. Consensus: {cons}")
        for pat in demo_patterns:
            p = prob_of_pattern_given_profile(pat, profile)
            print(f"Pr({pat} | profile) = {p:.10g}")
        return

    # Optional: set RNG seed
    if args.seed is not None:
        random.seed(args.seed)

    # 1) Probability of a specific pattern
    if args.pattern:
        cons = consensus_from_profile(profile)
        p = prob_of_pattern_given_profile(args.pattern, profile)
        print(f"Pattern: {args.pattern}")
        print(f"Consensus (from profile): {cons}")
        print(f"Pr(pattern | profile) = {p:.10g}")
        if args.show_log2:
            lp = log2_prob_of_pattern(args.pattern, profile)
            print(f"log2 Pr = {lp:.6f}")

    # 2) Most-probable k-mer in a text
    if args.most_probable_in:
        k = args.k
        if k is None:
            raise SystemExit("--k is required with --most-probable-in.")
        profile_width = len(next(iter(profile.values())))
        if k != profile_width:
            raise SystemExit(f"--k must equal profile width ({profile_width})")
        best, bestp = most_probable_kmer_in_text(args.most_probable_in, k, profile)
        print(f"Most-probable k-mer in text: {best}  (Pr={bestp:.10g})")

    # 3) Sampling
    for _ in range(max(0, args.sample)):
        print(sample_kmer_from_profile(profile))

if __name__ == "__main__":
    main()
