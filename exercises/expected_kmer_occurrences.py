# exercises/expected_kmer_occurrences.py
# run:
#   python expected_kmer_occurrences.py
#

import argparse
import random

def expected_occurrences(num_strings: int, n: int, k: int, alphabet: int = 4) -> float:
    """
    Analytical expectation of occurrences of a fixed k-mer across multiple strings.
    Counts overlapping occurrences.
    """
    if k > n:
        return 0.0
    positions = n - k + 1
    p = (1 / alphabet) ** k
    return num_strings * positions * p

def simulate_occurrences(num_strings: int, n: int, k: int, alphabet: int = 4, trials: int = 1000) -> float:
    """
    Monte Carlo simulation: generate random strings and count exact matches for a fixed k-mer.
    Counts overlapping occurrences.
    """
    # Choose an arbitrary target k-mer, e.g., 'A' * k
    target = 'A' * k
    letters = ['A', 'C', 'G', 'T'][:alphabet]

    total = 0
    for _ in range(trials):
        count_trial = 0
        for _ in range(num_strings):
            s = ''.join(random.choice(letters) for _ in range(n))
            # Naive overlapping-count scan
            for i in range(n - k + 1):
                if s[i:i+k] == target:
                    count_trial += 1
        total += count_trial
    return total / trials

def main():
    parser = argparse.ArgumentParser(description="Expected number of occurrences of a k-mer in random strings.")
    parser.add_argument("--m", "--num_strings", type=int, default=500, dest="m",
                        help="Number of strings (default: 500)")
    parser.add_argument("--n", type=int, default=1000,
                        help="Length of each string (default: 1000)")
    parser.add_argument("--k", type=int, default=9,
                        help="k-mer length (default: 9)")
    parser.add_argument("--alphabet", type=int, default=4,
                        help="Alphabet size (DNA=4) (default: 4)")
    parser.add_argument("--simulate", type=int, default=0,
                        help="Number of Monte Carlo trials for validation (0 = skip)")
    args = parser.parse_args()

    exp_val = expected_occurrences(args.m, args.n, args.k, args.alphabet)
    print(f"Analytical expectation: {exp_val:.12f}")

    if args.simulate > 0:
        sim_val = simulate_occurrences(args.m, args.n, args.k, args.alphabet, trials=args.simulate)
        print(f"Monte Carlo estimate ({args.simulate} trials): {sim_val:.12f}")

if __name__ == "__main__":
    main()
