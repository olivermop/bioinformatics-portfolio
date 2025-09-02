# exercises/find_hidden_motif.py
# run:
# #   python exercises/find_hidden_motif.py
# Finds a common 15-mer across multiple DNA strings by k-mer set intersection.
# Prints the motif(s) found and their 0-based and 1-based positions in each string.

from typing import List, Set

def kmers(s: str, k: int) -> Set[str]:
    return {s[i:i+k] for i in range(len(s) - k + 1)}

def positions(s: str, pat: str) -> List[int]:
    k = len(pat)
    return [i for i in range(len(s) - k + 1) if s[i:i+k] == pat]

def find_common_kmers(seqs: List[str], k: int) -> Set[str]:
    common = kmers(seqs[0], k)
    for s in seqs[1:]:
        common &= kmers(s, k)
        if not common:
            break
    return common

def main():
    seqs = [
        "atgaccgggatactgataaaaaaaagggggggggcgtacacattagataaacgtatgaagtacgttagactcggcgccgccg",
        "acccctattttttgagcagatttagtgacctggaaaaaaaatttgagtacaaaacttttccgaataaaaaaaaaggggggga",
        "tgagtatccctgggatgacttaaaaaaaagggggggtgctctcccgatttttgaatatgtaggatcattcgccagggtccga",
        "gctgagaattggatgaaaaaaaagggggggtccacgcaatcgcgaaccaacgcggacccaaaggcaagaccgataaaggaga",
        "tcccttttgcggtaatgtgccgggaggctggttacgtagggaagccctaacggacttaataaaaaaaagggggggcttatag",
        "gtcaatcatgttcttgtgaatggatttaaaaaaaaggggggggaccgcttggcgcacccaaattcagtgtgggcgagcgcaa",
        "cggttttggcccttgttagaggcccccgtaaaaaaaagggggggcaattatgagagagctaatctatcgcgtgcgtgttcat",
        "aacttgagttaaaaaaaagggggggctggggcacatacaagaggagtcttccttatcagttaatgctgtatgacactatgta",
        "ttggcccattggctaaaagcccaacttgacaaatggaagatagaatccttgcataaaaaaaagggggggaccgaaagggaag",
        "ctggtgagcaacgacagattcttacgtgcattagctcgcttccggggatctaatagcacgaagcttaaaaaaaaggggggga",
    ]

    k = 15
    common = find_common_kmers(seqs, k)
    if not common:
        print("No common 15-mer found.")
        return

    print(f"Common {k}-mer(s) found across all sequences:")
    for motif in sorted(common):
        print("  ", motif)

    print("\nPositions by sequence (0-based and 1-based):")
    for idx, s in enumerate(seqs, start=1):
        for motif in sorted(common):
            pos0 = positions(s, motif)
            pos1 = [p + 1 for p in pos0]
            print(f"Seq {idx:2d}: {motif} -> idx0={pos0}, idx1={pos1}")

if __name__ == "__main__":
    main()
