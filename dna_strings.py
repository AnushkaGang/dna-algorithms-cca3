"""
dna_strings.py — Part A, Q3: String Manipulation for Genomics.

Run:
    $ python dna_strings.py

Functions:
- to_upper(seq)
- to_lower(seq)
- clean_sequence(seq)            # keep only A/T/G/C
- split_codons(seq, frame=0, drop_incomplete=True)
- merge_fragments(fragments)
"""

from __future__ import annotations
from typing import Iterable, List, Tuple

VALID = set("ATGC")


def to_upper(seq: str) -> str:
    """Convert to uppercase. O(n)."""
    return seq.upper()


def to_lower(seq: str) -> str:
    """Convert to lowercase. O(n)."""
    return seq.lower()


def clean_sequence(seq: str) -> str:
    """
    Remove all non-nucleotide characters, keep only A/T/G/C (uppercase).
    Examples:
        'atg cAx\\n---tgtN' -> 'ATGCATGT'
    O(n) time, O(1) aux space (fixed alphabet).
    """
    seq = seq.upper()
    return "".join(ch for ch in seq if ch in VALID)


def split_codons(seq: str, frame: int = 0, drop_incomplete: bool = True) -> Tuple[List[str], str]:
    """
    Split a DNA string into codons (triplets) starting at frame 0/1/2.
    Returns (codons_list, leftover).

    - frame: 0 means start at first base, 1 means skip first base, 2 means skip first two.
    - drop_incomplete: if True, leftover (<3 bases) returned separately; if False, include leftover as last item.
    """
    if frame not in (0, 1, 2):
        raise ValueError("frame must be 0, 1, or 2")

    s = seq[frame:]
    codons = [s[i:i+3] for i in range(0, len(s) - len(s) % 3, 3)]
    leftover = s[len(codons) * 3:]
    if not drop_incomplete and leftover:
        codons.append(leftover)
        return codons, ""
    return codons, leftover


def merge_fragments(fragments: Iterable[str]) -> str:
    """
    Merge multiple DNA fragments into one cleaned uppercase sequence.
    Each fragment is cleaned (A/T/G/C only) before merge.
    """
    cleaned = (clean_sequence(f) for f in fragments)
    return "".join(cleaned)


# ---------- simple CLI for demo ----------
if __name__ == "__main__":
    print("Try samples you can paste:")
    print("  messy:  atg cAx\\n---tgtN")
    print("  frags:  ATG,  caa-tt,  g\n")

    while True:
        print("Menu:")
        print("  1) Uppercase / Lowercase")
        print("  2) Clean non-ATGC characters")
        print("  3) Split into codons")
        print("  4) Merge fragments")
        print("  0) Exit")
        choice = input("Enter 1/2/3/4/0: ").strip()

        if choice == "1":
            s = input("Enter sequence: ")
            print("UPPER:", to_upper(s))
            print("lower:", to_lower(s))
            print()
        elif choice == "2":
            s = input("Enter possibly messy sequence: ")
            print("CLEANED:", clean_sequence(s))
            print()
        elif choice == "3":
            s = input("Enter sequence (we will clean it first): ")
            cs = clean_sequence(s)
            try:
                frame = int(input("Frame (0/1/2): ").strip() or "0")
            except ValueError:
                frame = 0
            drop_flag = input("Drop incomplete codon? (y/n, default y): ").strip().lower() != "n"
            codons, leftover = split_codons(cs, frame=frame, drop_incomplete=drop_flag)
            print("CODONS:", codons)
            print("LEFTOVER:", leftover if leftover else "<none>")
            print()
        elif choice == "4":
            raw = input("Enter fragments separated by commas: ")
            frags = [f.strip() for f in raw.split(",")]
            merged = merge_fragments(frags)
            print("MERGED (cleaned):", merged)
            print()
        elif choice == "0":
            break
        else:
            print("Please enter 1/2/3/4/0.\n")
