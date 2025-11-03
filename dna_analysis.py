"""
dna_analysis.py — Part A, Q2: Nucleotide counting & analysis utilities.

Run:
    $ python dna_analysis.py

Try:
  • Single sequence: ATGCA
  • Compare: ATGCA  and  AAAAA
"""

from __future__ import annotations
from typing import Dict
from dna_core import DNA, InvalidNucleotideError


def count_nucleotides(seq: str) -> Dict[str, int]:
    """Count A/T/G/C. O(n)."""
    dna = DNA(seq)
    return dna.counts()


def nucleotide_frequencies(seq: str) -> Dict[str, float]:
    """Frequencies in percent. O(n)."""
    dna = DNA(seq)
    return dna.frequencies()


def analysis_report(seq: str, name: str | None = None) -> str:
    """Formatted report."""
    dna = DNA(seq, name=name)
    c = dna.counts()
    f = dna.frequencies()
    gc = dna.gc_content()
    at = dna.at_content()
    lines = [
        f"=== Nucleotide Analysis Report — {dna.name or '<unnamed>'} ===",
        f"Length: {len(dna)}",
        f"Counts : A={c['A']}  T={c['T']}  G={c['G']}  C={c['C']}",
        "Freq % : "
        + f"A={f['A']:.2f}%  T={f['T']:.2f}%  G={f['G']:.2f}%  C={f['C']:.2f}%",
        f"GC%   : {gc:.2f}%",
        f"AT%   : {at:.2f}%",
    ]
    return "\n".join(lines)


def compare_composition(seq1: str, seq2: str) -> Dict[str, float]:
    """
    Compare nucleotide composition between two sequences.
    Returns per-nucleotide % diffs (seq1 - seq2) + GC% diff + L1 distance.
    """
    d1, d2 = DNA(seq1), DNA(seq2)
    f1, f2 = d1.frequencies(), d2.frequencies()
    diffs = {
        "ΔA%": f1["A"] - f2["A"],
        "ΔT%": f1["T"] - f2["T"],
        "ΔG%": f1["G"] - f2["G"],
        "ΔC%": f1["C"] - f2["C"],
        "ΔGC%": d1.gc_content() - d2.gc_content(),
    }
    l1 = abs(diffs["ΔA%"]) + abs(diffs["ΔT%"]) + abs(diffs["ΔG%"]) + abs(diffs["ΔC%"])
    diffs["L1_distance%"] = l1
    return diffs


# -------- CLI --------
if __name__ == "__main__":
    print("Example inputs:")
    print("  Single: ATGCA")
    print("  Compare: ATGCA  and  AAAAA\n")

    while True:
        print("Choose:")
        print("  1) Analyze one sequence")
        print("  2) Compare two sequences")
        print("  0) Exit")
        choice = input("Enter 1/2/0: ").strip()

        if choice == "1":
            name = input("Name/label (optional): ").strip() or None
            seq = input("Enter DNA (A/T/G/C only) — e.g., ATGCA: ").strip()
            try:
                print("\n" + analysis_report(seq, name=name) + "\n")
            except InvalidNucleotideError as e:
                print("\nERROR:", e)
                if getattr(e, "invalid_positions", None):
                    print("Invalid positions:", e.invalid_positions)
                print()
        elif choice == "2":
            seq1 = input("Enter first DNA (A/T/G/C only): ").strip()
            seq2 = input("Enter second DNA (A/T/G/C only): ").strip()
            try:
                diffs = compare_composition(seq1, seq2)
                print("\n=== Composition Differences (seq1 - seq2) ===")
                for k in ["ΔA%", "ΔT%", "ΔG%", "ΔC%", "ΔGC%", "L1_distance%"]:
                    print(f"{k}: {diffs[k]:.2f}")
                print()
            except InvalidNucleotideError as e:
                print("\nERROR:", e)
                if getattr(e, "invalid_positions", None):
                    print("Invalid positions:", e.invalid_positions)
                print()
        elif choice == "0":
            break
        else:
            print("Please enter 1/2/0.\n")
