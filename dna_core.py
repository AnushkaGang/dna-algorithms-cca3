"""
dna_core.py — Part A, Question 1: DNA Representation and Basic Operations.

Run:
    $ python dna_core.py

Sample to try:
    ATGCA
"""

from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
from typing import Dict, List, Tuple

VALID_NUCLEOTIDES: Tuple[str, ...] = ("A", "T", "G", "C")


class InvalidNucleotideError(ValueError):
    """Raised when the DNA sequence contains invalid nucleotide characters."""
    def __init__(self, message: str, invalid_positions: List[int] | None = None):
        super().__init__(message)
        self.invalid_positions = invalid_positions or []


@dataclass(slots=True)
class DNA:
    """
    Represents a DNA sequence with basic utilities (A/T/G/C only for Part A).
    """
    sequence: str
    name: str | None = None

    def __post_init__(self) -> None:
        seq = self.sequence.upper()
        invalid_positions = [i + 1 for i, ch in enumerate(seq) if ch not in VALID_NUCLEOTIDES]
        if invalid_positions:
            bad_chars = [seq[i - 1] for i in invalid_positions]
            raise InvalidNucleotideError(
                "Invalid DNA sequence: found non-ATGC symbols at positions "
                f"{invalid_positions} -> {bad_chars}",
                invalid_positions,
            )
        self.sequence = seq

    def __len__(self) -> int:  # O(1)
        return len(self.sequence)

    def counts(self) -> Dict[str, int]:  # O(n)
        c = Counter(self.sequence)
        return {nuc: int(c.get(nuc, 0)) for nuc in VALID_NUCLEOTIDES}

    def frequencies(self) -> Dict[str, float]:  # O(n)
        n = len(self)
        if n == 0:
            return {nuc: 0.0 for nuc in VALID_NUCLEOTIDES}
        counts = self.counts()
        return {nuc: (counts[nuc] / n) * 100.0 for nuc in VALID_NUCLEOTIDES}

    def gc_content(self) -> float:  # O(n)
        n = len(self)
        if n == 0:
            return 0.0
        c = self.counts()
        return ((c["G"] + c["C"]) / n) * 100.0

    def at_content(self) -> float:  # O(n)
        n = len(self)
        if n == 0:
            return 0.0
        c = self.counts()
        return ((c["A"] + c["T"]) / n) * 100.0

    def basic_stats(self) -> Dict[str, object]:  # O(n)
        return {
            "name": self.name,
            "length": len(self),
            "counts": self.counts(),
            "frequencies": self.frequencies(),
            "GC%": self.gc_content(),
            "AT%": self.at_content(),
        }


def _print_basic_report(dna: DNA) -> None:
    c = dna.counts()
    f = dna.frequencies()
    print(f"\n=== DNA Basic Report — {dna.name or '<unnamed>'} ===")
    print(f"Length : {len(dna)}")
    print(f"Counts : A={c['A']}  T={c['T']}  G={c['G']}  C={c['C']}")
    print("Freq % : "
          f"A={f['A']:.2f}%  T={f['T']:.2f}%  G={f['G']:.2f}%  C={f['C']:.2f}%")
    print(f"GC%    : {dna.gc_content():.2f}%")
    print(f"AT%    : {dna.at_content():.2f}%\n")


if __name__ == "__main__":
    print("Sample input you can paste:  ATGCA")
    print("Note: In Part A we accept only A/T/G/C. (Cleaning comes in Q3.)\n")

    while True:
        name = input("Enter a name/label (or leave blank): ").strip() or None
        seq = input("Enter DNA sequence (A/T/G/C only) — e.g., ATGCA: ").strip()

        try:
            dna = DNA(seq, name=name)
            _print_basic_report(dna)
        except InvalidNucleotideError as e:
            print("\nERROR:", e)
            if getattr(e, "invalid_positions", None):
                print("Invalid positions:", e.invalid_positions)
            print("Please try again with only A/T/G/C.\n")

        again = input("Analyze another sequence? (y/n): ").strip().lower()
        if again != "y":
            break
