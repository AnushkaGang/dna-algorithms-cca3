"""
dna_transcription.py — Part B, Q4: DNA Transcription

Run:
    $ python dna_transcription.py

Features:
- DNA -> RNA (T->U)
- Supports 'coding' and 'template' strands
- Template orientation control (assume input given 5'->3' by default)
- Batch processing: comma-separated input or file path (one sequence per line)

Dependencies:
- dna_core.DNA  (validates ATGC)
"""

from __future__ import annotations
from typing import Iterable, List, Tuple
from dna_core import DNA, InvalidNucleotideError

_COMP = str.maketrans({"A": "T", "T": "A", "G": "C", "C": "G"})

def complement(seq: str) -> str:
    """DNA complement (uppercase). O(n)."""
    return seq.upper().translate(_COMP)

def reverse_complement(seq: str) -> str:
    """DNA reverse complement (uppercase). O(n)."""
    return complement(seq)[::-1]

def dna_to_rna(dna_seq: str) -> str:
    """Convert DNA (A/T/G/C) to RNA by T->U. O(n)."""
    return dna_seq.upper().replace("T", "U")

def transcribe(
    seq: str,
    strand: str = "coding",
    template_orientation: str = "5to3",
) -> str:
    """
    Transcribe DNA to RNA.

    Parameters
    ----------
    seq : str
        DNA sequence (A/T/G/C only for Q4).
    strand : {'coding','template'}
        - 'coding': RNA is same as coding strand with T->U.
        - 'template': RNA is complementary to template, oriented 5'->3'.
    template_orientation : {'5to3','3to5'}
        How the template string is provided.
        - If '5to3' (common in notes), RNA = reverse_complement(template) then T->U.
        - If '3to5', RNA = complement(template) then T->U.

    Returns
    -------
    str : RNA sequence (5'->3').

    Notes
    -----
    Validation uses DNA() and accepts only A/T/G/C in Q4.
    """
    dna = DNA(seq)  # validates
    s = dna.sequence
    if strand == "coding":
        # mRNA sequence matches coding (5'->3') except U instead of T
        return dna_to_rna(s)
    elif strand == "template":
        if template_orientation == "5to3":
            # If template is given 5'->3', mRNA is reverse-complement
            return dna_to_rna(reverse_complement(s))
        elif template_orientation == "3to5":
            # If template is given 3'->5', mRNA is direct complement
            return dna_to_rna(complement(s))
        else:
            raise ValueError("template_orientation must be '5to3' or '3to5'")
    else:
        raise ValueError("strand must be 'coding' or 'template'")

def transcribe_batch(
    seqs: Iterable[str],
    strand: str = "coding",
    template_orientation: str = "5to3",
) -> List[str]:
    """
    Batch transcribe an iterable of sequences. Raises on first invalid sequence.
    """
    out: List[str] = []
    for s in seqs:
        s = s.strip()
        if not s:
            continue
        out.append(transcribe(s, strand=strand, template_orientation=template_orientation))
    return out


# ---------------------- CLI ----------------------
def _menu() -> None:
    print("Examples you can paste:")
    print("  coding strand (5'->3'):  ATGGACT")
    print("  template 5'->3'       :  TACCTGA")
    print("  file mode: give a path to a .txt file (one sequence per line)")
    print()

    while True:
        print("Choose an option:")
        print("  1) Transcribe a SINGLE sequence")
        print("  2) Transcribe MULTIPLE sequences (comma-separated)")
        print("  3) Transcribe from FILE (one sequence per line)")
        print("  0) Exit")
        choice = input("Enter 1/2/3/0: ").strip()

        if choice == "0":
            break
        elif choice in {"1", "2"}:
            strand = (input("Strand (coding/template) [coding]: ").strip().lower() or "coding")
            if strand not in {"coding", "template"}:
                print("Invalid strand. Use 'coding' or 'template'.\n")
                continue

            template_orientation = "5to3"
            if strand == "template":
                template_orientation = (input("Template orientation (5to3/3to5) [5to3]: ").strip().lower() or "5to3")
                if template_orientation not in {"5to3", "3to5"}:
                    print("Invalid orientation. Use '5to3' or '3to5'.\n")
                    continue

            if choice == "1":
                s = input("Enter DNA sequence (A/T/G/C only): ").strip()
                try:
                    rna = transcribe(s, strand=strand, template_orientation=template_orientation)
                    print("\nRNA (5'->3'):", rna, "\n")
                except InvalidNucleotideError as e:
                    print("\nERROR:", e)
                    if getattr(e, "invalid_positions", None):
                        print("Invalid positions:", e.invalid_positions)
                    print()
            else:
                raw = input("Enter comma-separated DNA sequences: ").strip()
                seqs = [x.strip() for x in raw.split(",") if x.strip()]
                try:
                    out = transcribe_batch(seqs, strand=strand, template_orientation=template_orientation)
                    print("\nResults:")
                    for i, (s, r) in enumerate(zip(seqs, out), 1):
                        print(f"{i:>2}. DNA: {s}  ->  RNA: {r}")
                    print()
                except InvalidNucleotideError as e:
                    print("\nERROR:", e)
                    if getattr(e, "invalid_positions", None):
                        print("Invalid positions:", e.invalid_positions)
                    print()
        elif choice == "3":
            path = input("Enter file path (.txt, one sequence per line): ").strip()
            strand = (input("Strand (coding/template) [coding]: ").strip().lower() or "coding")
            if strand not in {"coding", "template"}:
                print("Invalid strand.\n"); continue
            template_orientation = "5to3"
            if strand == "template":
                template_orientation = (input("Template orientation (5to3/3to5) [5to3]: ").strip().lower() or "5to3")
                if template_orientation not in {"5to3", "3to5"}:
                    print("Invalid orientation.\n"); continue
            try:
                with open(path, "r", encoding="utf-8") as f:
                    seqs = [line.strip() for line in f if line.strip()]
                out = transcribe_batch(seqs, strand=strand, template_orientation=template_orientation)
                print("\nResults:")
                for i, (s, r) in enumerate(zip(seqs, out), 1):
                    print(f"{i:>2}. DNA: {s}  ->  RNA: {r}")
                print()
            except FileNotFoundError:
                print("File not found.\n")
            except InvalidNucleotideError as e:
                print("\nERROR:", e)
                if getattr(e, "invalid_positions", None):
                    print("Invalid positions:", e.invalid_positions)
                print()
        else:
            print("Please enter 1/2/3/0.\n")

if __name__ == "__main__":
    _menu()
