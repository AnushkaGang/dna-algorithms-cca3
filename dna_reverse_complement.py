"""
dna_reverse_complement.py — Part B, Q5: Reverse Complement with IUPAC degenerates.

Run:
    $ python dna_reverse_complement.py

What it supports
- Complement and Reverse-Complement
- IUPAC degenerate nucleotides:
  A, T, G, C, R, Y, S, W, K, M, B, D, H, V, N
- Input orientation hint (5'->3' or 3'->5') and desired output orientation
- Batch via comma-separated or file mode
- Clear validation with invalid positions

Complexity
- O(n) per sequence; translate() and slicing are C-optimized
"""

from __future__ import annotations
from typing import Iterable, List, Tuple

# Allowed IUPAC DNA symbols (uppercase)
IUPAC_SET = set("ATGCRYSWKMBDHVN")

# IUPAC complements (uppercase)
# A<->T, C<->G, R<->Y, S<->S, W<->W, K<->M, B<->V, D<->H, H<->D, V<->B, N<->N
_COMP_MAP = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
}
_COMP_TABLE = str.maketrans(_COMP_MAP)

class InvalidSymbolError(ValueError):
    """Raised when sequence contains invalid (non-IUPAC DNA) characters."""
    def __init__(self, message: str, invalid_positions: List[int] | None = None):
        super().__init__(message)
        self.invalid_positions = invalid_positions or []

def _validate_iupac(seq: str) -> str:
    """Uppercase and validate against IUPAC_SET. Raises on first invalid."""
    s = seq.upper()
    bad = [i + 1 for i, ch in enumerate(s) if ch not in IUPAC_SET]
    if bad:
        bad_chars = [s[i - 1] for i in bad]
        raise InvalidSymbolError(
            f"Invalid symbols at positions {bad} -> {bad_chars}. "
            "Allowed: A,T,G,C,R,Y,S,W,K,M,B,D,H,V,N",
            bad
        )
    return s

def complement_iupac(seq: str) -> str:
    """Return complement (same orientation as input). O(n)."""
    s = _validate_iupac(seq)
    return s.translate(_COMP_TABLE)

def reverse_complement_iupac(seq: str) -> str:
    """Return reverse complement (orientation reversed). O(n)."""
    return complement_iupac(seq)[::-1]

def rc_with_orientations(
    seq: str,
    op: str = "revcomp",           # 'comp' or 'revcomp'
    input_orientation: str = "5to3",
    output_orientation: str = "5to3",
) -> str:
    """
    Compute complement/reverse-complement with explicit orientation control.

    Parameters
    ----------
    seq : str
        DNA sequence (IUPAC symbols supported).
    op : {'comp','revcomp'}
        'comp'    -> complement only (keeps order)
        'revcomp' -> reverse-complement (reverses order)
    input_orientation : {'5to3','3to5'}
        How the given string is oriented. The internal result is normalized to 5'->3' first.
    output_orientation : {'5to3','3to5'}
        Desired orientation of the returned string.

    Notes
    -----
    - If input_orientation == '3to5', we reverse to interpret it as 5'->3' before applying op.
    - After op, we reverse again if output_orientation == '3to5'.
    """
    s = _validate_iupac(seq)

    # Normalize input to 5'->3'
    if input_orientation == "3to5":
        s = s[::-1]
    elif input_orientation != "5to3":
        raise ValueError("input_orientation must be '5to3' or '3to5'")

    # Apply operation in 5'->3' space
    if op == "comp":
        r = s.translate(_COMP_TABLE)
    elif op == "revcomp":
        r = s.translate(_COMP_TABLE)[::-1]
    else:
        raise ValueError("op must be 'comp' or 'revcomp'")

    # Adjust desired output orientation
    if output_orientation == "5to3":
        return r
    elif output_orientation == "3to5":
        return r[::-1]
    else:
        raise ValueError("output_orientation must be '5to3' or '3to5'")

def rc_batch(
    seqs: Iterable[str],
    op: str = "revcomp",
    input_orientation: str = "5to3",
    output_orientation: str = "5to3",
) -> List[str]:
    """Batch mode for rc_with_orientations."""
    out: List[str] = []
    for s in seqs:
        s = s.strip()
        if not s:
            continue
        out.append(rc_with_orientations(
            s, op=op,
            input_orientation=input_orientation,
            output_orientation=output_orientation,
        ))
    return out


# ---------------------- CLI ----------------------
def _menu() -> None:
    print("Examples you can paste:")
    print("  Simple:           ATGCRY")
    print("  With degenerates: ATGCNBDHV")
    print("  3'->5' input:     (reverse your 5'->3' string to simulate)\n")

    while True:
        print("Choose an option:")
        print("  1) Complement a SINGLE sequence")
        print("  2) Reverse-complement a SINGLE sequence")
        print("  3) Complement MULTIPLE (comma-separated)")
        print("  4) Reverse-complement MULTIPLE (comma-separated)")
        print("  5) File mode (one sequence per line) — choose comp/revcomp")
        print("  0) Exit")
        choice = input("Enter 1/2/3/4/5/0: ").strip()

        if choice == "0":
            break

        try:
            if choice in {"1", "2"}:
                op = "comp" if choice == "1" else "revcomp"
                seq = input("Enter DNA (IUPAC allowed): ").strip()
                iori = (input("Input orientation (5to3/3to5) [5to3]: ").strip().lower() or "5to3")
                oori = (input("Output orientation (5to3/3to5) [5to3]: ").strip().lower() or "5to3")
                res = rc_with_orientations(seq, op=op, input_orientation=iori, output_orientation=oori)
                print(f"\nResult ({op}):", res, "\n")

            elif choice in {"3", "4"}:
                op = "comp" if choice == "3" else "revcomp"
                raw = input("Enter comma-separated sequences: ").strip()
                seqs = [x.strip() for x in raw.split(",") if x.strip()]
                iori = (input("Input orientation (5to3/3to5) [5to3]: ").strip().lower() or "5to3")
                oori = (input("Output orientation (5to3/3to5) [5to3]: ").strip().lower() or "5to3")
                out = rc_batch(seqs, op=op, input_orientation=iori, output_orientation=oori)
                print("\nResults:")
                for i, (s, r) in enumerate(zip(seqs, out), 1):
                    print(f"{i:>2}. IN : {s}   ->   OUT: {r}")
                print()

            elif choice == "5":
                op_in = (input("Operation (comp/revcomp) [revcomp]: ").strip().lower() or "revcomp")
                iori = (input("Input orientation (5to3/3to5) [5to3]: ").strip().lower() or "5to3")
                oori = (input("Output orientation (5to3/3to5) [5to3]: ").strip().lower() or "5to3")
                path = input("Enter file path (.txt): ").strip()
                with open(path, "r", encoding="utf-8") as f:
                    seqs = [line.strip() for line in f if line.strip()]
                out = rc_batch(seqs, op=op_in, input_orientation=iori, output_orientation=oori)
                print("\nResults:")
                for i, (s, r) in enumerate(zip(seqs, out), 1):
                    print(f"{i:>2}. IN : {s}   ->   OUT: {r}")
                print()

            else:
                print("Please enter 1/2/3/4/5/0.\n")

        except InvalidSymbolError as e:
            print("\nERROR:", e)
            if getattr(e, "invalid_positions", None):
                print("Invalid positions:", e.invalid_positions)
            print()
        except FileNotFoundError:
            print("\nFile not found.\n")
        except ValueError as e:
            print("\nERROR:", e, "\n")

if __name__ == "__main__":
    _menu()
