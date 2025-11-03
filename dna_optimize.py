"""
dna_optimize.py — Part C, Q6: Algorithm Optimization & Profiling

Run:
    $ python dna_optimize.py

What you can do:
- Benchmark different nucleotide counting approaches
- Profile a chosen approach with cProfile
- Count nucleotides from a HUGE FASTA/TXT file in a memory-efficient stream

Depends on: dna_core.DNA
"""

from __future__ import annotations
import time
import timeit
import cProfile
import pstats
from collections import Counter
from typing import Dict, Iterable, Tuple
from dna_core import DNA, InvalidNucleotideError

VALID = ("A", "T", "G", "C")


# -------------------- Approaches to compare --------------------

def counts_counter(seq: str) -> Dict[str, int]:
    """Approach A: collections.Counter (fast, readable). O(n) time, O(1) aux."""
    dna = DNA(seq)
    c = Counter(dna.sequence)
    return {n: int(c.get(n, 0)) for n in VALID}

def counts_strcount(seq: str) -> Dict[str, int]:
    """Approach B: str.count per nucleotide (often very fast). O(4n) ~ O(n)."""
    s = DNA(seq).sequence
    return {n: s.count(n) for n in VALID}

def counts_manual_loop(seq: str) -> Dict[str, int]:
    """Approach C: manual loop (baseline). O(n)."""
    s = DNA(seq).sequence
    a = t = g = c = 0
    for ch in s:
        if ch == "A": a += 1
        elif ch == "T": t += 1
        elif ch == "G": g += 1
        elif ch == "C": c += 1
    return {"A": a, "T": t, "G": g, "C": c}


# -------------------- Memory-efficient streaming --------------------

def stream_counts(lines: Iterable[str]) -> Dict[str, int]:
    """
    Count A/T/G/C from an iterable of lines without loading full sequence.
    - Uppercases on the fly.
    - Ignores non-ATGC (Q3 cleaning idea).
    O(total chars) time, O(1) aux.
    """
    a = t = g = c = 0
    for line in lines:
        for ch in line.upper():
            if ch == "A": a += 1
            elif ch == "T": t += 1
            elif ch == "G": g += 1
            elif ch == "C": c += 1
            else:
                # ignore other characters (headers, Ns, gaps, whitespace)
                pass
    return {"A": a, "T": t, "G": g, "C": c}


# -------------------- CLI helpers --------------------

def _benchmark(seq: str, repeats: int = 5) -> None:
    targets = [
        ("Counter", lambda: counts_counter(seq)),
        ("str.count", lambda: counts_strcount(seq)),
        ("manual loop", lambda: counts_manual_loop(seq)),
    ]
    print(f"\nSequence length: {len(seq)} | Repeats per method: {repeats}\n")
    for name, fn in targets:
        t = timeit.timeit(fn, number=repeats)
        res = fn()
        print(f"{name:<12} time: {t:.6f}s   result: {res}")
    print()

def _profile(seq: str, method: str = "counter") -> None:
    mapping = {
        "counter": counts_counter,
        "strcount": counts_strcount,
        "manual": counts_manual_loop,
    }
    func = mapping.get(method.lower())
    if not func:
        print("Choose method from: counter | strcount | manual")
        return
    pr = cProfile.Profile()
    pr.enable()
    res = func(seq)
    pr.disable()
    print("Result:", res, "\n")
    ps = pstats.Stats(pr).strip_dirs().sort_stats(pstats.SortKey.CUMULATIVE)
    ps.print_stats(12)

def _read_file_stream(path: str) -> Iterable[str]:
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            yield line


# -------------------- CLI --------------------

def _menu() -> None:
    print("Examples you can paste:")
    print("  seq:  ATGCATGCATGC")
    print("  file: path to a .txt / .fasta (one or many lines)\n")

    while True:
        print("Choose an option:")
        print("  1) Benchmark counting methods on a sequence")
        print("  2) Profile one method on a sequence")
        print("  3) Stream-count from FILE (memory-efficient)")
        print("  0) Exit")
        choice = input("Enter 1/2/3/0: ").strip()

        if choice == "0":
            break
        elif choice == "1":
            seq = input("Enter DNA (A/T/G/C only recommended, others ignored in some paths): ").strip()
            try:
                # ensure validation path once (so times reflect pure counting thereafter)
                DNA(seq)
            except InvalidNucleotideError:
                # ok for benchmarking, but tell user; we’ll still test counts on uppercase
                print("Note: Invalid characters for strict DNA; benchmarking may still run.\n")
            repeats = input("Repeats [5]: ").strip()
            try:
                r = int(repeats) if repeats else 5
            except ValueError:
                r = 5
            _benchmark(seq, repeats=r)

        elif choice == "2":
            seq = input("Enter DNA (A/T/G/C only recommended): ").strip()
            method = input("Method (counter/strcount/manual) [counter]: ").strip() or "counter"
            _profile(seq, method)

        elif choice == "3":
            path = input("Enter file path: ").strip()
            t0 = time.time()
            try:
                res = stream_counts(_read_file_stream(path))
                t1 = time.time()
                print("\nStream counts:", res)
                print(f"Elapsed: {t1 - t0:.3f}s (memory-efficient)\n")
            except FileNotFoundError:
                print("File not found.\n")
        else:
            print("Please enter 1/2/3/0.\n")


if __name__ == "__main__":
    _menu()
