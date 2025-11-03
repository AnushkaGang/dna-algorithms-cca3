# Assignment 1: DNA Fundamentals and Basic Tools (CCA3) — Report


**Name:** Anushka Gamgwar  
**Course:** TY B.Tech CSE (Panel-A)  
**PRN:** 1032233324 (roll no 46)
**Language Used:** Python  
**Structure Used:** Python Classes + Functions + CLI

---

## Purpose of this Assignment
This assignment implements foundational computational genomics skills using Python.

Core goals achieved:

- Represent DNA using Python data structures
- Validate nucleotides
- Perform nucleotide analysis
- Perform DNA transcription (coding + template strand)
- Generate reverse complement (IUPAC supported)
- Optimize algorithms and analyse performance
- Build a complete test suite

---

## Summary of What is Implemented

| Part | Details |
|------|---------|
| **Part A** | DNA data structures, cleaning, counting, codon splitting, merging |
| **Part B** | DNA → RNA transcription, reverse complement with IUPAC degenerates |
| **Part C** | Optimization, profiling, + full testing suite (unit tests) |

All final implementations include **user-input CLI menus** to demonstrate usage live.

---

## Performance Notes / Complexity Summary

| Operation | Time Complexity | Space Complexity |
|-----------|-----------------|------------------|
| Counting nucleotides (Counter, str.count, manual) | O(n) | O(1) |
| Complement / Reverse Complement | O(n) | O(1) |
| Cleaning sequences | O(n) | O(1) |
| Codon splitting | O(n) | O(1) |
| Transcription | O(n) | O(1) |
| Streaming counts from huge files | O(total length) | O(1) |

I compared 3 methods of counting nucleotides:

- `Counter()` — clean + fast
- `str.count()` — extremely fast (CPython C-level)
- manual loop — baseline reference

Result: all are O(n).  
`str.count()` is fastest in CPython.  
Streaming memory method allows processing large FASTA without RAM spike.

---

## Test Suite

All features are unit tested:

- empty sequences
- single base sequences
- invalid symbols
- random DNA
- random IUPAC
- orientation logic
- transcription logic

Test file: `tests/test_dna.py`  




---

## Tools Used

| Tool | Purpose |
|------|---------|
| Python + venv | Language Runtime |
| VS Code | Editor |
| GitHub | Version Control |
| unittest | Testing Framework |

---

## Conclusion

All requirements for Assignment 1 (CCA3) were completed:

- Code quality maintained
- Documented functions
- Time & space complexity explained
- Testing done with edge cases
- GitHub repository properly maintained



