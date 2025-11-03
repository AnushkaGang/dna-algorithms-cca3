"""
tests/test_dna.py — Part C, Q7: Comprehensive Testing Suite

Run (either way):
    $ python tests/test_dna.py
    # or
    $ python -m unittest tests/test_dna.py

Includes:
- Unit tests for dna_core, dna_analysis, dna_strings, dna_transcription, dna_reverse_complement
- Edge cases and random generators
"""

import random
import string
import unittest

# Import modules under test (assumes they are in parent folder)
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from dna_core import DNA, InvalidNucleotideError
from dna_analysis import count_nucleotides, nucleotide_frequencies, analysis_report, compare_composition
from dna_strings import clean_sequence, split_codons, merge_fragments, to_upper, to_lower
from dna_transcription import transcribe, dna_to_rna
from dna_reverse_complement import complement_iupac, reverse_complement_iupac, rc_with_orientations


# ------------- Random data generators -------------

VALID = "ATGC"
IUPAC = "ATGCRYSWKMBDHVN"

def rand_dna(n: int) -> str:
    """Random DNA (A/T/G/C) of length n."""
    return "".join(random.choice(VALID) for _ in range(n))

def rand_iupac(n: int) -> str:
    """Random IUPAC DNA of length n (for degenerate RC tests)."""
    return "".join(random.choice(IUPAC) for _ in range(n))


# ------------------------- Tests -------------------------

class TestDNAClass(unittest.TestCase):
    def test_empty(self):
        d = DNA("")
        self.assertEqual(len(d), 0)
        self.assertEqual(d.counts(), {"A":0,"T":0,"G":0,"C":0})
        self.assertEqual(d.gc_content(), 0.0)
        self.assertEqual(d.at_content(), 0.0)

    def test_valid_counts(self):
        d = DNA("ATGCA")
        self.assertEqual(d.counts(), {"A":2,"T":1,"G":1,"C":1})

    def test_invalid_chars(self):
        with self.assertRaises(InvalidNucleotideError):
            DNA("ATN-X")

class TestAnalysis(unittest.TestCase):
    def test_frequencies_sum(self):
        s = "ATGCA"
        f = nucleotide_frequencies(s)
        self.assertAlmostEqual(sum(f.values()), 100.0, places=7)

    def test_compare_composition_symmetry(self):
        s1, s2 = "ATGCA", "AAAAA"
        diffs = compare_composition(s1, s2)
        diffs2 = compare_composition(s2, s1)
        self.assertAlmostEqual(diffs["L1_distance%"], diffs2["L1_distance%"], places=7)

class TestStrings(unittest.TestCase):
    def test_clean_sequence(self):
        messy = "atg cAx\n---tgtN"
        self.assertEqual(clean_sequence(messy), "ATGCATGT")

    def test_split_codons(self):
        s = "ATGCATG"
        codons, leftover = split_codons(s, frame=0, drop_incomplete=True)
        self.assertEqual(codons, ["ATG","CAT"])
        self.assertEqual(leftover, "G")

    def test_merge(self):
        frags = ["ATG", " caa-tt", "g"]
        self.assertEqual(merge_fragments(frags), "ATGCAATTG")

class TestTranscription(unittest.TestCase):
    def test_coding(self):
        # coding: mRNA = coding with T->U
        self.assertEqual(transcribe("ATGGACT", strand="coding"), "AUGGACU")

    def test_template_5to3(self):
        # template(5'->3') reverse-complement -> RNA(T->U)
        self.assertEqual(transcribe("TACCTGA", strand="template", template_orientation="5to3"), "AUGGACU")

    def test_dna_to_rna(self):
        self.assertEqual(dna_to_rna("TTT"), "UUU")

class TestReverseComplement(unittest.TestCase):
    def test_iupac_complement_pairs(self):
        self.assertEqual(complement_iupac("ATGC"), "TACG")
        self.assertEqual(complement_iupac("RYKMBVDHNWS"), "YRMKVBHDNW S".replace(" ", ""))  # spaces to keep line short

    def test_revcomp_involution(self):
        s = rand_iupac(100)
        rc = reverse_complement_iupac(s)
        # RC(RC(s)) == s
        rc2 = reverse_complement_iupac(rc)
        self.assertEqual(rc2, s)

    def test_orientation_handling(self):
        s = "ATGC"
        # If input given as 3'->5' (CGTA), revcomp to 5'->3' should be GCAT
        res = rc_with_orientations("CGTA", op="revcomp", input_orientation="3to5", output_orientation="5to3")
        self.assertEqual(res, "GCAT")

class TestRandomProperties(unittest.TestCase):
    def test_counts_sum_length(self):
        for _ in range(5):
            s = rand_dna(1000)
            c = count_nucleotides(s)
            self.assertEqual(sum(c.values()), len(s))

    def test_freqs_consistent(self):
        for _ in range(3):
            s = rand_dna(500)
            f = nucleotide_frequencies(s)
            # Sum ~ 100 within floating rounding tolerance
            self.assertTrue(abs(sum(f.values()) - 100.0) < 1e-9)

# --------------- CLI ---------------

if __name__ == "__main__":
    # allow quick run without unittest discovery
    unittest.main(verbosity=2)
