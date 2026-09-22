#!/usr/bin/env python3
"""Regression tests for collinear_block_merge (C3 acceptance gate)."""

import unittest

from collinear_block_merge import (
    can_merge_segments,
    merge_alignment_segments,
    normalize_blast_hsp,
)


def _seg(qs, qe, ts, te, strand="+", nmatch=1000, alen=1000):
    return {
        "qstart": qs,
        "qend": qe,
        "tstart": ts,
        "tend": te,
        "strand": strand,
        "nmatch": nmatch,
        "alen": alen,
    }


class TestCollinearBlockMerge(unittest.TestCase):
    max_gap = 200
    overlap_tol = 50

    def test_a_forward_translocation_must_split(self):
        """Same-strand query-adjacent + ref jump >> max_gap → 2 blocks."""
        records = [
            _seg(0, 5000, 0, 5000),
            _seg(5100, 10100, 56000, 61100),  # qgap=100, rgap=51000
        ]
        self.assertFalse(
            can_merge_segments(records[0], records[1], self.max_gap, self.overlap_tol)
        )
        blocks = merge_alignment_segments(
            records, max_gap=self.max_gap, overlap_tol=self.overlap_tol
        )
        self.assertEqual(len(blocks), 2)

    def test_b_collinear_adjacent_stays_one_block(self):
        records = [
            _seg(0, 5000, 0, 5000),
            _seg(5100, 10100, 5100, 10100),  # qgap=100, rgap=100
        ]
        self.assertTrue(
            can_merge_segments(records[0], records[1], self.max_gap, self.overlap_tol)
        )
        blocks = merge_alignment_segments(
            records, max_gap=self.max_gap, overlap_tol=self.overlap_tol
        )
        self.assertEqual(len(blocks), 1)
        self.assertEqual(blocks[0]["n_segments"], 2)

    def test_c_opposite_strand_inversion_must_split(self):
        records = [
            _seg(0, 5000, 0, 5000, strand="+"),
            _seg(5100, 10100, 5100, 10100, strand="-"),
        ]
        self.assertFalse(
            can_merge_segments(records[0], records[1], self.max_gap, self.overlap_tol)
        )
        blocks = merge_alignment_segments(
            records, max_gap=self.max_gap, overlap_tol=self.overlap_tol
        )
        self.assertEqual(len(blocks), 2)

    def test_d_micro_overlap_does_not_split(self):
        """Tiny overlap within overlap_tol → still one block (C1)."""
        records = [
            _seg(0, 5000, 0, 5000),
            _seg(4990, 9990, 4990, 9990),  # qgap=-10, rgap=-10
        ]
        self.assertTrue(
            can_merge_segments(records[0], records[1], self.max_gap, self.overlap_tol)
        )
        blocks = merge_alignment_segments(
            records, max_gap=self.max_gap, overlap_tol=self.overlap_tol
        )
        self.assertEqual(len(blocks), 1)

    def test_blast_minus_strand_normalization(self):
        norm = normalize_blast_hsp(100, 600, 800, 300, "minus")
        self.assertEqual(norm["strand"], "-")
        self.assertLess(norm["tstart"], norm["tend"])
        self.assertEqual(norm["tstart"], 300)
        self.assertEqual(norm["tend"], 800)


if __name__ == "__main__":
    unittest.main()
