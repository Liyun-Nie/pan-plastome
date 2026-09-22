"""
Dual-axis collinear block merging for PAF / BLAST alignment segments.

Aligned with mitoCart identify_hsp_blocks() (S2.10.3) semantics, with
overlap tolerance (C1) for minimap2 -c chained segments.

Merge iff same_strand AND both axes satisfy:
    -overlap_tol <= gap <= max_gap
"""

from __future__ import annotations


def ref_gap(cur: dict, rec: dict) -> int:
    """Signed gap on reference in strand-forward direction."""
    if cur["strand"] == "+":
        return rec["tstart"] - cur["tend"]
    return cur["tstart"] - rec["tend"]


def axis_gap_adjacent(gap: int, max_gap: int, overlap_tol: int) -> bool:
    return -overlap_tol <= gap <= max_gap


def can_merge_segments(
    cur: dict,
    rec: dict,
    max_gap: int = 200,
    overlap_tol: int = 50,
) -> bool:
    if rec["strand"] != cur["strand"]:
        return False
    qgap = rec["qstart"] - cur["qend"]
    rgap = ref_gap(cur, rec)
    return axis_gap_adjacent(qgap, max_gap, overlap_tol) and axis_gap_adjacent(
        rgap, max_gap, overlap_tol
    )


def merge_alignment_segments(
    records: list[dict],
    max_gap: int = 200,
    overlap_tol: int = 50,
) -> list[dict]:
    """Merge alignment segments into synteny blocks (dual-axis adjacency)."""
    if not records:
        return []

    blocks = []
    cur = _segment_to_block(records[0], block_id=0)

    for rec in records[1:]:
        if can_merge_segments(cur, rec, max_gap=max_gap, overlap_tol=overlap_tol):
            cur["qend"] = max(cur["qend"], rec["qend"])
            if cur["strand"] == "+":
                cur["tend"] = max(cur["tend"], rec["tend"])
            else:
                cur["tstart"] = min(cur["tstart"], rec["tstart"])
            cur["nmatch"] += rec.get("nmatch", 0)
            cur["alen"] += rec.get("alen", 0)
            cur["n_segments"] += 1
        else:
            blocks.append(_finalize_block(cur))
            cur = _segment_to_block(rec, block_id=len(blocks))

    blocks.append(_finalize_block(cur))
    for i, b in enumerate(blocks):
        b["block_id"] = i
    return blocks


def _segment_to_block(rec: dict, block_id: int) -> dict:
    return {
        "block_id": block_id,
        "qstart": rec["qstart"],
        "qend": rec["qend"],
        "tstart": rec["tstart"],
        "tend": rec["tend"],
        "strand": rec["strand"],
        "nmatch": rec.get("nmatch", 0),
        "alen": rec.get("alen", 0),
        "n_segments": 1,
    }


def _finalize_block(cur: dict) -> dict:
    span_q = max(1, cur["qend"] - cur["qstart"])
    span_t = max(1, abs(cur["tend"] - cur["tstart"]))
    cur["identity"] = round(100.0 * cur["nmatch"] / span_q, 2)
    cur["block_len"] = span_q
    return cur


def normalize_blast_hsp(
    qstart: int,
    qend: int,
    sstart: int,
    send: int,
    sstrand: str,
) -> dict:
    """
    Normalize BLAST outfmt6 coordinates for collinear_block_merge (C2).

    Ensures tstart < tend and strand in {'+', '-'} regardless of BLAST
    sstart/send ordering on minus strand.
    """
    strand = "+" if str(sstrand).strip() in ("+", "plus", "1") else "-"
    tstart, tend = sorted((int(sstart), int(send)))
    return {
        "qstart": int(qstart),
        "qend": int(qend),
        "tstart": tstart,
        "tend": tend,
        "strand": strand,
    }
