"""Independent stock support cannot come from co-movers or repeated aligner votes."""
import importlib.util
from pathlib import Path

import pysam


_spec = importlib.util.spec_from_file_location(
    "recount_junction_support", Path(__file__).resolve().parents[1] / "scripts/recount_junction_support.py")
audit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(audit)


def _bam(path, reads):
    with pysam.AlignmentFile(str(path), "wb", header={"SQ": [{"SN": "chrT", "LN": 1000}]}) as out:
        for name, start, cigar, flag in reads:
            rec = pysam.AlignedSegment(out.header)
            rec.query_name = name
            rec.reference_id = 0
            rec.reference_start = start
            rec.cigarstring = cigar
            rec.flag = flag
            rec.query_sequence = "A" * sum(n for op, n in rec.cigartuples if op in (0, 1, 4, 7, 8))
            out.write(rec)
    return path


def _target(name, start=100, end=200):
    return dict(read_id=name, chrom="chrT", intron_start=start, intron_end=end, note="review")


def test_co_movers_and_secondary_alignments_cannot_create_stock_support(tmp_path):
    stock = _bam(tmp_path / "stock.bam", [
        ("moving1", 80, "20M90N20M", 0),
        ("moving2", 80, "20M90N20M", 0),
        ("secondary", 80, "20M100N20M", 0x100),
        ("supplementary", 80, "20M100N20M", 0x800),
    ])
    rows = audit.recount_support([stock], [_target("moving1"), _target("moving2")])
    assert all(row["stock_other_qnames"] == 0 for row in rows)
    assert all(row["note"] == "review" for row in rows)


def test_unique_stock_qnames_across_arms_exclude_the_moving_read(tmp_path):
    records = [
        ("moving", 80, "20M100N20M", 0),
        ("peer", 80, "20M100N20M", 16),
        ("peer", 80, "20M100N20M", 0),  # duplicate primary must not add a vote
        ("short", 90, "10M100N20M", 0),
        ("adjacent_indel", 80, "20M1I100N20M", 0),
    ]
    paths = [_bam(tmp_path / f"arm{i}.bam", records) for i in range(2)]
    row = audit.recount_support(paths, [_target("moving")])[0]
    assert row["stock_support_qnames"] == 4
    assert row["stock_other_qnames"] == 3
    assert row["stock_other_anchored_qnames"] == 1
