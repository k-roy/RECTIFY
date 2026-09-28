"""An accepted SSP indel must not remove an mRNA base or retain bridge sequence."""
import random

import pysam
import pytest

from rectify.core.cdna._constants import ANCHOR_FWD, SSP_FWD
from rectify.core.cdna.consensus import pretrim_consensus
from rectify.core.cdna.io import write_stage1_fastq
from rectify.core.cdna.read_info import extract_read_info, revcomp


def primer(kind):
    if kind == "exact":
        return SSP_FWD
    if kind == "substitution":
        return SSP_FWD[:-1] + "C"
    operation, position = kind.split("_")
    i = int(position)
    return SSP_FWD[:i] + ("A" + SSP_FWD[i:] if operation == "insertion" else SSP_FWD[i + 1:])


@pytest.mark.parametrize("kind", ["exact", "substitution", "deletion_8", "deletion_10",
                                  "insertion_8", "insertion_12"])
@pytest.mark.parametrize("frame", ["fwd", "rev"])
def test_fuzzy_primer_span_preserves_known_body_in_emitted_fastq(tmp_path, kind, frame):
    rng = random.Random(812)
    body = "CT" + "".join(rng.choice("ACGT") for _ in range(394)) + "CGTC"
    # Unique bridge boundary; the original G-ending UMI forms a tied four-G
    # tract for the terminal substitution and is tested as an explicit refusal.
    umi = "ACGTTGCAACGTTGCAACGTTGCAACC"
    prefix = primer(kind) + umi + "GGG"
    tail = "A" * 30 + ANCHOR_FWD + "CTGCTCGTGC"
    seq = prefix + body + tail
    if frame == "rev":
        seq = revcomp(seq)
    h = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chrT", "LN": 2000}]})
    read = pysam.AlignedSegment(h)
    read.query_name = f"{kind}_{frame}"
    read.reference_id = 0
    read.reference_start = 500
    read.mapping_quality = 60
    read.flag = 16 if frame == "rev" else 0
    read.query_sequence = seq
    read.query_qualities = [30] * len(seq)
    read.cigarstring = (f"{len(prefix)}S400M{len(tail)}S" if frame == "fwd"
                        else f"{len(tail)}S400M{len(prefix)}S")
    info = extract_read_info(read)
    assert info.read_type == 1 and info.umi == umi
    pretrim = pretrim_consensus(seq, frame, info.read_type)
    assert pretrim.trim_5p == len(prefix)
    assert pretrim.seq == (body if frame == "fwd" else revcomp(body))
    bam, fastq = tmp_path / "input.bam", tmp_path / "emitted.fq"
    with pysam.AlignmentFile(str(bam), "wb", header=h) as handle:
        handle.write(read)
    write_stage1_fastq(bam, fastq, [[info]], {0: info.umi}, {0: info.xf_tier}, {0: 30},
                       reference=None)
    lines = fastq.read_text().splitlines()
    assert lines[1] == body
    assert len(lines[3]) == len(body)
    tags = dict(token.split(":", 2)[::2] for token in lines[0].split()[1:] if token.count(":") >= 2)
    assert tags["XN"] == "1"
    assert int(tags["XQ"]) == len(prefix)
