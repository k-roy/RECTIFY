"""Conservative native-continuation correction for tail-supported terminal introns.

The two placements compare the same non-tail query bases.  Tail bases provide
end context only, never positional support.  This is a deliberately narrow
policy: exact native continuation, an unannotated terminal junction, complete
query/tail context, and a composition-conditioned fixed-position null bound.
It does not estimate the error rate of an aligner's genome-wide search.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
import json
import math
from typing import Callable, Iterable, Optional

from rectify.config import MIN_POLYA_LENGTH
from .overhang_informativeness import DEFAULT_ALPHA, min_self_match_period

TAG = "Zt"
_QUERY = frozenset((0, 1, 4, 7, 8))
_REF = frozenset((0, 2, 3, 7, 8))
_MATCH = frozenset((0, 7, 8))
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


@dataclass(frozen=True)
class TerminalTailResult:
    applied: bool
    reason: str
    evidence: dict = field(default_factory=dict)


def _rc(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def conditional_exact_bound(prefix: str) -> float:
    """Union bound for TWO fixed positions under base-count permutations ending in non-A.

    The maximal terminal A run is conditioned on, not searched over.  A query
    partition selected from several tail splits would need another multiplicity
    factor and is not supported by this policy.  Real serial dependence and an
    aligner's preceding genome-wide search are not calibrated by this bound.
    """
    counts = Counter(prefix)
    log_arrangements = math.lgamma(len(prefix) + 1) - sum(
        math.lgamma(count + 1) for count in counts.values())
    # Conditioning on the maximal terminal A run also fixes prefix[-1] != A.
    # Only arrangements with a non-A final base preserve that partition.
    non_a_fraction = (len(prefix) - counts.get("A", 0)) / len(prefix) if prefix else 0
    if not non_a_fraction:
        return 1.0
    return min(1.0, 2.0 * math.exp(-log_arrangements) / non_a_fraction)


def _merge(ops):
    merged = []
    for op, length in ops:
        if length <= 0:
            continue
        if merged and merged[-1][0] == op:
            merged[-1] = (op, merged[-1][1] + length)
        else:
            merged.append((op, length))
    return merged


def _query_placements(read):
    """Yield query placements without materializing potentially megabase N runs."""
    qpos, rpos = 0, read.reference_start
    for op, length in read.cigartuples or ():
        if op in _MATCH:
            for offset in range(length):
                yield qpos + offset, rpos + offset
        elif op in (1, 4):
            for offset in range(length):
                yield qpos + offset, None
        if op in _QUERY:
            qpos += length
        if op in _REF:
            rpos += length


def refine_terminal_tail(
    read,
    fetch: Callable[[str, int, int], str],
    *,
    rna_strand: Optional[str],
    protocol: Optional[str],
    annotated_junctions: Optional[Iterable[tuple]],
) -> TerminalTailResult:
    """Mutate one read only when the narrow native-plus-tail model dominates.

    ``fetch`` uses genomic, zero-based half-open coordinates.  The caller must
    establish RNA strand and a long-RNA protocol (``drs`` or ``ont-cdna``).
    ``annotated_junctions`` must be an available annotation set, with keys
    ``(chrom,start,end)`` or ``(chrom,start,end,strand)``.  Missing annotation,
    unsupported context and every refusal leave the record byte-for-byte alone.
    A successful edit keeps SEQ/QUAL and all preceding RNA-body placements;
    placement-relative '=' bases are decoded BEFORE changing the CIGAR.
    """
    def refuse(reason, **evidence):
        return TerminalTailResult(False, reason, evidence)

    if protocol not in ("drs", "ont-cdna") or rna_strand not in ("+", "-"):
        return refuse("unsupported_rna_context")
    if annotated_junctions is None:
        return refuse("annotation_unavailable")
    if (read.is_unmapped or read.is_secondary or read.is_supplementary
            or read.has_tag("SA") or not read.query_sequence):
        return refuse("incomplete_or_chimeric_record")
    if read.has_tag(TAG):
        return refuse("provenance_tag_occupied")
    ops = list(read.cigartuples or ())
    indices = [i for i, (op, _) in enumerate(ops) if op == 3]
    if not indices:
        return refuse("no_terminal_intron")
    plus = rna_strand == "+"
    ni = indices[-1] if plus else indices[0]
    ns = read.reference_start + sum(n for op, n in ops[:ni] if op in _REF)
    ne = ns + ops[ni][1]
    chrom = read.reference_name
    if ((chrom, ns, ne) in annotated_junctions
            or (chrom, ns, ne, "+") in annotated_junctions
            or (chrom, ns, ne, "-") in annotated_junctions):
        return refuse("annotated_junction", junction=[ns, ne])
    distal_ops = ops[ni + 1:] if plus else ops[:ni]
    body_ops = ops[:ni] if plus else ops[ni + 1:]
    if (not body_ops or (body_ops[-1] if plus else body_ops[0])[0] not in _MATCH
            or any(op not in (0, 1, 2, 4, 7, 8) for op, _ in distal_ops)):
        return refuse("unanchored_or_incomplete_terminal_block")
    # A soft clip can occur only at the RNA 3' terminus of this block.
    if any(op == 4 for op, _ in (distal_ops[:-1] if plus else distal_ops[1:])):
        return refuse("internal_softclip")
    q_boundary = sum(n for op, n in ops[:ni] if op in _QUERY)
    sequence = read.query_sequence.upper()
    pairs = None
    if "=" in sequence:
        decoded = list(sequence)
        pairs = list(_query_placements(read))
        try:
            for qp, rp in pairs:
                if qp is not None and rp is not None and decoded[qp] == "=":
                    base = fetch(chrom, rp, rp + 1).upper()
                    if len(base) != 1:
                        return refuse("reference_unavailable")
                    decoded[qp] = base
        except (ValueError, KeyError, IndexError):
            return refuse("reference_unavailable")
        sequence = "".join(decoded)
    if "=" in sequence:
        return refuse("unresolved_equals_sequence")
    terminal = sequence[q_boundary:] if plus else _rc(sequence[:q_boundary])
    prefix = terminal.rstrip("A")
    tail_length = len(terminal) - len(prefix)
    if tail_length < MIN_POLYA_LENGTH:
        return refuse("insufficient_terminal_tail", tail_length=tail_length)
    if not prefix or set(prefix) - set("ACGT"):
        return refuse("noninformative_prefix")
    bound = conditional_exact_bound(prefix)
    if bound > DEFAULT_ALPHA or min_self_match_period(prefix) is not None:
        return refuse("conditional_null_not_supported", prefix=prefix, null_bound=bound)
    n = len(prefix)
    start, end = (ns, ns + n) if plus else (ne - n, ne)
    if start < 0:
        return refuse("contig_edge")
    try:
        native = fetch(chrom, start, end).upper()
    except (ValueError, KeyError, IndexError):
        return refuse("reference_unavailable")
    if not plus:
        native = _rc(native)
    if len(native) != n:
        return refuse("contig_edge")
    if native != prefix:
        return refuse("native_not_exact", prefix=prefix, native=native)
    # Verify the body base adjacent to the disappearing N, not just an op edge.
    anchor_q, anchor_r = ((q_boundary - 1, ns - 1) if plus
                          else (q_boundary, ne))
    try:
        if sequence[anchor_q] != fetch(chrom, anchor_r, anchor_r + 1).upper():
            return refuse("body_anchor_mismatch")
    except (ValueError, KeyError, IndexError):
        return refuse("reference_unavailable")
    qlo, qhi = ((q_boundary, q_boundary + n) if plus
                else (q_boundary - n, q_boundary))
    # Actual base placements only: leading I/D cannot masquerade as an anchor.
    if pairs is None:
        pairs = list(_query_placements(read))
    current_matches = 0
    try:
        for qp, rp in pairs:
            if qp is not None and qlo <= qp < qhi and rp is not None:
                current_matches += sequence[qp] == fetch(chrom, rp, rp + 1).upper()
    except (ValueError, KeyError, IndexError):
        return refuse("reference_unavailable")
    if current_matches == n:
        return refuse("no_positional_gain", prefix=prefix)
    # Native is gap-free and exact, the best possible placement under any
    # positive mismatch/gap cost.  The observed placement loses query matches.
    evidence = dict(v=1, model="native_plus_terminal_tail", protocol=protocol,
                    strand=rna_strand, junction=[ns, ne], old_start=read.reference_start,
                    old_cigar=read.cigarstring, prefix=prefix, prefix_query=[qlo, qhi],
                    tail_length=tail_length, current_matches=current_matches,
                    native_matches=n, native_interval=[start, end],
                    null_bound=bound, fixed_positions=2, tail_splits_searched=1)
    new_ops = (_merge(body_ops + [(0, n), (4, tail_length)]) if plus
               else _merge([(4, tail_length), (0, n)] + body_ops))
    if sum(length for op, length in new_ops if op in _QUERY) != len(sequence):
        return refuse("query_consumption_mismatch")
    qualities = read.query_qualities
    read.query_sequence = sequence
    if qualities is not None:
        read.query_qualities = qualities
    read.cigartuples = new_ops
    if not plus:
        read.reference_start = start
    # These tags describe the old alignment, not the unchanged molecule.
    for tag in ("MD", "NM", "AS", "ms", "cs", "de", "dv", "SA"):
        if read.has_tag(tag):
            read.set_tag(tag, None)
    read.set_tag(TAG, json.dumps(evidence, separators=(",", ":")), value_type="Z")
    return TerminalTailResult(True, "exact_native_dominance", evidence)


def modern_cdna_rna_strand(read) -> Optional[str]:
    """Recognize the stage-1 RNA-sense contract, not legacy XO frame guesses.

    XN:i:1 is emitted with XT:i:{1,2} and XR:Z:source-read-membership by
    ``cdna.io.write_consensus_fastq``.  XN alone is insufficient to activate a
    cDNA policy; typed companions prevent accidental cross-protocol activation.
    """
    for tag, types in (("XN", "cCsSiI"), ("XT", "cCsSiI"), ("XR", "Z")):
        if not read.has_tag(tag):
            return None
        _, value_type = read.get_tag(tag, with_value_type=True)
        if value_type not in types:
            return None
    if (read.get_tag("XN") != 1 or read.get_tag("XT") not in (1, 2)
            or not read.get_tag("XR").strip()):
        return None
    return "-" if read.is_reverse else "+"


def load_annotation_junctions(annotation_path: str) -> set:
    """Use the same GFF/GTF annotation interpretation as aligner hints."""
    from rectify.utils.junction_bed import parse_gff_introns, parse_gtf_introns
    name = str(annotation_path).lower()
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith((".gff", ".gff3")):
        rows = parse_gff_introns(str(annotation_path))
    elif name.endswith(".gtf"):
        rows = parse_gtf_introns(str(annotation_path))
    else:
        raise ValueError(f"Unsupported junction annotation: {annotation_path}")
    return {(chrom, start, end) for chrom, start, end, _, _ in rows}


def run_terminal_tail_postpass(input_bam, genome_path, output_bam,
                               annotated_junctions, context=None):
    """Stream one alignment arm; only modern cDNA records can change.

    Returns ``(selected_bam_path, counters)``.  An unchanged arm retains its
    original file. Changed output is sorted because minus-strand edits change
    POS; the original BAM is retained. No additional aligner arm is introduced.
    """
    from pathlib import Path
    import tempfile
    import pysam
    from contextlib import ExitStack
    from ..align.reference import open_alignment_reference

    context = context if context is not None else terminal_tail_context()
    counters = Counter()
    output = Path(output_bam)
    # Detect a real edit before writing any full-size output. DRS and legacy
    # cDNA arms without the modern frame contract incur only a streaming scan.
    needs_output = False
    with pysam.AlignmentFile(str(input_bam), "rb") as source:
        reference = None
        reference_stack = ExitStack()
        try:
            for read in source:
                counters["records"] += 1
                strand = modern_cdna_rna_strand(read)
                if strand is None:
                    counters["unsupported_rna_context"] += 1
                    continue
                if not any(op == 3 for op, _ in read.cigartuples or ()):
                    counters["no_terminal_intron"] += 1
                    continue
                if callable(annotated_junctions):
                    annotated_junctions = annotated_junctions()
                if reference is None:
                    reference = reference_stack.enter_context(open_alignment_reference(genome_path, output.parent))
                result = refine_terminal_tail(read, reference.fetch, rna_strand=strand,
                    protocol="ont-cdna", annotated_junctions=annotated_junctions)
                counters[result.reason] += 1
                if result.applied:
                    needs_output = True
                    break
        finally:
            reference_stack.close()
    if not needs_output:
        # A new no-op decision invalidates an earlier successful receipt even
        # when the input itself is unchanged (e.g. a changed annotation).
        receipt = Path(str(output) + ".terminal_tail.json")
        if receipt.exists():
            receipt.unlink()
        return str(input_bam), dict(counters)
    counters.clear()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".terminal-tail-", dir=output.parent) as tmp:
        unsorted = str(Path(tmp) / "edits.bam")
        with open_alignment_reference(genome_path, output.parent) as ref, \
                pysam.AlignmentFile(str(input_bam), "rb") as source, \
                pysam.AlignmentFile(unsorted, "wb", header=source.header) as sink:
            for read in source:
                counters["records"] += 1
                strand = modern_cdna_rna_strand(read)
                if strand is None:
                    counters["unsupported_rna_context"] += 1
                else:
                    result = refine_terminal_tail(read, ref.fetch, rna_strand=strand,
                        protocol="ont-cdna", annotated_junctions=annotated_junctions)
                    counters[result.reason] += 1
                    counters["applied"] += int(result.applied)
                sink.write(read)
        sorted_bam = str(Path(tmp) / "sorted.bam")
        pysam.sort("-@", "1", "-o", sorted_bam, unsorted)
        pysam.index(sorted_bam)
        Path(sorted_bam).replace(output)
        Path(sorted_bam + ".bai").replace(Path(str(output) + ".bai"))
        manifest = dict(version=1, context=context, input=_file_identity(input_bam),
                        output=_file_identity(output), counters=dict(counters))
        manifest_temp = Path(tmp) / "manifest.json"
        manifest_temp.write_text(json.dumps(manifest, indent=2) + "\n")
        manifest_temp.replace(Path(str(output) + ".terminal_tail.json"))
    return str(output), dict(counters)


def _file_identity(path):
    from pathlib import Path
    p = Path(path)
    stat = p.stat()
    return dict(path=str(p.resolve()), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def terminal_tail_output_for(input_bam, output_bam, context=None):
    """Return a completed post-pass only when its exact input/output still match.

    Run-all discovers arms from disk after ``run_align`` returns. This receipt
    keeps that path consistent with the in-process substitution without reusing
    stale post-pass output after a source BAM has been replaced.
    """
    from pathlib import Path
    output = Path(output_bam)
    context = context if context is not None else terminal_tail_context()
    try:
        manifest = json.loads(Path(str(output) + ".terminal_tail.json").read_text())
        if (manifest.get("version") == 1
                and manifest.get("context") == context
                and manifest.get("input") == _file_identity(input_bam)
                and manifest.get("output") == _file_identity(output)
                and Path(str(output) + ".bai").is_file()):
            return output
    except (OSError, ValueError, TypeError):
        pass
    return None


def terminal_tail_context(genome_path=None, annotation_path=None):
    """Fingerprint implementation, fixed policy and optional run references."""
    import hashlib
    from pathlib import Path
    return dict(implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                protocol="modern_ont_cdna_XN1_XT12_XR", min_tail=MIN_POLYA_LENGTH,
                alpha=DEFAULT_ALPHA,
                genome=_file_identity(genome_path) if genome_path else None,
                annotation=_file_identity(annotation_path) if annotation_path else None)


def write_terminal_selection_receipt(output_bam, selected_arms, context):
    """Bind the cached consensus to the exact arm placements and RNA policy."""
    from pathlib import Path
    output = Path(output_bam)
    arms = {}
    for name, path in selected_arms.items():
        record = {"file": _file_identity(path)}
        receipt = Path(str(path) + ".terminal_tail.json")
        if receipt.exists():
            record["terminal_receipt"] = _file_identity(receipt)
            record["terminal_input"] = json.loads(receipt.read_text())["input"]
        arms[name] = record
    payload = dict(version=1, context=context, output=_file_identity(output), arms=arms)
    marker = Path(str(output) + ".terminal_selection.json")
    temporary = Path(str(marker) + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    temporary.replace(marker)


def terminal_selection_matches(output_bam, context):
    """A raw-aligner trust override never authorizes stale post-pass selection."""
    from pathlib import Path
    marker = Path(str(output_bam) + ".terminal_selection.json")
    if not marker.exists():
        return context is None
    try:
        payload = json.loads(marker.read_text())
        if (payload.get("version") != 1 or payload.get("context") != context
                or payload.get("output") != _file_identity(output_bam)):
            return False
        for record in payload["arms"].values():
            for key in ("file", "terminal_receipt", "terminal_input"):
                expected = record.get(key)
                if expected is not None and expected != _file_identity(expected["path"]):
                    return False
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False
