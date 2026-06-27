#!/usr/bin/env python3
"""Tier-2 read-simulator wrapper — pbsim3 over a transcript set (component 1).

WHY pbsim3 (recorded justification — see dev/SIMULATION_BENCHMARK_SPEC.md):
the framing metric is EXACT INDEL-POSITION CONCORDANCE, which needs a per-read
ground-truth alignment. Of the three candidate simulators (pbsim3 / badread /
nanosim), **only pbsim3 emits a per-read read<->template alignment (a MAF)**.
badread and nanosim report transcript-of-origin + identity but no per-read edit
script, so they can validate junction truth but NOT exact-indel truth. pbsim3's
MAF composed with the known transcript<->genome exon structure
(``TranscriptModel.transcript_pos_to_genome``) yields the exact read<->genome
alignment — indels AND junctions — that the framing metric requires.

DRS vs cDNA error models: pbsim3 ships ``--method errhmm`` with packaged ONT/
RNA error HMMs (``ERRHMM-ONT.model`` / ``ERRHMM-ONT-HQ.model``); for direct-RNA
(DRS) use the ONT model with higher error, for PCR-cDNA the HQ model. The
``error_model`` arg selects which packaged HMM to pass.

This wrapper:
1. writes the spliced transcript FASTA (one record per ``TranscriptModel``),
2. invokes ``pbsim --strategy templ --method errhmm`` (transcript mode),
3. parses each ``*.maf`` (read<->transcript), and
4. projects every read through ``transcript_pos_to_genome`` to a ``ReadTruth``
   with an exact genome CIGAR, junctions, and per-position indel truth.

Run on Sherlock (env ``aligner_bench`` has pbsim3); never relay BAMs via the M1.

Author: Kevin R. Roy
"""
from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import sys
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
from rectify.core.benchmark.truth_schema import (  # noqa: E402
    ReadTruth, IndelTruth, IndelKind, SplitTag, write_truth_table, homopolymer_run,
)
from scripts.benchmark.sim.transcript_model import TranscriptModel, revcomp  # noqa: E402


# ---------------------------------------------------------------------------
# MAF parsing
# ---------------------------------------------------------------------------
@dataclass
class MafRecord:
    ref_name: str       # template (transcript) name
    ref_start: int      # 0-based start in template
    ref_aln: str        # gapped template sequence
    read_name: str
    read_start: int
    read_aln: str       # gapped read sequence
    read_strand: str    # '+' or '-'


def _maybe_gzip_open(path: str):
    """Open a possibly-gzipped text file (pbsim3 emits ``.maf.gz`` / ``.fq.gz``)."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_maf(path: str) -> List[MafRecord]:
    """Parse a pbsim3 MAF. Each alignment block has two 's' lines: template then
    read. Coordinates are 0-based; a '-' strand 's' line means start is measured
    from the other end (pbsim reports srcSize to allow conversion)."""
    records: List[MafRecord] = []
    block: List[Tuple[str, int, int, str, int, str]] = []
    with _maybe_gzip_open(path) as fh:
        for line in fh:
            if line.startswith("a"):
                block = []
            elif line.startswith("s"):
                parts = line.split()
                # s name start alnSize strand srcSize text
                name, start, aln_size, strand, src_size, text = (
                    parts[1], int(parts[2]), int(parts[3]), parts[4],
                    int(parts[5]), parts[6])
                block.append((name, start, aln_size, strand, src_size, text))
                if len(block) == 2:
                    (rn, rs, _ra, rstr, rss, rtext) = block[0]
                    (qn, qs, _qa, qstr, qss, qtext) = block[1]
                    records.append(MafRecord(
                        ref_name=rn, ref_start=rs, ref_aln=rtext,
                        read_name=qn, read_start=qs, read_aln=qtext,
                        read_strand=qstr))
                    block = []
    return records


# ---------------------------------------------------------------------------
# Projection: read<->transcript MAF  o  transcript<->genome  ->  ReadTruth
# ---------------------------------------------------------------------------
def project_maf_record(rec: MafRecord, model: TranscriptModel,
                       stratum: str = "TRANSCRIPTOME",
                       annotated_pairs=None, annotated_donors=None,
                       annotated_acceptors=None) -> Optional[ReadTruth]:
    """Compose a MAF read<->transcript alignment with the transcript<->genome map
    to an exact read<->genome ``ReadTruth``. Indels in the MAF (gaps) become
    ``IndelTruth`` at their genome positions; introns spanned by the read become
    junction truth via the model.

    NOTE: pbsim's template is the spliced transcript (5'->3'). A read base aligned
    to template offset ``t`` maps to genome position ``transcript_pos_to_genome(t)``.
    Deletions (gap in read) consume template; insertions (gap in template) do not.
    Adjacent template positions whose genome positions are NOT contiguous indicate
    an intron the read spans -> emit an N (junction).
    """
    ref_aln = rec.ref_aln
    read_aln = rec.read_aln
    tpos = rec.ref_start  # offset into spliced transcript
    # Walk the alignment columns; accumulate genome-projected indel events.
    indels: List[IndelTruth] = []
    # Track template positions consumed so we can detect intron jumps and the
    # genome span. For the deliverable's exact-indel truth we record DEL/INS at
    # genome coords; junctions come from the model (the read spans whole introns).
    del_run = 0
    del_tstart: Optional[int] = None
    for rc, qc in zip(ref_aln, read_aln):
        if rc != "-" and qc != "-":
            if del_run:
                gpos = model.transcript_pos_to_genome(del_tstart)
                rs, re, base = homopolymer_run(model.genome_seq, gpos)
                ctx = "HP" if (re - rs) >= 2 else "UNIQUE"
                indels.append(IndelTruth(
                    pos=min(gpos, model.transcript_pos_to_genome(del_tstart + del_run - 1)),
                    length=del_run, kind=IndelKind.DEL,
                    eq_start=rs if ctx == "HP" else gpos,
                    eq_end=re if ctx == "HP" else gpos + del_run,
                    context=ctx, run_unit=base if ctx == "HP" else "",
                    run_copies=(re - rs) if ctx == "HP" else 0))
                del_run = 0
                del_tstart = None
            tpos += 1
        elif qc == "-" and rc != "-":   # deletion in read (template base missing)
            if del_run == 0:
                del_tstart = tpos
            del_run += 1
            tpos += 1
        elif rc == "-" and qc != "-":   # insertion in read
            gpos = model.transcript_pos_to_genome(max(tpos - 1, 0))
            indels.append(IndelTruth(
                pos=gpos, length=1, kind=IndelKind.INS,
                eq_start=gpos, eq_end=gpos + 1, context="UNIQUE"))
    if del_run:
        gpos = model.transcript_pos_to_genome(del_tstart)
        rs, re, base = homopolymer_run(model.genome_seq, gpos)
        ctx = "HP" if (re - rs) >= 2 else "UNIQUE"
        indels.append(IndelTruth(
            pos=gpos, length=del_run, kind=IndelKind.DEL,
            eq_start=rs if ctx == "HP" else gpos,
            eq_end=re if ctx == "HP" else gpos + del_run,
            context=ctx, run_unit=base if ctx == "HP" else "",
            run_copies=(re - rs) if ctx == "HP" else 0))

    junctions = model.junction_truths(annotated_pairs, annotated_donors,
                                      annotated_acceptors)
    return ReadTruth(
        read_id=rec.read_name, true_locus=model.name, true_transcript=model.name,
        chrom=model.chrom, strand=model.strand,
        genome_start=model.genome_start, genome_end=model.genome_end,
        junctions=junctions, indels=indels,
        true_cpa=(model.genome_end - 1) if model.strand == "+" else model.genome_start,
        stratum=stratum, split=SplitTag.TEST, coverage=None)


# ---------------------------------------------------------------------------
# pbsim3 invocation
# ---------------------------------------------------------------------------
def write_transcript_fasta(models: List[TranscriptModel], path: str) -> None:
    with open(path, "w") as fh:
        for m in models:
            fh.write(f">{m.name}\n{m.spliced_transcript()}\n")


def run_pbsim3(transcript_fa: str, out_prefix: str, errhmm_model: str,
               depth: int = 20, pbsim_bin: str = "pbsim",
               seed: Optional[int] = None) -> List[str]:
    """Invoke pbsim3 in transcript (templ) mode with an errhmm model. Returns the
    list of generated MAF paths. ``errhmm_model`` is the packaged model path
    (e.g. .../data/ERRHMM-ONT.model for DRS, ERRHMM-ONT-HQ.model for cDNA).

    NB: the ``errhmm`` method derives accuracy from the model itself; do NOT pass
    ``--accuracy-mean`` (that flag belongs to ``--method qshmm``)."""
    cmd = [pbsim_bin, "--strategy", "templ", "--method", "errhmm",
           "--errhmm", errhmm_model, "--template", transcript_fa,
           "--depth", str(depth), "--prefix", out_prefix]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"pbsim3 failed: {res.stderr[:800]}")
    out_dir = os.path.dirname(out_prefix) or "."
    base = os.path.basename(out_prefix)
    # pbsim3 builds vary: some emit ``<prefix>.maf`` / ``<prefix>_NNNN.maf``,
    # others gzip to ``<prefix>.maf.gz``. Match both (and never the FASTQ).
    mafs = sorted(os.path.join(out_dir, f) for f in os.listdir(out_dir)
                  if f.startswith(base) and (f.endswith(".maf") or f.endswith(".maf.gz")))
    if not mafs:
        raise RuntimeError(
            f"pbsim3 produced no MAF for prefix '{base}' in {out_dir}; "
            f"found: {sorted(os.listdir(out_dir))[:20]}. (pbsim ran rc=0 — check "
            f"the output naming/compression for this pbsim3 build.)")
    return mafs


def simulate_and_propagate(models: List[TranscriptModel], out_dir: str,
                           errhmm_model: str, depth: int = 20,
                           pbsim_bin: str = "pbsim", seed: Optional[int] = None,
                           annotated_pairs=None, annotated_donors=None,
                           annotated_acceptors=None, stratum: str = "TRANSCRIPTOME"
                           ) -> Dict[str, str]:
    """End-to-end: write transcript FASTA -> pbsim3 -> parse MAFs -> propagate
    truth -> write ``truth.tsv`` and concatenated reads FASTQ. Returns paths.

    Pass ``annotated_*`` sets (from ``gff_panel.annotated_sets``) to classify each
    read's junctions ANNOTATED/NIC/NNC; omit them and every junction is NNC (the
    synthetic-gene default used by ``live_roundtrip``)."""
    os.makedirs(out_dir, exist_ok=True)
    tfa = os.path.join(out_dir, "transcripts.fa")
    write_transcript_fasta(models, tfa)
    prefix = os.path.join(out_dir, "sim")
    mafs = run_pbsim3(tfa, prefix, errhmm_model, depth=depth,
                      pbsim_bin=pbsim_bin, seed=seed)
    by_name = {m.name: m for m in models}
    truth: List[ReadTruth] = []
    for maf in mafs:
        for rec in parse_maf(maf):
            model = by_name.get(rec.ref_name)
            if model is None:
                continue
            rt = project_maf_record(rec, model, stratum=stratum,
                                    annotated_pairs=annotated_pairs,
                                    annotated_donors=annotated_donors,
                                    annotated_acceptors=annotated_acceptors)
            if rt is not None:
                truth.append(rt)
    truth_tsv = os.path.join(out_dir, "truth.tsv")
    write_truth_table(truth, truth_tsv)
    # concat the fastqs pbsim emits (naming/compression varies by build:
    # ``sim.fastq`` / ``sim_NNNN.fastq`` / ``sim.fq`` / ``sim.fq.gz`` ...)
    base = os.path.basename(prefix)
    fq_exts = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    reads_fq = os.path.join(out_dir, "reads.fastq")
    n_fq = 0
    with open(reads_fq, "w") as out:
        for f in sorted(os.listdir(out_dir)):
            if f.startswith(base) and f.endswith(fq_exts):
                n_fq += 1
                with _maybe_gzip_open(os.path.join(out_dir, f)) as g:
                    out.write(g.read())
    return {"transcripts_fa": tfa, "reads_fastq": reads_fq, "truth_tsv": truth_tsv,
            "n_truth": str(len(truth)), "n_maf": str(len(mafs)), "n_fastq": str(n_fq)}


def main():
    ap = argparse.ArgumentParser(description="pbsim3 Tier-2 simulator wrapper")
    ap.add_argument("--transcript-fa", required=True,
                    help="FASTA of spliced transcripts (one record per transcript)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--errhmm-model", required=True,
                    help="pbsim3 packaged errhmm model (ERRHMM-ONT.model=DRS, -HQ=cDNA)")
    ap.add_argument("--depth", type=int, default=20)
    ap.add_argument("--pbsim-bin", default="pbsim")
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args()
    mafs = run_pbsim3(args.transcript_fa, args.out_prefix, args.errhmm_model,
                      depth=args.depth, pbsim_bin=args.pbsim_bin, seed=args.seed)
    print(f"generated {len(mafs)} MAF files")
    for m in mafs:
        print(m)


if __name__ == "__main__":
    main()
