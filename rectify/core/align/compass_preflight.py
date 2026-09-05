"""Pre-flight for the COMPASS short-read panel and the deSALT arm.

🔴 **Nothing in RECTIFY builds these indices.** ``_compass_index_paths`` only *derives* where they
would be, and ``--Scer`` bundles a genome and a GTF but no index of any kind. The consequences,
all measured on Hoffman2 (planning 861 §S0c/§S3, 862):

* With ``--Scer`` all five COMPASS index paths are absent and ``run-all --short-read``
  **deadlocked for six hours** rather than degrading, because hisat2's Perl wrapper re-spawns
  itself and the orphan holds the stdout pipe open (mitigated by ``_require_compass_index``).
* An arm whose index is missing is *dropped*, and a dropped arm still exits 0 — the panel silently
  shrinks and no one is told which arms actually ran.
* 🔴 **The subtlest one: the bundled GTF has NO ``exon`` lines.** Its feature census is
  ``mRNA 11,125 · CDS 7,072 · gene 6,613 · intron 378``. STAR's ``--sjdbGTFfile`` defaults to
  ``--sjdbGTFfeatureExon exon`` and ``hisat2_extract_splice_sites.py`` reads ``exon`` too, so an
  "annotated" index built from it is **silently unannotated** — it builds without error and aligns
  at a normal rate with zero annotated junctions.

This module is diagnosis plus one repair, and both are pure functions over text:

* :func:`gtf_feature_census` / :func:`gtf_exon_status` — what the annotation actually contains.
* :func:`synthesize_exon_gtf` — write a copy with ``exon`` lines derived from ``CDS`` (grouped by
  transcript), so an annotated index CAN be built from an exon-less GTF. ⚠️ CDS-derived exons carry
  no UTRs, so the transcript ENDS are wrong; this is correct for splice-junction annotation (yeast
  introns are essentially all inside the CDS, and the GTF's own 378 ``intron`` records agree) and
  must NOT be used as a transcript model. Callers that need real transcript bounds should supply a
  proper GTF instead.
* :func:`compass_preflight` — one report naming every index the requested panel needs, whether it
  is present, and the exact command that builds it.
"""
from __future__ import annotations

import gzip
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

__all__ = [
    'ArmStatus',
    'PreflightReport',
    'compass_preflight',
    'gtf_exon_status',
    'gtf_feature_census',
    'synthesize_exon_gtf',
]

_TX_RE = re.compile(r'transcript_id\s+"([^"]+)"')
_GENE_RE = re.compile(r'gene_id\s+"([^"]+)"')
_PARENT_RE = re.compile(r'(?:^|;)\s*Parent=([^;]+)')
_ID_RE = re.compile(r'(?:^|;)\s*ID=([^;]+)')


def _open_text(path):
    p = Path(path)
    return gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p, 'rt')


def gtf_feature_census(path) -> Counter:
    """Count feature types (column 3) in a GTF/GFF, transparently handling ``.gz``.

    Comment and short lines are skipped. Returns an empty Counter for an unreadable file rather
    than raising — a pre-flight must never be the thing that stops a run.
    """
    counts: Counter = Counter()
    try:
        with _open_text(path) as fh:
            for line in fh:
                if not line or line[0] == '#':
                    continue
                f = line.rstrip('\n').split('\t')
                if len(f) >= 9:
                    counts[f[2]] += 1
    except OSError:
        return Counter()
    return counts


def gtf_exon_status(path) -> Tuple[str, Counter]:
    """Classify an annotation for annotated-index building.

    Returns ``(status, census)`` where status is:

    * ``'exon'``     — has ``exon`` features; STAR/HISAT2 defaults work as documented.
    * ``'cds_only'`` — no ``exon`` but has ``CDS``; an annotated index built with the DEFAULT
      feature name is silently unannotated. Repairable with :func:`synthesize_exon_gtf`.
    * ``'unusable'`` — neither; nothing can derive splice sites from it.
    * ``'missing'``  — no annotation given, or the file could not be read.
    """
    if path is None:
        return 'missing', Counter()
    census = gtf_feature_census(path)
    if not census:
        return 'missing', census
    if census.get('exon', 0) > 0:
        return 'exon', census
    if census.get('CDS', 0) > 0:
        return 'cds_only', census
    return 'unusable', census


def _transcript_key(attrs: str) -> Optional[str]:
    """Transcript id from GTF ``transcript_id "X"`` or GFF3 ``Parent=X`` (first parent only)."""
    m = _TX_RE.search(attrs)
    if m:
        return m.group(1)
    m = _PARENT_RE.search(attrs)
    if m:
        return m.group(1).split(',')[0]
    m = _GENE_RE.search(attrs) or _ID_RE.search(attrs)
    return m.group(1) if m else None


def synthesize_exon_gtf(src, dst) -> Dict[str, int]:
    """Write *dst*: *src* plus one ``exon`` line per ``CDS`` line, keyed by transcript.

    The synthesized lines copy the CDS coordinates, strand and attributes verbatim and only change
    column 3 to ``exon``, so a downstream ``hisat2_extract_splice_sites.py`` or
    ``STAR --sjdbGTFfile`` sees the same intron structure the CDS records imply.

    ⚠️ **CDS-derived exons have no UTRs.** Transcript starts and ends in *dst* are ORF bounds, not
    transcript bounds. That is correct for splice-site extraction and wrong for anything measuring
    a distance from a transcript terminus (the SGD ``gene``-vs-``CDS`` trap in reverse).

    Returns ``{'cds_in': N, 'exon_written': N, 'transcripts': N, 'exon_already': N}``. When *src*
    already has ``exon`` features nothing is synthesized and *dst* is a byte copy, so calling this
    unconditionally is safe.
    """
    status, census = gtf_exon_status(src)
    out_lines: List[str] = []
    n_cds = 0
    tx = set()
    with _open_text(src) as fh:
        for line in fh:
            out_lines.append(line if line.endswith('\n') else line + '\n')
            if not line or line[0] == '#':
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS' or status == 'exon':
                continue
            n_cds += 1
            key = _transcript_key(f[8])
            if key:
                tx.add(key)
            g = list(f)
            g[2] = 'exon'
            out_lines.append('\t'.join(g) + '\n')
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(''.join(out_lines))
    return {
        'cds_in': census.get('CDS', 0),
        'exon_written': n_cds,
        'transcripts': len(tx),
        'exon_already': census.get('exon', 0),
    }


@dataclass
class ArmStatus:
    """One aligner arm's index: where it should be, whether it is there, how to build it."""
    arm: str
    path: Path
    present: bool
    build_cmd: str
    note: str = ''


@dataclass
class PreflightReport:
    arms: List[ArmStatus] = field(default_factory=list)
    annotation_status: str = 'missing'
    annotation_census: Counter = field(default_factory=Counter)
    annotation_path: Optional[Path] = None

    @property
    def missing(self) -> List[ArmStatus]:
        return [a for a in self.arms if not a.present]

    @property
    def ok(self) -> bool:
        return not self.missing and self.annotation_status in ('exon', 'missing')

    def render(self) -> str:
        """A single readable block — the point is that a dropped arm is never silent again."""
        lines = ['COMPASS / deSALT index pre-flight']
        for a in self.arms:
            mark = 'OK     ' if a.present else 'MISSING'
            lines.append(f'  [{mark}] {a.arm:<20} {a.path}')
            if not a.present:
                lines.append(f'            build: {a.build_cmd}')
            if a.note:
                lines.append(f'            note:  {a.note}')
        if self.annotation_status == 'cds_only':
            lines.append(
                '  [WARNING] the annotation has NO "exon" features '
                f'({dict(self.annotation_census)}). STAR --sjdbGTFfile and '
                'hisat2_extract_splice_sites.py both read "exon" by default, so an annotated index '
                'built from it is SILENTLY UNANNOTATED — it builds without error and aligns '
                'normally with zero annotated junctions. Fix: build from a GTF with exon lines '
                '(rectify can synthesize one from CDS), or pass STAR '
                '--sjdbGTFfeatureExon CDS.')
        elif self.annotation_status == 'unusable':
            lines.append(
                '  [WARNING] the annotation has neither "exon" nor "CDS" features '
                f'({dict(self.annotation_census)}); no splice sites can be derived from it.')
        if self.missing:
            lines.append(f'  => {len(self.missing)} arm(s) would be DROPPED. A dropped arm still '
                         'exits 0; pass --require-compass-index to make this fatal.')
        return '\n'.join(lines)


def compass_preflight(
    genome_path,
    annotation_path=None,
    *,
    read_length: int = 150,
    arms: Optional[List[str]] = None,
    desalt: bool = False,
) -> PreflightReport:
    """Check every index the requested arms need and report it in one block.

    *arms* names the COMPASS arms to check (default: the full short-read panel). Index locations
    come from ``multi_aligner._compass_index_paths`` so this can never drift from what the aligners
    actually open.
    """
    from .multi_aligner import _compass_index_paths

    p = _compass_index_paths(genome_path, read_length=read_length)
    g = Path(genome_path)
    want = arms or ['STAR', 'HISAT2', 'splice_sites', 'magicblast', 'gsnap']
    rep = PreflightReport(annotation_path=Path(annotation_path) if annotation_path else None)
    rep.annotation_status, rep.annotation_census = gtf_exon_status(annotation_path)
    gtf = str(annotation_path) if annotation_path else '<annotation.gtf>'

    if 'STAR' in want:
        rep.arms.append(ArmStatus(
            'STAR', p.star_dir, (p.star_dir / 'SAindex').exists(),
            f'STAR --runMode genomeGenerate --genomeDir {p.star_dir} --genomeFastaFiles {g} '
            f'--sjdbGTFfile {gtf} --sjdbOverhang {read_length - 1} --genomeSAindexNbases 11',
            note='SAindex is the completion sentinel; a half-built dir looks present to a naive check.'))
    if 'HISAT2' in want:
        rep.arms.append(ArmStatus(
            'HISAT2', Path(str(p.hisat2_index) + '.1.ht2'),
            Path(str(p.hisat2_index) + '.1.ht2').exists(),
            f'hisat2-build -p 4 {g} {p.hisat2_index}'))
    if 'splice_sites' in want:
        rep.arms.append(ArmStatus(
            'HISAT2 splice sites', p.splice_sites, Path(p.splice_sites).exists(),
            f'hisat2_extract_splice_sites.py {gtf} > {p.splice_sites}',
            note='hisat2 accepts a MISSING --known-splicesite-infile silently: exit 0, normal '
                 'alignment rate, zero annotated junctions.'))
    if 'magicblast' in want:
        rep.arms.append(ArmStatus(
            'magicblast (BLAST db)', Path(str(p.blast_index) + '.nsq'),
            Path(str(p.blast_index) + '.nsq').exists(),
            f'makeblastdb -in {g} -dbtype nucl -parse_seqids -out {p.blast_index}'))
    if 'gsnap' in want:
        rep.arms.append(ArmStatus(
            'gsnap', p.gsnap_dir, p.gsnap_dir.is_dir() and any(p.gsnap_dir.iterdir())
            if p.gsnap_dir.is_dir() else False,
            f'gmap_build -D {p.gsnap_dir.parent} -d {p.genome_version} {g}'))
    if desalt:
        d = g.parent / 'desalt_index'
        rep.arms.append(ArmStatus(
            'deSALT', d, (d / 'ref.seq').exists() and (d / 'ref.seq').stat().st_size > 0
            if (d / 'ref.seq').exists() else False,
            f'deSALT index {g} {d}',
            note='ref.seq is the completion sentinel; an empty desalt_index/ placeholder reads as '
                 'present to a directory check.'))
    return rep
