#!/usr/bin/env python3
"""
RECTIFY 'tag-polya' command implementation.

Adds poly(A) evidence tags to an already-aligned BAM file:

  ps:f  RECTIFY poly(A) confidence score (0–1).  Always written.
  pt:i  Sequence-based poly(A) tail length estimate.  Written only when the
        read does not already carry a dorado pt:i tag.  Existing dorado values
        are never overwritten.

Intended use-cases:
  1. Retroactive annotation of roadblocks BAMs produced before dorado pt:i
     support (or processed through pipelines that strip aux tags).
  2. Any pipeline that produces aligned BAMs without pt:i (FASTQ → minimap2
     without the -y tag-passthrough flag).

Author: Kevin R. Roy
Date: 2026-04-11
"""

import logging
import sys
from pathlib import Path
from typing import Optional

import pysam

from ..utils.genome import reverse_complement
from ..utils.alignment import extract_soft_clips
from .polya_model import PolyAModel, load_model as load_polya_model, get_default_model


logger = logging.getLogger(__name__)


def _score_read(
    read: pysam.AlignedSegment,
    strand: str,
    polya_model: PolyAModel,
) -> tuple:
    """
    Score the 3' soft-clip of a read using the poly(A) model.

    Returns:
        (polya_score, soft_clip_a_count)
        polya_score: float confidence 0-1, or None if no soft-clip
        soft_clip_a_count: number of A's in the 3' soft-clip (sequence-based
            tail length estimate), or None if no soft-clip
    """
    soft_clips = extract_soft_clips(read)

    three_prime_clip = None
    for clip in soft_clips:
        if (strand == '+' and clip['side'] == 'right') or \
                (strand == '-' and clip['side'] == 'left'):
            three_prime_clip = clip
            break

    if three_prime_clip is None or not three_prime_clip.get('seq'):
        return None, None

    clip_seq = three_prime_clip['seq']
    if strand == '-':
        clip_seq = reverse_complement(clip_seq)

    clip_seq = clip_seq.upper()
    score_result = polya_model.score_sequence(clip_seq)
    polya_score = round(score_result['confidence'], 4)
    soft_clip_a_count = clip_seq.count('A')

    return polya_score, soft_clip_a_count


def run(args) -> int:
    """
    Execute the 'tag-polya' command.

    Args:
        args: Parsed argparse namespace with fields:
            bam, output, polya_model, threads, verbose

    Returns:
        Exit code (0 = success, 1 = error)
    """
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    bam_path = Path(args.bam)
    output_path = Path(args.output)
    threads = args.threads

    if not bam_path.exists():
        logger.error(f"Input BAM not found: {bam_path}")
        return 1

    if not output_path.parent.exists():
        logger.error(f"Output directory does not exist: {output_path.parent}")
        return 1

    # Load poly(A) model
    model_path = getattr(args, 'polya_model', None)
    if model_path is not None:
        polya_model = load_polya_model(Path(model_path))
        if polya_model is None:
            logger.error(f"Failed to load poly(A) model from {model_path}")
            return 1
        logger.info(f"Loaded poly(A) model from {model_path}")
    else:
        polya_model = get_default_model()
        logger.info("Using built-in default poly(A) model")

    logger.info(f"Input BAM:  {bam_path}")
    logger.info(f"Output BAM: {output_path}")

    n_total = 0
    n_scored = 0
    n_pt_written = 0
    n_pt_preserved = 0

    try:
        with pysam.AlignmentFile(str(bam_path), 'rb') as bam_in:
            with pysam.AlignmentFile(
                str(output_path), 'wb',
                header=bam_in.header,
                threads=threads,
            ) as bam_out:
                for read in bam_in:
                    n_total += 1

                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        bam_out.write(read)
                        continue

                    strand = '-' if read.is_reverse else '+'

                    # Check for existing dorado pt:i tag — never overwrite
                    try:
                        _existing_pt = read.get_tag('pt')
                        has_pt_tag = True
                        n_pt_preserved += 1
                    except KeyError:
                        _existing_pt = None
                        has_pt_tag = False

                    # Score the 3' soft-clip
                    polya_score, soft_clip_a_count = _score_read(read, strand, polya_model)

                    if polya_score is not None:
                        n_scored += 1
                        # ps:f — RECTIFY poly(A) confidence (always written)
                        read.set_tag('ps', polya_score, value_type='f')

                        # pt:i — sequence-based tail length estimate (only when absent)
                        if not has_pt_tag and soft_clip_a_count is not None:
                            read.set_tag('pt', soft_clip_a_count, value_type='i')
                            n_pt_written += 1

                    bam_out.write(read)

                    if n_total % 500000 == 0:
                        logger.info(f"  Processed {n_total:,} reads...")

    except Exception as e:
        logger.error(f"Error processing BAM: {e}")
        return 1

    logger.info(f"Done. {n_total:,} reads processed.")
    logger.info(f"  Scored (ps:f written): {n_scored:,}")
    logger.info(f"  pt:i written (seq-based): {n_pt_written:,}")
    logger.info(f"  pt:i preserved (dorado): {n_pt_preserved:,}")

    # Index output BAM
    try:
        pysam.sort('-o', str(output_path) + '.sorted.bam', str(output_path))
        import os
        os.rename(str(output_path) + '.sorted.bam', str(output_path))
        pysam.index(str(output_path), f'-@ {threads}')
        logger.info(f"Indexed output BAM: {output_path}.bai")
    except Exception as e:
        logger.warning(f"Could not index output BAM: {e} — run 'samtools index' manually")

    return 0
