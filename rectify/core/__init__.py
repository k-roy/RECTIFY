"""
RECTIFY core package.

Organized into subpackages by concern:

- ``align/``       Multi-aligner orchestration, local aligner shims, FASTQ
                   preprocessing, and mapPacBio read splitting.
- ``bam/``         BAM I/O — per-read 3'-end extraction, BAM writer,
                   processing stats, NET-seq BAM processor.
- ``consensus/``   Consensus selection from multi-aligner output:
                   ``consensus`` (selector), ``corrected_consensus``,
                   ``chimeric_consensus``.
- ``correct/``     3'-end correction core — protocol-agnostic walkback
                   (``walkback``), protocol wrappers (``protocols/``),
                   and softclip / variant-aware rescue (``indel_corrector``).
- ``splice/``      Splice junction handling — refiner, validator,
                   splice-aware 5' end, terminal exon refiner, false
                   junction filter, junction-overhang calibration.
- ``polya/``       Poly(A)-related modules — model, trimmer,
                   AG-mispriming screen, A-tract detector.
- ``netseq/``      NET-seq refinement — refiner, deconvolution, output.
                   (NET-seq BAM I/O lives in ``bam/`` alongside the other
                   BAM processors.)
- ``analyze/``     Downstream analysis modules (clustering, motif
                   discovery, APA detection, junction analysis, etc.).
- ``aggregate/``   Multi-sample aggregation.
- ``classify/``    Full-length classifier.
- ``commands/``    CLI subcommand entry points (``*_command.py``).

Cross-cutting modules at this level:
    ``exclusion_regions``, ``position_index``, ``spikein_filter``,
    ``unified_record``.
"""
