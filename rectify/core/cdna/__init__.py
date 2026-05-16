"""ONT PCR-cDNA pipeline subpackage (PCB114.24 chemistry).

Algorithms extracted from the historical
``rectify/core/commands/cdna_correct_command.py`` monolith. The command file
now contains only CLI wiring + the ``run()`` orchestration; everything
substantive lives here:

  * :mod:`._constants` — chemistry + adapter constants (SSP, UMI, ANCHOR, polyA regex)
  * :mod:`.gff`        — rDNA + gene-tree loaders; overlap-based XS classification
  * :mod:`.walkback`   — cDNA-specific walkback wrappers + bridge-G TSS correction
  * :mod:`.read_info`  — :class:`ReadInfo` + per-read UMI/orient/anchor extraction
  * :mod:`.umi`        — connected-components + umi_tools directional clustering
  * :mod:`.cluster`    — Stage-1 anchor-bucket + UMI/position clustering
  * :mod:`.consensus`  — pileup / POA / strand-aware POA + adapter pretrim
  * :mod:`.isoforms`   — Stage-2 isoform grouping + Type-1↔Type-2 pair reconciliation
  * :mod:`.io`         — BAM streaming + Stage-1 consensus FASTQ emission
"""
