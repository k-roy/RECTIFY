"""
RECTIFY core modules.

This package contains the main correction algorithms:
- bam_processor: BAM reading and 3' end extraction
- atract_detector: A-tract ambiguity detection (universal)
- ag_mispriming: AG-richness screening (oligo-dT methods)
- polya_trimmer: Poly(A) tail modeling and trimming
- indel_corrector: Indel artifact removal
- netseq_refiner: NET-seq peak matching and refinement
- spikein_filter: Spike-in RNA detection and filtering
- output_writer: Unified output format

NET-seq processing modules:
- exclusion_regions: rDNA/Pol III region exclusion
- netseq_bam_processor: NET-seq BAM processing
- netseq_deconvolution: A-tract deconvolution
- netseq_output: Output generation (parquet, bedgraph)
- netseq_command: CLI command

And CLI command implementations:
- correct_command: Main correction workflow
- train_polya_command: Poly(A) model training
- validate_command: Validation against NET-seq
- netseq_command: NET-seq processing
"""

__all__ = [
    "bam_processor",
    "atract_detector",
    "ag_mispriming",
    "polya_trimmer",
    "indel_corrector",
    "netseq_refiner",
    "spikein_filter",
    "output_writer",
    "correct_command",
    "train_polya_command",
    "validate_command",
    # NET-seq modules
    "exclusion_regions",
    "netseq_bam_processor",
    "netseq_deconvolution",
    "netseq_output",
    "netseq_command",
]
