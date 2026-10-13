"""Report-only flag for strand-switching-primer (SSP) strand invasion in ONT PCR-cDNA (ISSUE-086).

During reverse transcription the SSP can anneal to an internal stretch of the RNA that matches its own sequence; the
RT then jumps onto it and the molecule begins at that stretch. Its aligned 5' end then sits on a genomic copy of part
of the SSP, while the rest of the primer and the UMI stay soft-clipped. CDC48 (YDL126C) is the case that showed it:
233 Ysh1-AA molecules piled their 5' ends at chrIV:237,413, where the genome reads GGTGCTGATATTGC in RNA sense,
SSP positions 9-22; direct RNA has no pile there (Chanfreau planning/1036, Kevin's review card C09).

This module only measures. ``ssp_match_length`` is the length of the longest prefix of the genome's RNA-sense
sequence at a molecule's aligned 5' end that occurs in the SSP. Dropping or down-weighting flagged molecules, and the
match length that defines "flagged", are Kevin's decisions; the default length is a reporting threshold, chosen so
that a uniform random genome reaches it about once per million molecules (``chance_rate``).
"""

from typing import Optional

from ._constants import SSP_FWD

SSP_FLAG_MIN_NT = 12        # report-only threshold; see chance_rate(12) ~= 7e-7 per molecule
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def ssp_match_length(rna_5p_seq: str, ssp: str = SSP_FWD) -> int:
    """Longest prefix of *rna_5p_seq* (RNA sense, starting at the aligned 5' end) found anywhere in *ssp*."""
    seq = (rna_5p_seq or "").upper()
    best = 0
    for n in range(1, min(len(seq), len(ssp)) + 1):
        if seq[:n] in ssp:
            best = n
        else:
            break
    return best


def chance_rate(n: int, ssp: str = SSP_FWD) -> float:
    """Chance that a uniform random genome starts with an n-mer found in *ssp* (distinct n-mers / 4**n)."""
    if n <= 0:
        return 1.0
    kmers = {ssp[i:i + n] for i in range(len(ssp) - n + 1)}
    return len(kmers) / float(4 ** n)


def rna_5p_genome(chrom_seq: str, aln_start: int, aln_end: int, rna_reverse: bool,
                  k: int = len(SSP_FWD)) -> Optional[str]:
    """The genome's first *k* bases in RNA sense from a molecule's aligned 5' end.

    *rna_reverse* is True when the RNA runs on the minus strand of the reference (poly(T) on the left, the
    'rev' frame): the 5' end is then the last aligned base and the sequence is reverse-complemented.
    Returns None when the window leaves the contig.
    """
    if rna_reverse:
        lo = aln_end - k
        if lo < 0 or aln_end > len(chrom_seq):
            return None
        return chrom_seq[lo:aln_end].translate(_COMPLEMENT)[::-1].upper()
    hi = aln_start + k
    if aln_start < 0 or hi > len(chrom_seq):
        return None
    return chrom_seq[aln_start:hi].upper()
