"""
Unit tests for rectify.visualize.metagene_loaders.

Tests cover:
- iupac_to_regex: IUPAC code conversion
- loci_from_tsv: loading loci from TSV files
- loci_from_bed: BED file loading with different center modes
- loci_from_gff: GFF3 annotation loading with coordinate conversion
- loci_from_motif_scan: motif scanning on both strands
- loci_from_pickle: PKL cache loading
- position_index_from_tsv: PositionIndex construction

Strand correctness is verified for motif_scan by checking that loci produced
for + and - strands have correctly defined center positions.

Author: Kevin R. Roy
"""

import io
import pickle
import textwrap
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from rectify.visualize.metagene_loaders import (
    iupac_to_regex,
    loci_from_tsv,
    loci_from_bed,
    loci_from_gff,
    loci_from_motif_scan,
    loci_from_pickle,
    position_index_from_tsv,
    _reverse_complement,
    _center_from_match,
)
from rectify.visualize.metagene import PositionIndex


# ---------------------------------------------------------------------------
# iupac_to_regex
# ---------------------------------------------------------------------------

def test_iupac_to_regex_simple():
    """Unambiguous IUPAC bases pass through unchanged."""
    assert iupac_to_regex('ACGT') == 'ACGT'


def test_iupac_to_regex_ambiguous():
    """Ambiguous IUPAC codes expand to character classes."""
    assert iupac_to_regex('N') == '[ACGT]'
    assert iupac_to_regex('R') == '[AG]'
    assert iupac_to_regex('Y') == '[CT]'


def test_iupac_to_regex_mixed():
    """Mixed IUPAC + literals converted correctly."""
    assert iupac_to_regex('TATAAW') == 'TATAA[AT]'


def test_iupac_to_regex_passthrough_quantifiers():
    """Quantifiers and other non-IUPAC chars are preserved."""
    assert iupac_to_regex('T{5,}') == 'T{5,}'


# ---------------------------------------------------------------------------
# loci_from_tsv
# ---------------------------------------------------------------------------

def test_loci_from_tsv_basic(tmp_path):
    """Basic TSV loading with default column names."""
    tsv = tmp_path / "loci.tsv"
    tsv.write_text("chrom\tstrand\tcenter\nchrI\t+\t1000\nchrII\t-\t2000\n")

    loci = loci_from_tsv(tsv)
    assert len(loci) == 2
    assert loci[0] == {'chrom': 'chrI', 'strand': '+', 'center': 1000}
    assert loci[1] == {'chrom': 'chrII', 'strand': '-', 'center': 2000}


def test_loci_from_tsv_additional_cols(tmp_path):
    """Additional columns are included in locus dicts."""
    tsv = tmp_path / "loci.tsv"
    tsv.write_text("chrom\tstrand\tcenter\tgene\nchrI\t+\t1000\tACT1\n")

    loci = loci_from_tsv(tsv, additional_cols=['gene'])
    assert loci[0]['gene'] == 'ACT1'


def test_loci_from_tsv_custom_cols(tmp_path):
    """Custom column names work correctly."""
    tsv = tmp_path / "loci.tsv"
    tsv.write_text("chr\tstr\tpos\nchrI\t+\t500\n")

    loci = loci_from_tsv(tsv, chrom_col='chr', strand_col='str', center_col='pos')
    assert loci[0]['center'] == 500


# ---------------------------------------------------------------------------
# loci_from_bed
# ---------------------------------------------------------------------------

def test_loci_from_bed_start_center(tmp_path):
    """BED center='start': plus→bed_start, minus→bed_end-1."""
    bed = tmp_path / "features.bed"
    bed.write_text("chrI\t100\t200\tfeat1\t0\t+\nchrI\t100\t200\tfeat2\t0\t-\n")

    loci = loci_from_bed(bed, center='start')
    plus = next(l for l in loci if l['strand'] == '+')
    minus = next(l for l in loci if l['strand'] == '-')
    assert plus['center'] == 100   # bed_start
    assert minus['center'] == 199  # bed_end - 1


def test_loci_from_bed_end_center(tmp_path):
    """BED center='end': plus→bed_end-1, minus→bed_start."""
    bed = tmp_path / "features.bed"
    bed.write_text("chrI\t100\t200\tfeat1\t0\t+\nchrI\t100\t200\tfeat2\t0\t-\n")

    loci = loci_from_bed(bed, center='end')
    plus = next(l for l in loci if l['strand'] == '+')
    minus = next(l for l in loci if l['strand'] == '-')
    assert plus['center'] == 199  # bed_end - 1
    assert minus['center'] == 100  # bed_start


def test_loci_from_bed_midpoint(tmp_path):
    """BED center='midpoint' is strand-independent."""
    bed = tmp_path / "features.bed"
    bed.write_text("chrI\t100\t200\tfeat1\t0\t+\n")

    loci = loci_from_bed(bed, center='midpoint')
    assert loci[0]['center'] == 149  # (100 + 199) // 2


def test_loci_from_bed_no_strand_col(tmp_path):
    """3-column BED (no strand) uses default_strand."""
    bed = tmp_path / "regions.bed"
    bed.write_text("chrI\t500\t600\n")

    loci = loci_from_bed(bed, default_strand='-')
    assert loci[0]['strand'] == '-'


# ---------------------------------------------------------------------------
# loci_from_gff
# ---------------------------------------------------------------------------

GFF_CONTENT = textwrap.dedent("""\
    ##gff-version 3
    chrI\tSGD\tgene\t1001\t2000\t.\t+\t.\tID=gene:YAL001;Name=ACT1
    chrI\tSGD\tgene\t3001\t4000\t.\t-\t.\tID=gene:YAL002;Name=TUB1
    chrI\tSGD\tmRNA\t1001\t2000\t.\t+\t.\tID=mRNA:YAL001;Parent=gene:YAL001
""")


def test_loci_from_gff_plus_start(tmp_path):
    """GFF center='start': plus strand → 0-based feature start."""
    gff = tmp_path / "annot.gff"
    gff.write_text(GFF_CONTENT)

    loci = loci_from_gff(gff, feature_type='gene', center='start')
    plus = next(l for l in loci if l['strand'] == '+')
    # GFF start=1001 (1-based) → 0-based = 1000
    assert plus['center'] == 1000


def test_loci_from_gff_minus_start(tmp_path):
    """GFF center='start' on minus strand → 0-based feature end."""
    gff = tmp_path / "annot.gff"
    gff.write_text(GFF_CONTENT)

    loci = loci_from_gff(gff, feature_type='gene', center='start')
    minus = next(l for l in loci if l['strand'] == '-')
    # GFF end=4000 (1-based inclusive) → 0-based = 3999
    assert minus['center'] == 3999


def test_loci_from_gff_end_center(tmp_path):
    """GFF center='end': plus strand → 0-based feature end."""
    gff = tmp_path / "annot.gff"
    gff.write_text(GFF_CONTENT)

    loci = loci_from_gff(gff, feature_type='gene', center='end')
    plus = next(l for l in loci if l['strand'] == '+')
    # GFF end=2000 (1-based) → 0-based = 1999
    assert plus['center'] == 1999


def test_loci_from_gff_feature_type_filter(tmp_path):
    """feature_type filter excludes non-matching rows."""
    gff = tmp_path / "annot.gff"
    gff.write_text(GFF_CONTENT)

    gene_loci = loci_from_gff(gff, feature_type='gene')
    mrna_loci = loci_from_gff(gff, feature_type='mRNA')
    assert len(gene_loci) == 2
    assert len(mrna_loci) == 1


def test_loci_from_gff_name_attribute(tmp_path):
    """Name attribute is included in locus dict."""
    gff = tmp_path / "annot.gff"
    gff.write_text(GFF_CONTENT)

    loci = loci_from_gff(gff, feature_type='gene', name_attr='Name')
    names = {l['name'] for l in loci}
    assert names == {'ACT1', 'TUB1'}


# ---------------------------------------------------------------------------
# loci_from_motif_scan
# ---------------------------------------------------------------------------

def test_motif_scan_plus_strand_simple():
    """Plus strand: T-tract motif, center='start' → first T position."""
    # sequence: ..TTTTT.. at positions 10-14
    seq = 'A' * 10 + 'TTTTT' + 'A' * 10
    seqs = {'chrI': seq}

    loci = loci_from_motif_scan(seqs, r'T{5}', strand='+', center='start')
    assert len(loci) == 1
    assert loci[0]['chrom'] == 'chrI'
    assert loci[0]['strand'] == '+'
    assert loci[0]['center'] == 10  # first T
    assert loci[0]['match_start'] == 10
    assert loci[0]['match_end'] == 14


def test_motif_scan_plus_strand_center_end():
    """Plus strand center='end' → last position of match."""
    seq = 'A' * 10 + 'TTTTT' + 'A' * 10
    loci = loci_from_motif_scan({'chrI': seq}, r'T{5}', strand='+', center='end')
    assert loci[0]['center'] == 14  # last T


def test_motif_scan_plus_strand_midpoint():
    """Plus strand center='midpoint' → geometric midpoint."""
    seq = 'A' * 10 + 'TTTTTT' + 'A' * 10  # 6 T's at 10-15
    loci = loci_from_motif_scan({'chrI': seq}, r'T{6}', strand='+', center='midpoint')
    assert loci[0]['center'] == 12  # (10 + 15) // 2


def test_motif_scan_minus_strand_center_start():
    """
    Minus strand: A-tract on plus = T-tract on minus.
    center='start' (5' end of minus strand) → highest genomic coord of match.
    """
    # AAAAA at positions 10-14 on + strand = T-tract on - strand
    # Transcription on - goes right-to-left: 5' end is at position 14
    seq = 'T' * 10 + 'AAAAA' + 'T' * 10  # length=25
    loci = loci_from_motif_scan({'chrI': seq}, r'T{5}', strand='-', center='start')
    # RC of seq = TTTTTAAAAATTTTTTTTTT (RC of the 25-char seq)
    # T{5} matches TTTTT at RC positions 0-4 (the complement of AAAAA at 10-14)
    # Genomic 5' end = 25 - 1 - 0 = 24? Wait, let me re-check...
    # seq = 'TTTTTTTTTTAAAAATTTTTTTTTTT' (10 T's, 5 A's, 10 T's) = length 25
    # RC = complement+reverse:
    #   complement: AAAAAAAAAAAAAAATTTTTTTTTTT... wait
    #   seq:       TTTTTTTTTTEEEEEAAAAAAAAAA  (T=10, A=5, T=10)
    # Actually seq = 'T'*10 + 'AAAAA' + 'T'*10
    # complement: 'A'*10 + 'TTTTT' + 'A'*10
    # reverse:    'A'*10 + 'TTTTT' + 'A'*10 (palindrome!)
    # RC = 'A'*10 + 'TTTTT' + 'A'*10
    # T{5} match in RC at positions 10-14
    # Genomic: 5' (start) = chr_len - 1 - rc_start = 25 - 1 - 10 = 14
    assert len(loci) == 1
    assert loci[0]['strand'] == '-'
    assert loci[0]['center'] == 14  # 5' end of minus strand T-tract


def test_motif_scan_minus_strand_center_end():
    """Minus strand center='end' → lowest genomic coord (3' end)."""
    seq = 'T' * 10 + 'AAAAA' + 'T' * 10  # length=25
    loci = loci_from_motif_scan({'chrI': seq}, r'T{5}', strand='-', center='end')
    # 3' end (lowest genomic coord) = chr_len - rc_end = 25 - (10+5) = 10
    assert loci[0]['center'] == 10


def test_motif_scan_both_strands():
    """strand='both' finds matches on both strands."""
    # T-tract on plus: positions 5-9; A-tract (T-tract on minus): positions 15-19
    seq = 'A' * 5 + 'TTTTT' + 'A' * 5 + 'AAAAA' + 'A' * 5
    loci = loci_from_motif_scan({'chrI': seq}, r'T{5}', strand='both', center='start')
    plus_loci = [l for l in loci if l['strand'] == '+']
    minus_loci = [l for l in loci if l['strand'] == '-']
    assert len(plus_loci) >= 1
    assert len(minus_loci) >= 1


def test_motif_scan_iupac():
    """IUPAC motif conversion and matching."""
    # TATAAA: TATA box
    seq = 'A' * 10 + 'TATAAA' + 'A' * 10
    loci = loci_from_motif_scan(
        {'chrI': seq}, 'TATAAA', iupac=True, strand='+', center='start'
    )
    assert len(loci) == 1
    assert loci[0]['center'] == 10


def test_motif_scan_min_max_length():
    """min_length and max_length filter variable-length matches."""
    # T-tracts of varying length
    seq = 'AAA' + 'TTT' + 'AAA' + 'TTTTTTT' + 'AAA'
    loci_all = loci_from_motif_scan({'chrI': seq}, r'T+', strand='+', center='start')
    loci_min5 = loci_from_motif_scan(
        {'chrI': seq}, r'T+', strand='+', center='start', min_length=5
    )
    assert len(loci_all) == 2
    assert len(loci_min5) == 1
    assert loci_min5[0]['match_length'] == 7


def test_motif_scan_chrom_filter():
    """chroms filter limits which chromosomes are scanned."""
    seqs = {
        'chrI':   'A' * 5 + 'TTTTT' + 'A' * 5,
        'chrII':  'A' * 5 + 'TTTTT' + 'A' * 5,
    }
    loci = loci_from_motif_scan(seqs, r'T{5}', strand='+', chroms=['chrI'])
    assert all(l['chrom'] == 'chrI' for l in loci)


def test_motif_scan_match_seq_minus_strand():
    """match_seq for minus strand loci is returned in genomic orientation."""
    # AAAAA at positions 10-14 is the complement of TTTTT on minus strand
    seq = 'G' * 10 + 'AAAAA' + 'G' * 10
    loci = loci_from_motif_scan({'chrI': seq}, r'T{5}', strand='-', center='start')
    # match_seq should be AAAAA (genomic orientation of the minus strand T-tract)
    assert loci[0]['match_seq'] == 'AAAAA'


# ---------------------------------------------------------------------------
# loci_from_pickle
# ---------------------------------------------------------------------------

def test_loci_from_pickle_list_format(tmp_path):
    """Pickle with list of dicts loads correctly."""
    cache = [
        {'chrom': 'chrI', 'strand': '+', 'center': 100},
        {'chrom': 'chrI', 'strand': '-', 'center': 200},
    ]
    pkl = tmp_path / "test.pkl"
    with open(pkl, 'wb') as f:
        pickle.dump(cache, f)

    loci = loci_from_pickle(pkl)
    assert len(loci) == 2
    assert loci[0]['center'] == 100


def test_loci_from_pickle_dict_with_loci_key(tmp_path):
    """Pickle with {'loci': [...]} format."""
    cache = {
        'loci': [{'chrom': 'chrI', 'strand': '+', 'center': 500}],
        'metadata': 'some info',
    }
    pkl = tmp_path / "cache.pkl"
    with open(pkl, 'wb') as f:
        pickle.dump(cache, f)

    loci = loci_from_pickle(pkl)
    assert len(loci) == 1
    assert loci[0]['center'] == 500


def test_loci_from_pickle_additional_keys(tmp_path):
    """additional_keys are included in result dicts."""
    cache = [{'chrom': 'chrI', 'strand': '+', 'center': 100, 'gene': 'ACT1'}]
    pkl = tmp_path / "test.pkl"
    with open(pkl, 'wb') as f:
        pickle.dump(cache, f)

    loci = loci_from_pickle(pkl, additional_keys=['gene'])
    assert loci[0]['gene'] == 'ACT1'


# ---------------------------------------------------------------------------
# position_index_from_tsv
# ---------------------------------------------------------------------------

def test_position_index_from_tsv_single_file(tmp_path):
    """Single TSV file builds a PositionIndex."""
    tsv = tmp_path / "positions.tsv"
    tsv.write_text("chrom\tstrand\tposition\nchrI\t+\t1000\nchrI\t+\t1001\n")

    index, n_reads = position_index_from_tsv(tsv)
    assert n_reads == 2
    assert index.count_at('chrI', '+', 1000) == 1
    assert index.count_at('chrI', '+', 1001) == 1


def test_position_index_from_tsv_multiple_files(tmp_path):
    """Multiple TSV files are concatenated."""
    tsv1 = tmp_path / "rep1.tsv"
    tsv2 = tmp_path / "rep2.tsv"
    tsv1.write_text("chrom\tstrand\tposition\nchrI\t+\t1000\n")
    tsv2.write_text("chrom\tstrand\tposition\nchrI\t+\t1000\nchrI\t-\t2000\n")

    index, n_reads = position_index_from_tsv([tsv1, tsv2])
    assert n_reads == 3
    assert index.count_at('chrI', '+', 1000) == 2
    assert index.count_at('chrI', '-', 2000) == 1


def test_position_index_from_tsv_custom_position_col(tmp_path):
    """Custom position column name is respected."""
    tsv = tmp_path / "ends.tsv"
    tsv.write_text("chrom\tstrand\tcorrected_3prime\nchrI\t+\t5000\n")

    index, n_reads = position_index_from_tsv(tsv, position_col='corrected_3prime')
    assert index.count_at('chrI', '+', 5000) == 1


def test_position_index_from_tsv_missing_file_skipped(tmp_path):
    """Missing files are skipped with a warning."""
    tsv_exists = tmp_path / "exists.tsv"
    tsv_exists.write_text("chrom\tstrand\tposition\nchrI\t+\t100\n")
    tsv_missing = tmp_path / "missing.tsv"

    # Should succeed, skipping the missing file
    index, n_reads = position_index_from_tsv([tsv_exists, tsv_missing])
    assert n_reads == 1
