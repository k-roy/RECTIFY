#!/usr/bin/env python3
"""Multi-aligner consensus tags (Xa/Xc/Xn/Xt) → four corrected-TSV columns.

Ledger item D3.  ``rectify correct`` reads the align-stage
``<prefix>.multialigned.bam`` (written by ``run_consensus_selection``,
``align_command.py:1212``; handed to correct by ``run/stages.py:79-89``), whose
winning records carry ``Xa`` (best aligner), ``Xc`` (confidence), ``Xn``
(aligners agreeing) and ``Xt`` (tied aligners) from
``consensus/consensus.py:971-975``.  None of it reached the per-read TSV.

The contract pinned here:

* four SEPARATE columns — ``consensus_aligner``, ``consensus_confidence``,
  ``consensus_n_agree``, ``consensus_tied`` — never a derived concordance
  class, because a downstream browser consumes the raw values live;
* appended LAST, so every pre-existing column keeps its index;
* empty string, never a placeholder, when the tag is absent — a per-aligner
  BAM (the correct-first path) carries none and must still work unchanged.
"""
from __future__ import annotations

from io import StringIO

import pysam
import pytest
from unittest.mock import MagicMock

from rectify.core.bam.output import (
    CONSENSUS_TAG_COLUMNS,
    CORRECTION_TSV_HEADER,
    consensus_tag_fields,
    correction_result_to_tsv_row,
)


CHROM = "chrTest"
CHROM_LEN = 6000


def _header():
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": CHROM_LEN}]}
    )


def _read(
    start: int = 1000,
    length: int = 200,
    *,
    query_name: str = "r",
    tags: dict | None = None,
    header=None,
):
    """Gapless plus-strand segment spanning [start, start+length)."""
    a = pysam.AlignedSegment(header or _header())
    a.query_name = query_name
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = 60
    a.flag = 0
    a.query_sequence = "C" * length
    a.query_qualities = pysam.qualitystring_to_array("I" * length)
    a.cigartuples = [(0, length)]
    for name, value in (tags or {}).items():
        a.set_tag(name, value)
    return a


def _flat_genome():
    """Non-A/T genome so no walkback or A-tract module shifts the terminus."""
    return {CHROM: "C" * CHROM_LEN}


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------
class TestSchema:
    def test_column_names_and_order(self):
        assert CONSENSUS_TAG_COLUMNS == (
            "consensus_aligner",
            "consensus_confidence",
            "consensus_n_agree",
            "consensus_tied",
        )

    def test_columns_are_a_contiguous_block_after_strand_evidence(self):
        """Appended as one block right after the last pre-D3 column, so every
        pre-existing column keeps its index. (Integration 2026-09-05: a later
        schema addition — ISSUE-002's ``five_prime_rescue_refused`` — now sits
        AFTER this block; the pin is the block's position, not the header's tail.)"""
        i = CORRECTION_TSV_HEADER.index("strand_evidence")
        assert tuple(CORRECTION_TSV_HEADER[i + 1:i + 5]) == CONSENSUS_TAG_COLUMNS

    def test_header_has_no_duplicates(self):
        assert len(set(CORRECTION_TSV_HEADER)) == len(CORRECTION_TSV_HEADER)

    def test_strand_evidence_kept_its_index(self):
        """The column that used to be last must not have moved — pinned by its
        ABSOLUTE index (38), which is what positional consumers depend on."""
        assert CORRECTION_TSV_HEADER[38] == "strand_evidence"


# ---------------------------------------------------------------------------
# consensus_tag_fields
# ---------------------------------------------------------------------------
class TestConsensusTagFields:
    def test_all_four_tags(self):
        read = _read(tags={"Xa": "minimap2", "Xc": "high", "Xn": 3,
                           "Xt": "deSALT,minimap2"})
        assert consensus_tag_fields(read) == {
            "consensus_aligner": "minimap2",
            "consensus_confidence": "high",
            "consensus_n_agree": "3",       # Xn is an int tag → decimal string
            "consensus_tied": "deSALT,minimap2",
        }

    def test_untied_read_has_empty_tied_column(self):
        """``Xt`` is set only when the vote tied (consensus.py:974-975)."""
        read = _read(tags={"Xa": "uLTRA", "Xc": "medium", "Xn": 1})
        fields = consensus_tag_fields(read)
        assert fields["consensus_aligner"] == "uLTRA"
        assert fields["consensus_n_agree"] == "1"
        assert fields["consensus_tied"] == ""

    def test_untagged_read_is_all_empty_strings(self):
        fields = consensus_tag_fields(_read())
        assert fields == {c: "" for c in CONSENSUS_TAG_COLUMNS}
        assert all(isinstance(v, str) for v in fields.values())

    def test_never_raises_when_get_tag_raises(self):
        """Exception-free by contract: ``correct_read_3prime`` calls this next
        to a ``try/except KeyError`` that would swallow an escaping error."""
        read = MagicMock()
        read.has_tag.return_value = True
        read.get_tag.side_effect = KeyError("boom")
        assert consensus_tag_fields(read) == {c: "" for c in CONSENSUS_TAG_COLUMNS}

    def test_never_raises_on_a_record_without_tag_accessors(self):
        assert consensus_tag_fields(object()) == {
            c: "" for c in CONSENSUS_TAG_COLUMNS
        }


# ---------------------------------------------------------------------------
# Row serialization
# ---------------------------------------------------------------------------
class TestRowSerialization:
    def _row_by_column(self, result):
        fields = correction_result_to_tsv_row(result)
        assert len(fields) == len(CORRECTION_TSV_HEADER)
        return dict(zip(CORRECTION_TSV_HEADER, fields))

    def _minimal_result(self, **extra):
        result = {
            "read_id": "r1",
            "chrom": CHROM,
            "strand": "+",
            "original_3prime": 100,
            "corrected_3prime": 100,
            "ambiguity_min": 100,
            "ambiguity_max": 100,
            "ambiguity_range": 0,
            "correction_applied": [],
            "confidence": "high",
            "qc_flags": [],
        }
        result.update(extra)
        return result

    def test_values_land_in_their_own_columns(self):
        by_col = self._row_by_column(self._minimal_result(
            consensus_aligner="deSALT",
            consensus_confidence="high",
            consensus_n_agree="4",
            consensus_tied="",
        ))
        assert by_col["consensus_aligner"] == "deSALT"
        assert by_col["consensus_confidence"] == "high"
        assert by_col["consensus_n_agree"] == "4"
        assert by_col["consensus_tied"] == ""

    def test_missing_keys_serialize_as_empty_strings(self):
        """A result dict from an older/partial producer must still round-trip."""
        by_col = self._row_by_column(self._minimal_result())
        for column in CONSENSUS_TAG_COLUMNS:
            assert by_col[column] == ""

    def test_none_serializes_as_empty_string_not_the_word_none(self):
        by_col = self._row_by_column(self._minimal_result(
            consensus_aligner=None, consensus_tied=None,
        ))
        assert by_col["consensus_aligner"] == ""
        assert by_col["consensus_tied"] == ""

    def test_int_n_agree_is_stringified(self):
        by_col = self._row_by_column(self._minimal_result(consensus_n_agree=2))
        assert by_col["consensus_n_agree"] == "2"

    def test_no_tabs_or_newlines_injected(self):
        fields = correction_result_to_tsv_row(self._minimal_result(
            consensus_aligner="minimap2", consensus_tied="deSALT,minimap2",
        ))
        assert all("\t" not in f and "\n" not in f for f in fields)

    def test_write_results_chunk_places_the_columns(self):
        from rectify.core.bam.parallel import _write_results_chunk

        out = StringIO()
        _write_results_chunk(out, [self._minimal_result(
            consensus_aligner="minimap2",
            consensus_confidence="low",
            consensus_n_agree="2",
            consensus_tied="gapmm2,minimap2",
        )])
        fields = out.getvalue().rstrip("\n").split("\t")
        assert len(fields) == len(CORRECTION_TSV_HEADER)
        by_col = dict(zip(CORRECTION_TSV_HEADER, fields))
        assert by_col["consensus_aligner"] == "minimap2"
        assert by_col["consensus_confidence"] == "low"
        assert by_col["consensus_n_agree"] == "2"
        assert by_col["consensus_tied"] == "gapmm2,minimap2"


# ---------------------------------------------------------------------------
# correct_read_3prime — main path and the chimeric short-circuit
# ---------------------------------------------------------------------------
class TestCorrectReadWiring:
    def test_main_path_carries_all_four(self):
        from rectify.core.bam.bam_processor import correct_read_3prime

        read = _read(tags={"Xa": "minimap2", "Xc": "high", "Xn": 5,
                           "Xt": "deSALT,minimap2"})
        rows = correct_read_3prime(read, _flat_genome(), apply_3ss_rescue=False)
        assert len(rows) == 1
        assert rows[0]["consensus_aligner"] == "minimap2"
        assert rows[0]["consensus_confidence"] == "high"
        assert rows[0]["consensus_n_agree"] == "5"
        assert rows[0]["consensus_tied"] == "deSALT,minimap2"

    def test_main_path_untied_leaves_tied_empty(self):
        from rectify.core.bam.bam_processor import correct_read_3prime

        read = _read(tags={"Xa": "gapmm2", "Xc": "medium", "Xn": 1})
        rows = correct_read_3prime(read, _flat_genome(), apply_3ss_rescue=False)
        assert rows[0]["consensus_aligner"] == "gapmm2"
        assert rows[0]["consensus_tied"] == ""

    def test_single_aligner_bam_is_unchanged(self):
        """No consensus stage upstream → four empty strings, nothing else moves."""
        from rectify.core.bam.bam_processor import correct_read_3prime

        read = _read()
        rows = correct_read_3prime(read, _flat_genome(), apply_3ss_rescue=False)
        assert len(rows) == 1
        for column in CONSENSUS_TAG_COLUMNS:
            assert rows[0][column] == ""

    def test_chimeric_untagged_still_short_circuits(self):
        """Pins the KeyError hazard: ``consensus_tag_fields`` is read ABOVE the
        chimeric ``try/except KeyError`` (bam_processor.py:388-395).  Were it
        read inside and were it to raise, the read would silently fall through
        to the non-chimeric path — a wrong row, not a missing column."""
        from rectify.core.bam.bam_processor import correct_read_3prime

        read = _read(tags={"Xz": 1})
        rows = correct_read_3prime(read, _flat_genome(), apply_3ss_rescue=False)
        assert len(rows) == 1
        assert rows[0]["confidence"] == "chimeric"
        for column in CONSENSUS_TAG_COLUMNS:
            assert rows[0][column] == ""

    def test_chimeric_comma_list_aligner_passes_through_verbatim(self):
        """``chimeric_consensus.py:986-989`` joins the per-segment winners into
        a comma list; the column carries it as-is (the TSV is tab-separated,
        and correction_applied/qc_flags are already comma-joined)."""
        from rectify.core.bam.bam_processor import correct_read_3prime

        read = _read(tags={"Xz": 1, "Xa": "minimap2,deSALT", "Xc": "medium"})
        rows = correct_read_3prime(read, _flat_genome(), apply_3ss_rescue=False)
        assert rows[0]["confidence"] == "chimeric"
        assert rows[0]["consensus_aligner"] == "minimap2,deSALT"
        assert rows[0]["consensus_confidence"] == "medium"
        # Chimeric records get no Xn/Xt (chimeric_consensus.py sets neither).
        assert rows[0]["consensus_n_agree"] == ""
        assert rows[0]["consensus_tied"] == ""

    def test_netseq_style_row_copies_inherit(self):
        """Multi-peak rows are ``dict(result)`` shallow copies
        (bam_processor.py:1187-1206) — they must carry the columns too."""
        from rectify.core.bam.bam_processor import correct_read_3prime

        read = _read(tags={"Xa": "uLTRA", "Xc": "high", "Xn": 2})
        rows = correct_read_3prime(read, _flat_genome(), apply_3ss_rescue=False)
        for row in rows:
            assert row["consensus_aligner"] == "uLTRA"


# ---------------------------------------------------------------------------
# End-to-end: a small tagged BAM → corrected TSV
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def tagged_tsv(tmp_path_factory):
    """Stream a 3-read BAM (tagged / partially tagged / untagged) through
    ``process_bam_streaming`` and return the parsed corrected TSV rows."""
    from rectify.core.bam.parallel import process_bam_streaming

    tmp = tmp_path_factory.mktemp("consensus_tags")

    genome_path = tmp / "genome.fa"
    genome_path.write_text(f">{CHROM}\n" + "C" * CHROM_LEN + "\n")

    bam_path = tmp / "tagged.bam"
    header = _header()
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        out.write(_read(1000, 200, query_name="tied_read", header=header,
                        tags={"Xa": "minimap2", "Xc": "high", "Xn": 3,
                              "Xt": "deSALT,minimap2"}))
        out.write(_read(2000, 200, query_name="untied_read", header=header,
                        tags={"Xa": "deSALT", "Xc": "low", "Xn": 1}))
        out.write(_read(3000, 200, query_name="untagged_read", header=header))
    pysam.index(str(bam_path))

    out_tsv = tmp / "corrected_reads.tsv"
    process_bam_streaming(
        bam_path=str(bam_path),
        genome_path=str(genome_path),
        output_path=str(out_tsv),
        show_progress=False,
        apply_3ss_rescue=False,
    )

    lines = out_tsv.read_text().rstrip("\n").split("\n")
    header_cols = lines[0].split("\t")
    rows = {}
    for line in lines[1:]:
        fields = line.split("\t")
        assert len(fields) == len(header_cols)
        by_col = dict(zip(header_cols, fields))
        rows[by_col["read_id"]] = by_col
    return header_cols, rows


class TestEndToEndTaggedBam:
    def test_header_carries_the_four_columns_last(self, tagged_tsv):
        header_cols, _ = tagged_tsv
        assert header_cols == CORRECTION_TSV_HEADER

    def test_all_three_reads_survive(self, tagged_tsv):
        _, rows = tagged_tsv
        assert set(rows) == {"tied_read", "untied_read", "untagged_read"}

    def test_tied_read(self, tagged_tsv):
        _, rows = tagged_tsv
        row = rows["tied_read"]
        assert row["consensus_aligner"] == "minimap2"
        assert row["consensus_confidence"] == "high"
        assert row["consensus_n_agree"] == "3"
        assert row["consensus_tied"] == "deSALT,minimap2"

    def test_untied_read_has_empty_tied_column(self, tagged_tsv):
        _, rows = tagged_tsv
        row = rows["untied_read"]
        assert row["consensus_aligner"] == "deSALT"
        assert row["consensus_confidence"] == "low"
        assert row["consensus_n_agree"] == "1"
        assert row["consensus_tied"] == ""

    def test_untagged_read_has_four_empty_columns(self, tagged_tsv):
        _, rows = tagged_tsv
        row = rows["untagged_read"]
        for column in CONSENSUS_TAG_COLUMNS:
            assert row[column] == ""

    def test_pandas_reads_it_without_shifting_columns(self, tagged_tsv):
        """An empty trailing field must not be read as a missing column —
        the failure mode that shifts every name one position right."""
        import pandas as pd

        _, rows = tagged_tsv
        assert rows  # fixture built
        df = pd.DataFrame(rows).T
        assert list(df.columns) == CORRECTION_TSV_HEADER
        assert df.loc["untagged_read", "consensus_aligner"] == ""


# ---------------------------------------------------------------------------
# merge_corrected_tsvs must carry the columns through, not synthesize them
# ---------------------------------------------------------------------------
_MERGE_COLUMNS = [
    'read_id', 'chrom', 'strand',
    'original_3prime', 'corrected_3prime',
    'five_prime_position', 'five_prime_rescued',
    'alignment_start', 'alignment_end',
    'ambiguity_min', 'ambiguity_max', 'ambiguity_range',
    'polya_length', 'aligned_a_length', 'soft_clip_a_length',
    'junctions', 'n_junctions',
    'five_prime_soft_clip_length', 'three_prime_soft_clip_length',
    'mapq', 'correction_applied', 'confidence', 'qc_flags', 'fraction',
] + list(CONSENSUS_TAG_COLUMNS)


def _merge_row(read_id, chrom, pos, *, span=100, **consensus):
    row = {
        'read_id': read_id, 'chrom': chrom, 'strand': '+',
        'original_3prime': pos, 'corrected_3prime': pos,
        'five_prime_position': 0, 'five_prime_rescued': False,
        'alignment_start': 0, 'alignment_end': span,
        'ambiguity_min': pos, 'ambiguity_max': pos, 'ambiguity_range': 0,
        'polya_length': 0, 'aligned_a_length': 0, 'soft_clip_a_length': 0,
        'junctions': '', 'n_junctions': 0,
        'five_prime_soft_clip_length': 0, 'three_prime_soft_clip_length': 0,
        'mapq': 60, 'correction_applied': '', 'confidence': 'high',
        'qc_flags': '', 'fraction': 1.0,
    }
    row.update({c: '' for c in CONSENSUS_TAG_COLUMNS})
    row.update(consensus)
    return row


class TestMergeCarriesColumnsThrough:
    def _merge(self, tmp_path, per_aligner_rows):
        import pandas as pd
        from rectify.core.consensus.corrected_consensus import merge_corrected_tsvs

        per_aligner_tsvs = {}
        for aligner, rows in per_aligner_rows.items():
            path = tmp_path / f"{aligner}.tsv"
            pd.DataFrame(rows, columns=_MERGE_COLUMNS).to_csv(
                path, sep='\t', index=False
            )
            per_aligner_tsvs[aligner] = path
        out_tsv = tmp_path / "merged.tsv"
        merge_corrected_tsvs(
            per_aligner_tsvs=per_aligner_tsvs,
            output_tsv=out_tsv,
            per_aligner_corrected_bams=None,
        )
        return pd.read_csv(out_tsv, sep='\t', keep_default_na=False)

    def test_winner_keeps_its_consensus_columns(self, tmp_path):
        """A merge is a row SELECTION — the four columns ride along with the
        winning aligner's row and are not recomputed from the merge's own vote
        (that decision is already published as ``winning_aligner``)."""
        merged = self._merge(tmp_path, {
            'minimap2': [_merge_row(
                'r1', 'chrI', 500, span=200,
                consensus_aligner='minimap2', consensus_confidence='high',
                consensus_n_agree='3', consensus_tied='',
            )],
            'deSALT': [_merge_row(
                'r1', 'chrI', 500, span=100,
                consensus_aligner='minimap2', consensus_confidence='high',
                consensus_n_agree='3', consensus_tied='',
            )],
        })
        assert len(merged) == 1
        for column in CONSENSUS_TAG_COLUMNS:
            assert column in merged.columns, f"{column} dropped by the merge"
        row = merged.iloc[0]
        assert row['consensus_aligner'] == 'minimap2'
        assert row['consensus_confidence'] == 'high'
        assert str(row['consensus_n_agree']) == '3'
        assert row['consensus_tied'] == ''
        # winning_aligner stays the merge's OWN answer — a separate column.
        assert 'winning_aligner' in merged.columns

    def test_empty_columns_survive_the_correct_first_path(self, tmp_path):
        """Per-aligner correct (the correct-first order) sees no consensus
        tags; the merged TSV must still carry four empty columns."""
        merged = self._merge(tmp_path, {
            'minimap2': [_merge_row('r1', 'chrI', 500, span=200)],
            'deSALT': [_merge_row('r1', 'chrI', 500, span=100)],
        })
        assert len(merged) == 1
        for column in CONSENSUS_TAG_COLUMNS:
            assert merged.iloc[0][column] == ''


# ---------------------------------------------------------------------------
# The merge-time join: fill the columns from the align-stage multialigned BAM
# ---------------------------------------------------------------------------
def _write_consensus_bam(path, reads):
    """*reads* = [(query_name, {tag: value}), ...] → a tiny sorted+indexed BAM."""
    header = _header()
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for offset, (name, tags) in enumerate(reads):
            out.write(_read(1000 + 300 * offset, 200, query_name=name,
                            tags=tags, header=header))
    pysam.index(str(path))
    return path


class TestMergeJoinsConsensusBam:
    """`correct` runs on PER-ALIGNER BAMs, which carry no Xa/Xc/Xn/Xt, so the
    four columns reach the merge empty. The align stage wrote those tags one
    step earlier onto `<sample>.multialigned.bam`; the merge joins them back by
    read name. Raw tags across a join — not a computed concordance class."""

    def _merge(self, tmp_path, per_aligner_rows, consensus_bam=None):
        import pandas as pd
        from rectify.core.consensus.corrected_consensus import merge_corrected_tsvs

        per_aligner_tsvs = {}
        for aligner, rows in per_aligner_rows.items():
            path = tmp_path / f"{aligner}.tsv"
            pd.DataFrame(rows, columns=_MERGE_COLUMNS).to_csv(
                path, sep='\t', index=False
            )
            per_aligner_tsvs[aligner] = path
        out_tsv = tmp_path / "merged.tsv"
        merge_corrected_tsvs(
            per_aligner_tsvs=per_aligner_tsvs,
            output_tsv=out_tsv,
            per_aligner_corrected_bams=None,
            consensus_bam=consensus_bam,
        )
        return pd.read_csv(out_tsv, sep='\t', keep_default_na=False)

    def _two_empty_arms(self, read_id='r1'):
        return {
            'minimap2': [_merge_row(read_id, 'chrI', 500, span=200)],
            'deSALT': [_merge_row(read_id, 'chrI', 500, span=100)],
        }

    def test_tagged_bam_fills_the_four_columns(self, tmp_path):
        bam = _write_consensus_bam(tmp_path / "s.multialigned.bam", [
            ('r1', {'Xa': 'minimap2', 'Xc': 'high', 'Xn': 3,
                    'Xt': 'deSALT,minimap2'}),
        ])
        merged = self._merge(tmp_path, self._two_empty_arms(), consensus_bam=bam)
        assert len(merged) == 1
        row = merged.iloc[0]
        assert row['consensus_aligner'] == 'minimap2'
        assert row['consensus_confidence'] == 'high'
        assert str(row['consensus_n_agree']) == '3'
        assert row['consensus_tied'] == 'deSALT,minimap2'
        # The merge's own verdict is untouched and still separate.
        assert row['winning_aligner'] in ('minimap2', 'deSALT')

    def test_no_bam_leaves_them_empty(self, tmp_path):
        merged = self._merge(tmp_path, self._two_empty_arms(), consensus_bam=None)
        for column in CONSENSUS_TAG_COLUMNS:
            assert merged.iloc[0][column] == ''

    def test_untagged_bam_leaves_them_empty(self, tmp_path):
        bam = _write_consensus_bam(tmp_path / "s.multialigned.bam", [('r1', {})])
        merged = self._merge(tmp_path, self._two_empty_arms(), consensus_bam=bam)
        for column in CONSENSUS_TAG_COLUMNS:
            assert merged.iloc[0][column] == ''

    def test_existing_values_are_not_overwritten(self, tmp_path):
        """A per-aligner TSV whose own BAM really did carry the tags is the
        more direct evidence; the join must not clobber it."""
        rows = {
            'minimap2': [_merge_row(
                'r1', 'chrI', 500, span=200,
                consensus_aligner='uLTRA', consensus_confidence='low',
                consensus_n_agree='1', consensus_tied='',
            )],
            'deSALT': [_merge_row(
                'r1', 'chrI', 500, span=100,
                consensus_aligner='uLTRA', consensus_confidence='low',
                consensus_n_agree='1', consensus_tied='',
            )],
        }
        bam = _write_consensus_bam(tmp_path / "s.multialigned.bam", [
            ('r1', {'Xa': 'minimap2', 'Xc': 'high', 'Xn': 3, 'Xt': 'a,b'}),
        ])
        merged = self._merge(tmp_path, rows, consensus_bam=bam)
        row = merged.iloc[0]
        assert row['consensus_aligner'] == 'uLTRA'
        assert row['consensus_confidence'] == 'low'
        assert str(row['consensus_n_agree']) == '1'
        # Only the genuinely empty cell is filled.
        assert row['consensus_tied'] == 'a,b'

    def test_existing_n_agree_survives_a_disagreeing_bam(self, tmp_path):
        """Sharper than the case above, which used '1' on both sides so int-vs-str
        dtype inference could not be told apart. An all-numeric column is read
        back from the arms as int64, so the no-overwrite path has to survive
        `.astype(str)` — assert a value the BAM actively disagrees with."""
        rows = {
            'minimap2': [_merge_row(
                'r1', 'chrI', 500, span=200,
                consensus_aligner='uLTRA', consensus_n_agree='3',
            )],
            'deSALT': [_merge_row(
                'r1', 'chrI', 500, span=100,
                consensus_aligner='uLTRA', consensus_n_agree='3',
            )],
        }
        bam = _write_consensus_bam(tmp_path / "s.multialigned.bam", [
            ('r1', {'Xa': 'minimap2', 'Xc': 'high', 'Xn': 5, 'Xt': 'a,b'}),
        ])
        merged = self._merge(tmp_path, rows, consensus_bam=bam)
        row = merged.iloc[0]
        assert str(row['consensus_n_agree']) == '3'   # NOT the BAM's 5
        assert row['consensus_aligner'] == 'uLTRA'    # NOT the BAM's minimap2
        # The two cells that really were empty do fill.
        assert row['consensus_confidence'] == 'high'
        assert row['consensus_tied'] == 'a,b'

    def test_read_absent_from_the_bam_stays_empty(self, tmp_path):
        bam = _write_consensus_bam(tmp_path / "s.multialigned.bam", [
            ('someone_else', {'Xa': 'minimap2', 'Xc': 'high', 'Xn': 2}),
        ])
        merged = self._merge(tmp_path, self._two_empty_arms(), consensus_bam=bam)
        for column in CONSENSUS_TAG_COLUMNS:
            assert merged.iloc[0][column] == ''

    def test_qname_suffix_is_normalized_before_the_join(self, tmp_path):
        """The join key must be `_normalize_bam_read_name` on BOTH sides —
        `_load_tsv` already normalizes the TSV `read_id`. An unnormalized BAM
        key would match nothing and fail silently as four empty columns."""
        bam = _write_consensus_bam(tmp_path / "s.multialigned.bam", [
            ('r1_pt:i:25', {'Xa': 'mapPacBio', 'Xc': 'medium', 'Xn': 2}),
        ])
        merged = self._merge(tmp_path, self._two_empty_arms(), consensus_bam=bam)
        assert merged.iloc[0]['consensus_aligner'] == 'mapPacBio'
        assert str(merged.iloc[0]['consensus_n_agree']) == '2'

    def test_columns_are_created_when_the_arms_predate_them(self, tmp_path):
        """Per-aligner TSVs written before this feature have no such columns."""
        import pandas as pd
        from rectify.core.consensus.corrected_consensus import merge_corrected_tsvs

        legacy_cols = [c for c in _MERGE_COLUMNS
                       if c not in CONSENSUS_TAG_COLUMNS]
        per_aligner_tsvs = {}
        for aligner, span in (('minimap2', 200), ('deSALT', 100)):
            path = tmp_path / f"{aligner}.tsv"
            pd.DataFrame([_merge_row('r1', 'chrI', 500, span=span)],
                         columns=legacy_cols).to_csv(path, sep='\t', index=False)
            per_aligner_tsvs[aligner] = path
        bam = _write_consensus_bam(tmp_path / "s.multialigned.bam", [
            ('r1', {'Xa': 'deSALT', 'Xc': 'high', 'Xn': 4}),
        ])
        out_tsv = tmp_path / "merged.tsv"
        merge_corrected_tsvs(
            per_aligner_tsvs=per_aligner_tsvs,
            output_tsv=out_tsv,
            per_aligner_corrected_bams=None,
            consensus_bam=bam,
        )
        merged = pd.read_csv(out_tsv, sep='\t', keep_default_na=False)
        assert merged.iloc[0]['consensus_aligner'] == 'deSALT'
        assert merged.iloc[0]['consensus_tied'] == ''

    def test_only_tagged_reads_are_retained_in_memory(self, tmp_path):
        """The BAM is streamed once into four short strings per TAGGED read;
        untagged reads cost nothing."""
        from rectify.core.consensus.corrected_consensus import (
            _load_consensus_tags_from_bam,
        )

        bam = _write_consensus_bam(tmp_path / "s.multialigned.bam", [
            ('tagged', {'Xa': 'minimap2', 'Xc': 'high', 'Xn': 3}),
            ('untagged', {}),
        ])
        tags = _load_consensus_tags_from_bam(bam)
        assert set(tags) == {'tagged'}
        assert tags['tagged'] == ('minimap2', 'high', '3', '')


class TestResolveConsensusTagBam:
    def _paths(self, tmp_path, name):
        path = tmp_path / name
        path.write_bytes(b'')  # existence is all the resolver checks
        return path

    def test_explicit_path_wins(self, tmp_path):
        from rectify.core.commands.run.helpers import _resolve_consensus_tag_bam

        explicit = self._paths(tmp_path, 'explicit.bam')
        self._paths(tmp_path, 's.multialigned.bam')
        assert _resolve_consensus_tag_bam(str(explicit), 's', (tmp_path,)) == explicit

    def test_missing_explicit_path_returns_none(self, tmp_path):
        from rectify.core.commands.run.helpers import _resolve_consensus_tag_bam

        self._paths(tmp_path, 's.multialigned.bam')
        assert _resolve_consensus_tag_bam(
            str(tmp_path / 'nope.bam'), 's', (tmp_path,)) is None

    def test_auto_discovers_multialigned(self, tmp_path):
        from rectify.core.commands.run.helpers import _resolve_consensus_tag_bam

        expected = self._paths(tmp_path, 's.multialigned.bam')
        assert _resolve_consensus_tag_bam(None, 's', (tmp_path,)) == expected

    def test_auto_discovers_legacy_consensus_name(self, tmp_path):
        from rectify.core.commands.run.helpers import _resolve_consensus_tag_bam

        expected = self._paths(tmp_path, 's.consensus.bam')
        assert _resolve_consensus_tag_bam(None, 's', (tmp_path,)) == expected

    def test_searches_later_dirs_and_skips_none(self, tmp_path):
        from rectify.core.commands.run.helpers import _resolve_consensus_tag_bam

        second = tmp_path / 'second'
        second.mkdir()
        expected = self._paths(second, 's.multialigned.bam')
        assert _resolve_consensus_tag_bam(
            None, 's', (tmp_path, None, second)) == expected

    def test_falls_back_to_the_corrected_input_bam(self, tmp_path):
        from rectify.core.commands.run.helpers import _resolve_consensus_tag_bam

        fallback = self._paths(tmp_path, 'input.bam')
        assert _resolve_consensus_tag_bam(
            None, 's', (tmp_path,), fallback=fallback) == fallback

    def test_returns_none_when_nothing_exists(self, tmp_path):
        from rectify.core.commands.run.helpers import _resolve_consensus_tag_bam

        assert _resolve_consensus_tag_bam(None, 's', (tmp_path,)) is None
