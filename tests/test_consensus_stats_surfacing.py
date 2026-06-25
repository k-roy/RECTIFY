"""Tests that consensus aligner stats (tied_score, by_aligner_combo) are
surfaced in both the stats TSV and the HTML report (dev/TODO.md)."""

from rectify.core.bam.processing_stats import write_consensus_stats_tsv
from rectify.core.analyze.summary import generate_consensus_html_report


def _sample_stats():
    return {
        'total_reads': 100,
        'consensus_high': 80,
        'consensus_medium': 15,
        'consensus_low': 5,
        '5prime_rescued': 7,
        'tied_score': 12,
        'by_aligner': {'minimap2': 60, 'gapmm2': 25, 'uLTRA': 15},
        'by_aligner_combo': {
            frozenset({'minimap2', 'gapmm2', 'uLTRA'}): 70,
            frozenset({'minimap2', 'gapmm2'}): 20,
            frozenset({'minimap2'}): 10,
        },
    }


def test_stats_tsv_includes_combo_and_tied(tmp_path):
    out = tmp_path / "stats.tsv"
    write_consensus_stats_tsv(_sample_stats(), str(out))
    text = out.read_text()
    assert "tied_score\t12\t12.00" in text
    # Combo rows, '+'-joined and sorted, descending by count.
    assert "aligner_combo_gapmm2+minimap2+uLTRA\t70\t70.00" in text
    assert "aligner_combo_gapmm2+minimap2\t20\t20.00" in text
    assert "aligner_combo_minimap2\t10\t10.00" in text
    # Descending order: the 3-aligner combo precedes the singleton.
    assert text.index("gapmm2+minimap2+uLTRA") < text.index("aligner_combo_minimap2\t")


def test_html_report_includes_combo_section(tmp_path):
    out = tmp_path / "report.html"
    generate_consensus_html_report(_sample_stats(), str(out), sample_name="s1")
    html = out.read_text()
    assert "Aligner Combination Breakdown" in html
    assert "gapmm2+minimap2+uLTRA" in html
    assert "Tiebreaker used" in html


def test_html_combo_section_absent_when_no_combo_data(tmp_path):
    stats = _sample_stats()
    del stats['by_aligner_combo']
    out = tmp_path / "report.html"
    generate_consensus_html_report(stats, str(out), sample_name="s1")
    html = out.read_text()
    assert "Aligner Combination Breakdown" not in html
    # Core report still renders.
    assert "Aligner Win Distribution" in html


def test_stats_tsv_no_combo_key_is_safe(tmp_path):
    stats = _sample_stats()
    del stats['by_aligner_combo']
    out = tmp_path / "stats.tsv"
    write_consensus_stats_tsv(stats, str(out))
    text = out.read_text()
    assert "aligner_combo_" not in text
    assert "tied_score" in text
