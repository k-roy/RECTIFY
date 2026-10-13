"""`rectify analyze` draws a logo of each scan's top motif (the deliverable spec asks for plots of the motifs).

Until 2026-10-13 the motif step left only the STREME/MEME HTML reports (gap ii of the 2026-09-11 inventory).
"""

import pytest

pytest.importorskip("matplotlib")

from rectify.core.analyze.motif_discovery import parse_meme_minimal, plot_motif_logos  # noqa: E402

STREME_TXT = """MEME version 5.5.4

ALPHABET= ACGT

strands: + -

MOTIF 1-TATATA STREME-1
letter-probability matrix: alength= 4 w= 6 nsites= 40 E= 2.5e-008
 0.05 0.05 0.05 0.85
 0.85 0.05 0.05 0.05
 0.05 0.05 0.05 0.85
 0.85 0.05 0.05 0.05
 0.05 0.05 0.05 0.85
 0.85 0.05 0.05 0.05

MOTIF 2-GGG STREME-2
letter-probability matrix: alength= 4 w= 3 nsites= 12 S= 1.0e-002
 0.1 0.1 0.7 0.1
 0.1 0.1 0.7 0.1
 0.1 0.1 0.7 0.1
"""


def test_parser_reads_every_motif_in_order(tmp_path):
    p = tmp_path / "streme.txt"
    p.write_text(STREME_TXT)
    motifs = parse_meme_minimal(str(p))
    assert [m["alt"] for m in motifs] == ["STREME-1", "STREME-2"]
    assert len(motifs[0]["pwm"]) == 6 and motifs[0]["nsites"] == 40
    assert motifs[0]["evalue"] == pytest.approx(2.5e-8)
    assert motifs[1]["evalue"] == pytest.approx(1.0e-2)


def test_logo_png_for_the_scans_that_found_motifs(tmp_path):
    d = tmp_path / "enriched_upstream_streme_results"
    d.mkdir()
    (d / "streme.txt").write_text(STREME_TXT)
    results = {"enriched_upstream": {"success": True, "output_dir": str(d)},
               "depleted_upstream": {"success": False, "error": "skipped"}}
    out = plot_motif_logos(results, str(tmp_path / "plots" / "motif_logos_x.png"), "x")
    assert out and (tmp_path / "plots" / "motif_logos_x.png").stat().st_size > 1000


def test_no_motifs_writes_nothing(tmp_path):
    out = plot_motif_logos({"enriched_upstream": {"success": False}}, str(tmp_path / "none.png"))
    assert out is None and not (tmp_path / "none.png").exists()
