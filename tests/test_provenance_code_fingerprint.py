#!/usr/bin/env python3
"""The provenance sidecar must actually pin the code that produced an output.

Regression guard for a silent, total failure found 2026-08-02: ``__git_sha__``
was READ by ``sidecar.py`` but never SET anywhere in the codebase, so every
provenance record ever written — on every machine — said
``rectify_git_sha: "unknown"``. The field designed to pin the version was inert,
which is how a ~2-month-stale rectify on Hoffman2 (a plain non-git copy of the
repo) produced the MS2 3'-end tables with nothing flagging the drift.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import rectify.core.provenance.sidecar as sc


def _reset():
    sc._CODE_FINGERPRINT = None


def test_fingerprint_is_never_the_useless_unknown():
    """The whole point: it must resolve to something that identifies the code."""
    _reset()
    fp = sc.resolve_code_fingerprint()
    _reset()
    assert fp != "unknown"
    assert fp.startswith(("git:", "pkghash:")), fp


def test_fingerprint_is_cached():
    _reset()
    a = sc.resolve_code_fingerprint()
    b = sc.resolve_code_fingerprint()
    _reset()
    assert a == b


def test_injected_sha_wins(monkeypatch):
    import rectify
    _reset()
    monkeypatch.setattr(rectify, "__git_sha__", "injected-abc123", raising=False)
    assert sc.resolve_code_fingerprint() == "injected-abc123"
    _reset()


def test_pkghash_changes_when_any_py_file_changes(tmp_path):
    """The case that actually bit us: a deployment that is a COPY of the repo,
    with no .git, must still get a fingerprint that MOVES when the code moves.

    Mirrors the helper's fallback (3) exactly; if that formula is changed, this
    test should be updated deliberately rather than silently drifting.
    """
    def pkghash(pkg: Path) -> str:
        h = hashlib.sha256()
        for f in sorted(pkg.rglob("*.py")):
            h.update(str(f.relative_to(pkg)).encode())
            h.update(f.read_bytes())
        return "pkghash:" + h.hexdigest()[:12]

    pkg = tmp_path / "rectify"
    (pkg / "core").mkdir(parents=True)
    (pkg / "__init__.py").write_text("__version__ = '0.9.0'\n")
    (pkg / "core" / "mod.py").write_text("X = 1\n")

    before = pkghash(pkg)
    (pkg / "core" / "mod.py").write_text("X = 2\n")   # one-character change
    after = pkghash(pkg)

    assert before != after, "a non-git deployment must still detect code drift"
    assert before.startswith("pkghash:") and len(before) == len("pkghash:") + 12


def test_record_carries_the_fingerprint(tmp_path):
    """End-to-end: a sidecar built without an explicit sha must not say unknown."""
    _reset()
    bam = tmp_path / "in.bam"
    bam.write_bytes(b"x")
    out = tmp_path / "out.tsv"
    out.write_text("a\n")
    rec = sc.ProvenanceRecord.from_components(
        stage="correct",
        sample_id="s1",
        sample_output_dir=tmp_path,
        inputs={"bam": bam},
        outputs={"tsv": out},
        started_at="2026-08-02T00:00:00Z",
        completed_at="2026-08-02T00:00:01Z",
        exit_status=0,
        argv=["rectify", "correct", str(bam)],
    )
    _reset()
    assert rec.rectify_git_sha != "unknown"
    assert rec.rectify_git_sha.startswith(("git:", "pkghash:"))


def test_trim_cdna_polya_writes_a_sidecar_with_a_real_sha(tmp_path, monkeypatch):
    """`trim-cdna-polya` rewrites the reads, so it MUST pin its code version.

    It wrote no provenance at all until 2026-08-02 — which is how the FASTQs
    behind the MS2 cDNA arm ended up with no record of the trim code that made
    them, in the very stage where two bugs were later found.
    """
    import argparse
    import gzip
    import json
    from rectify.core.commands import cdna_trim_command as ct

    # PortablePath rightly refuses to record non-transient paths under $TMPDIR
    # (ephemeral). pytest's tmp_path lives there, so root this run in $SCRATCH,
    # which is a stable env root and wins on longest-prefix match.
    monkeypatch.delenv("TMPDIR", raising=False)
    monkeypatch.delenv("L_SCRATCH", raising=False)
    monkeypatch.setenv("SCRATCH", str(tmp_path.resolve()))

    fq = tmp_path / "s.fastq.gz"
    with gzip.open(fq, "wt") as fh:
        for i in range(20):
            fh.write("@r%d\n%s\n+\n%s\n" % (i, "ACGT" * 30 + "A" * 20, "I" * 140))
    out = tmp_path / "trim"

    args = argparse.Namespace(
        input_fastq=str(fq), output_dir=str(out), prefix="s", tsv=True,
        adapter_3p=ct._CRTA, adapter_5p=ct._CRTA_RC,
        seed_len=ct._SEED_LEN, seed_max_mismatches=ct._SEED_MAX_MISMATCHES,
        max_error_rate=0.0, max_consecutive_non_a=1,
        adapter_window=ct._ADAPTER_WINDOW, trim_5p_polyt=True,
        min_polya=1, min_polyt=1,
    )
    assert ct.run(args) == 0

    sidecars = list(out.rglob("*provenance*.json"))
    assert sidecars, "trim-cdna-polya wrote NO provenance sidecar"
    d = json.loads(sidecars[0].read_text())
    assert d["stage"] == "trim_cdna_polya"
    assert d["rectify_git_sha"] != "unknown", "sidecar did not pin the code version"
    assert d["rectify_git_sha"].startswith(("git:", "pkghash:", "injected"))
