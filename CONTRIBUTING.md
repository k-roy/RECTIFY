# Contributing to RECTIFY

Thanks for your interest in improving RECTIFY. This file is the entry point
for contributors; the detailed engineering walkthrough lives at
[docs/contributing.md](docs/contributing.md).

## Quick checklist

1. **Open an issue first** for anything larger than a typo or a small bug fix,
   so we can agree on scope before you write code.
2. **Fork → branch → PR.** Branch from `master` and use a descriptive
   slash-separated name, e.g. `fix/walkback-genomic-A-skip` or
   `feat/multisample-strand-export`.
3. **Run the test suite** (`pytest -m "not slow"`, ~1,640 tests, ~1 minute) and
   confirm it is green before pushing. Add tests for any new behavior.
4. **Match the coding style** — `black` (line length 100) and `flake8`. The
   `dev` extras (`pip install -e ".[dev]"`) install both.
5. **Update [CHANGELOG.md](CHANGELOG.md)** under `## [Unreleased]` describing
   user-visible changes (Added / Changed / Fixed / Removed / Deprecated).
6. **Sign your commits** with a descriptive message — see existing commits for
   tone and granularity.

## Development environment

```bash
git clone https://github.com/k-roy/RECTIFY
cd RECTIFY
pip install -e ".[dev,visualize]"
pytest -m "not slow"
```

The fast suite (~1,640 tests, ~1 minute) runs without optional binary aligners.
Slow integration tests (marked `@pytest.mark.slow`) require `minimap2` and
`mapPacBio` on `PATH`; install via conda:

```bash
conda install -c bioconda minimap2 bbmap
```

## Reporting bugs

Open a GitHub issue with:

- RECTIFY version (`rectify --version`)
- Python version, OS
- Minimal command that reproduces the bug
- Full traceback or unexpected output
- Sample BAM / FASTQ (≤10 MB) if relevant, or a public dataset accession

For sensitive data or security issues, see [SECURITY.md](SECURITY.md).

## Pull request review

Reviews typically respond within a week. The reviewer will check:

- Tests cover the new behavior and the suite is green
- Public-facing changes are documented in CHANGELOG and (where relevant) under
  `docs/`
- No regression in the bundled-data smoke check
  (`rectify test`)
- Coding style passes `black` and `flake8`

Once approved, the reviewer squash-merges the PR.

## Code of conduct

Participation in this project is governed by our
[CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md). Contributors agree to follow it.

## Detailed contributor guide

The full contributor guide — architecture, module boundaries, test layout,
performance pitfalls, the bundled validation set, and how to add a new
aligner — is in [docs/contributing.md](docs/contributing.md).
