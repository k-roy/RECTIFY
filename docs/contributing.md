# Contributing

## Development setup

```bash
git clone https://github.com/k-roy/RECTIFY
cd RECTIFY
pip install -e ".[dev]"
```

This installs RECTIFY in editable mode with development dependencies: `pytest`, `pytest-cov`, `black`, `flake8`, `mypy`.

---

## Running tests

```bash
# Run all tests
pytest

# Run with coverage report
pytest --cov=rectify --cov-report=html

# Run a specific test file
pytest tests/test_indel_corrector.py -v

# Run tests matching a pattern
pytest -k "atract" -v
```

The test suite has 435 tests. All should pass before submitting a PR.

---

## Code style

RECTIFY uses `black` for formatting with a line length of 100:

```bash
black rectify/ tests/
```

Linting with `flake8`:

```bash
flake8 rectify/ tests/
```

---

## Architecture notes

### Two `run-all` implementations

There are two `run-all` implementations. **Do not confuse them.**

| File | Status | Purpose |
|------|--------|---------|
| `rectify/core/run_command.py` | **Active** (wired to CLI) | Production pipeline |
| `rectify/core/run_all_command.py` | Inactive | Experimental redesign |

Always edit `run_command.py` for pipeline changes.

### Coordinate convention

All internal coordinates are **0-based, half-open** (BED/pysam convention). GFF files use 1-based coordinates — subtract 1 on load.

See [Coordinate System](coordinate_system.md) for details.

### Memory-efficient streaming

The multi-sample analysis pipeline is designed to run on datasets of any size. Key constraints:
- Never load > 1 sample's data at once
- Position columns only; drop all other columns immediately
- Use position index files (`_index.bed.gz`) when available

See [Multi-Sample Analysis](user_guide/multi_sample.md) for details.

---

## Adding a new correction module

1. Create `rectify/core/my_module.py`
2. Add the correction function signature:
   ```python
   def my_correction(read: pysam.AlignedSegment, strand: str, genome: dict, **kwargs) -> Optional[dict]:
       """
       Correct position using my algorithm.

       Returns: dict with corrected_pos, confidence, or None if no correction
       """
   ```
3. Wire into `bam_processor.py` → `correct_read_3prime()`
4. Add tests in `tests/test_my_module.py`
5. Add a flag to `cli.py` if it should be user-controllable (e.g., `--skip-my-module`)

---

## Adding a new analysis module

1. Create `rectify/core/analyze/my_analysis.py`
2. Implement functions with standard input signatures:
   ```python
   def run_my_analysis(count_matrix: pd.DataFrame, metadata: pd.DataFrame, **kwargs) -> pd.DataFrame:
       """Run my analysis."""
   ```
3. Wire into `analyze_command.py`
4. Add a flag to `cli.py` (e.g., `--run-my-analysis`)
5. Add tests in `tests/test_my_analysis.py`

---

## Building the documentation locally

```bash
pip install "rectify-rna[docs]"
mkdocs serve    # Live reload at http://127.0.0.1:8000
mkdocs build    # Build static site in site/
```

---

## Submitting changes

1. Fork the repository on GitHub
2. Create a feature branch: `git checkout -b my-feature`
3. Make your changes and add tests
4. Ensure all tests pass: `pytest`
5. Format code: `black rectify/ tests/`
6. Open a pull request against `main`

---

## Reporting bugs

Please open an issue at [https://github.com/k-roy/RECTIFY/issues](https://github.com/k-roy/RECTIFY/issues) with:

- RECTIFY version (`rectify --version`)
- Python version and OS
- Minimal reproducible example
- Full traceback if applicable
