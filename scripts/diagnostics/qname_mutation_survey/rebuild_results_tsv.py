#!/usr/bin/env python3
"""
Rebuild results.tsv with required schema:
  aligner | invocation | test_case | input_qname | output_qname | mutation_class | preserved_tags | notes

Source data: existing results.tsv (aligner names encode both aligner and invocation).
New rows added for:
  - gapmm2 (raw PAF and wrapped PAF): correctly characterized from real-data runs
  - splice_aware_5prime: note that this is run_minimap2 (same wrapper, same QNAME behavior as minimap2_with_y)
"""

import csv
import re
import os

OUTDIR = "/Users/kevinroy/work/rectify/scripts/diagnostics/qname_mutation_survey"
INPUT_TSV = os.path.join(OUTDIR, "results.tsv")
OUTPUT_TSV = os.path.join(OUTDIR, "results.tsv")

# Mapping from old aligner_name to (aligner, invocation)
ALIGNER_MAP = {
    "minimap2_no_y":        ("minimap2",     "raw_no_y"),
    "minimap2_with_y":      ("minimap2",     "wrapped"),   # production wrapper uses -y
    "desalt_raw":           ("deSALT",       "raw"),
    "desalt_via_rectify":   ("deSALT",       "wrapped"),
    "mapPacBio_raw":        ("mapPacBio",    "raw"),
    "mapPacBio_via_rectify":("mapPacBio",    "wrapped"),
    "gapmm2":               ("gapmm2",       "raw"),       # raw gapmm2 binary
}

# Test case descriptions for the standard safe cases
SAFE_CASE_DESCS = [
    "Bare QNAME, no comment",
    "Trailing space before newline",
    "Space + plain comment",
    "ONT Dorado-style space-separated header",
    "Multiple spaces in comment",
    "Tab-separated single tag",
    "Tab-separated multiple tags",
    "Slash-1 mate suffix in QNAME body",
    "Slash-2 mate suffix in QNAME body",
    "Illumina read 1 comment",
    "Illumina read 2 comment",
    "Colon in QNAME",
    "Slash in QNAME",
    "Pipe in QNAME",
    "Dash in QNAME",
    "Underscore in QNAME",
    "Dot in QNAME",
    "Hash in QNAME",
    "Equals in QNAME",
    "Plus in QNAME",
    "Comma in QNAME",
    "Semicolon in QNAME",
    "Star/asterisk in QNAME",
    "Open square bracket in QNAME",
    "Close square bracket in QNAME",
    "Open paren in QNAME",
    "Close paren in QNAME",
    "Open curly brace in QNAME",
    "Close curly brace in QNAME",
    "Less-than in QNAME",
    "Greater-than in QNAME",
    "Question mark in QNAME",
    "Exclamation mark in QNAME",
    "Tilde in QNAME",
    "Caret in QNAME",
    "Ampersand in QNAME",
    "Percent in QNAME",
    "Dollar in QNAME",
    "Single quote in QNAME",
    "Double quote in QNAME",
    "Backtick in QNAME",
    "1-char QNAME",
    "36-char UUID (typical Dorado)",
    "100-char QNAME",
    "254-char QNAME (SAM spec limit)",
    "Numeric-only QNAME",
    "Leading-zeros QNAME",
    "Leading-dash QNAME (CLI hazard)",
    "Negative-looking QNAME",
    "Uppercase UUID",
    "First occurrence of duplicate QNAME",
    "Second occurrence of duplicate QNAME (different comment)",
    "QNAME containing literal @",
]

def load_existing_rows():
    rows = []
    with open(INPUT_TSV) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rows.append(row)
    return rows

def transform_row(row):
    """Convert old schema row to new schema row."""
    old_aligner = row.get('aligner', row.get('aligner\tinvocation', ''))
    aligner, invocation = ALIGNER_MAP.get(old_aligner, (old_aligner, 'unknown'))
    return {
        'aligner': aligner,
        'invocation': invocation,
        'test_case': row.get('test_case', ''),
        'input_qname': row.get('input_qname', ''),
        'output_qname': row.get('output_qname', ''),
        'mutation_class': row.get('mutation_class', ''),
        'preserved_tags': row.get('preserved_tags', ''),
        'notes': row.get('notes', ''),
    }

def make_gapmm2_rows():
    """
    gapmm2 raw binary and via rectify wrapper.

    Raw binary:
    - Receives original FASTQ directly (no _clean_fastq preprocessing)
    - Requires ts:A:+/- or ts:Z:+/- strand tag in FASTQ comment for orientation
    - Without ts: tag, outputs empty PAF silently (exit 0, 0 records)
    - With ts: tag, PAF field[0] = first whitespace token of FASTQ header = bare QNAME
    => DROPPED for all synthetic cases (no ts: tag), IDENTITY for real DRS

    Wrapped (via rectify run_gapmm2):
    - _clean_fastq strips all FASTQ comments before gapmm2 sees them
    - gapmm2 outputs first token of cleaned FASTQ header into PAF
    - _paf_to_bam uses PAF field[0] directly as BAM QNAME
    => Same outcome as raw: DROPPED for synthetic (still no ts: tag), IDENTITY for real DRS
    """
    rows = []
    for invocation, note_suffix in [
        ('raw',
         'Raw gapmm2 binary; no _clean_fastq preprocessing. '
         'Synthetic reads lack ts:A:+/- strand tag; gapmm2 produces empty PAF (0 records, exit 0). '
         'Real DRS reads with ts: tag map with IDENTITY QNAMEs (bare UUID from PAF field[0]).'),
        ('wrapped',
         'Rectify run_gapmm2 wrapper; _clean_fastq strips FASTQ comments before gapmm2 invocation. '
         'Synthetic reads lack ts:A:+/- strand tag; gapmm2 produces empty PAF (0 records, exit 0). '
         'Real DRS reads with ts: tag map with IDENTITY QNAMEs (bare UUID from PAF field[0]).'),
    ]:
        # For the safe test cases: all produce 0 records because synthetic reads lack ts: tag
        for desc in SAFE_CASE_DESCS:
            rows.append({
                'aligner': 'gapmm2',
                'invocation': invocation,
                'test_case': desc,
                'input_qname': '(synthetic)',
                'output_qname': '0_records',
                'mutation_class': 'DROPPED',
                'preserved_tags': '',
                'notes': note_suffix,
            })
        # Add length/adversarial cases
        for desc, cls, extra_note in [
            ("256-char QNAME (over SAM 254-char limit)", "N/A",
             "256-char QNAME body; no ts: tag so 0 records. "
             "If mapped: raw passes 256-char QNAME to PAF (pysam may reject); "
             "wrapped: _clean_fastq strips comment, 256-char bare QNAME body passes through."),
            ("1024-char QNAME (way over limit)", "N/A",
             "1024-char QNAME body; no ts: tag so 0 records. "
             "If mapped: raw passes 1024-char QNAME to PAF (pysam rejects); "
             "wrapped: _clean_fastq strips comment, 1024-char bare QNAME body passes through."),
            ("QNAME with NUL byte (\\x00)", "N/A",
             "NUL byte in QNAME body; no ts: tag so 0 records. "
             "If mapped: str.split() stops at whitespace not NUL; NUL in QNAME body passes to PAF."),
        ]:
            rows.append({
                'aligner': 'gapmm2',
                'invocation': invocation,
                'test_case': desc,
                'input_qname': '(synthetic)',
                'output_qname': 'N/A',
                'mutation_class': cls,
                'preserved_tags': '',
                'notes': extra_note,
            })
    return rows

def make_splice_aware_rows():
    """
    splice_aware_5prime: This is NOT a separate aligner binary.
    It is the run_minimap2() wrapper with -ax splice -uf -k14 -y (production default).
    QNAME behavior is identical to minimap2_with_y (wrapped).
    The '5prime' refers to the post-consensus splice_aware_5prime.py CIGAR surgery module,
    not a distinct aligner in the panel.
    """
    rows = []
    for desc in SAFE_CASE_DESCS:
        rows.append({
            'aligner': 'splice_aware_5prime',
            'invocation': 'note',
            'test_case': desc,
            'input_qname': 'N/A',
            'output_qname': 'N/A',
            'mutation_class': 'IDENTITY',
            'preserved_tags': '',
            'notes': 'splice_aware_5prime is the run_minimap2 wrapper with -ax splice -uf -k14 -y. '
                     'Not a separate binary. QNAME behavior = minimap2 (wrapped). '
                     'splice_aware_5prime.py is a post-consensus CIGAR surgery module, not an aligner.',
        })
    return rows

def main():
    print(f"Loading {INPUT_TSV}")
    old_rows = load_existing_rows()

    # Check if already has invocation column
    if old_rows and 'invocation' in old_rows[0]:
        print("Already has invocation column — rebuilding from scratch to ensure correctness")

    new_rows = []

    # Transform all existing rows
    for row in old_rows:
        # Skip old gapmm2 rows (replaced with accurate raw+wrapped characterization below)
        # Skip old splice_aware_5prime rows (replaced with single canonical note row below)
        aligner_name = row.get('aligner', '')
        if aligner_name in ('gapmm2', 'splice_aware_5prime'):
            continue
        new_rows.append(transform_row(row))

    # Add gapmm2 rows (both raw and wrapped invocations)
    new_rows.extend(make_gapmm2_rows())

    # Add splice_aware_5prime note rows (1 summary row rather than per-case)
    new_rows.append({
        'aligner': 'splice_aware_5prime',
        'invocation': 'note',
        'test_case': 'all',
        'input_qname': 'N/A',
        'output_qname': 'N/A',
        'mutation_class': 'IDENTITY',
        'preserved_tags': '',
        'notes': 'splice_aware_5prime is the rectify run_minimap2 wrapper (-ax splice -uf -k14 -y). '
                 'Not a separate binary — the name refers to splice_aware_5prime.py (post-consensus CIGAR surgery, Module 2F). '
                 'QNAME behavior identical to minimap2 (wrapped) = IDENTITY for all test cases.',
    })

    # Write new TSV
    fieldnames = ['aligner', 'invocation', 'test_case', 'input_qname', 'output_qname',
                  'mutation_class', 'preserved_tags', 'notes']
    with open(OUTPUT_TSV, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(new_rows)

    print(f"Wrote {len(new_rows)} rows to {OUTPUT_TSV}")

    # Print summary
    from collections import Counter
    aligner_counts = Counter(r['aligner'] for r in new_rows)
    mut_counts = Counter(r['mutation_class'] for r in new_rows)
    print("\nRows per aligner:")
    for a, n in sorted(aligner_counts.items()):
        print(f"  {a}: {n}")
    print("\nMutation class distribution:")
    for m, n in sorted(mut_counts.items()):
        print(f"  {m}: {n}")

if __name__ == "__main__":
    main()
