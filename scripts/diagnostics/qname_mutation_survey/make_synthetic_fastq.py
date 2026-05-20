#!/usr/bin/env python3
"""
Generate the synthetic adversarial FASTQ for the QNAME mutation survey.

Each read tests one specific QNAME edge case using a 201bp sequence from
chrI:10000-10200 of S288C R64-5-1 that maps with minimap2 -ax map-ont.

Usage:
    python3 make_synthetic_fastq.py [outdir]

Writes:
    <outdir>/synthetic_safe.fastq  -- all "safe" test cases in one file
    <outdir>/synthetic_adversarial.fastq -- potentially crashing cases separately
    <outdir>/test_cases.tsv  -- mapping from read_id -> description
"""
import os
import sys
import string

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "/Users/kevinroy/work/rectify/scripts/diagnostics/qname_mutation_survey"

# 201bp chrI:10000-10200 (S288C R64-5-1), maps with minimap2 -ax map-ont
SEQ = "AGAATGAATCGGATATACTATTTTTGATCATAATGACGGACATCATGATATAATAACGTTATACGGATAACTTTATTTCAAAAGCACCATCATGTTATCTCTTGTAAAAAGAAGTATTCTTCATTCAATACCAATTACTCGTCACATTCTTCCAATCCAATTAATATTGGTTAAAATGAACCATGTGCAAATCAGAAACAT"
QUAL = "I" * len(SEQ)

# Safe test cases: (read_id, fastq_header_comment, description)
# fastq_header = "@" + read_id [+ " " + comment]
SAFE_CASES = [
    # --- Whitespace axis ---
    ("bare_id", "", "Bare QNAME, no comment"),
    ("id_trailing_space", " ", "Trailing space before newline"),
    ("id_space_comment", " some plain comment text", "Space + plain comment"),
    ("id_dorado_style", " runid=abc ch=42 start_time=2024 flow_cell_id=XYZ", "ONT Dorado-style space-separated header"),
    ("id_multi_space", "  one  two  three", "Multiple spaces in comment"),
    # tab-separated tags (cDNA-pipeline shape)
    ("id_tab_XA", "\tXA:Z:foo", "Tab-separated single tag"),
    ("id_tab_multi", "\tXA:Z:foo\tXC:i:42\tXF:Z:bar", "Tab-separated multiple tags"),

    # --- Mate / Illumina suffixes ---
    ("id_mate_slash_1", "/1", "Slash-1 mate suffix in QNAME body"),
    ("id_mate_slash_2", "/2", "Slash-2 mate suffix in QNAME body"),
    ("id_illumina_1", " 1:N:0:ATCGAT", "Illumina read 1 comment"),
    ("id_illumina_2", " 2:N:0:ATCGAT", "Illumina read 2 comment"),

    # --- Special chars in QNAME body (no space or tab in QNAME body itself) ---
    ("id_colon:char", "", "Colon in QNAME"),
    ("id_slash/char", "", "Slash in QNAME"),
    ("id_pipe|char", "", "Pipe in QNAME"),
    ("id_dash-char", "", "Dash in QNAME"),
    ("id_underscore_char", "", "Underscore in QNAME"),
    ("id_dot.char", "", "Dot in QNAME"),
    ("id_hash#char", "", "Hash in QNAME"),
    ("id_equals=char", "", "Equals in QNAME"),
    ("id_plus+char", "", "Plus in QNAME"),
    ("id_comma,char", "", "Comma in QNAME"),
    ("id_semicolon;char", "", "Semicolon in QNAME"),
    ("id_star*char", "", "Star/asterisk in QNAME"),
    ("id_opensquare[char", "", "Open square bracket in QNAME"),
    ("id_closesquare]char", "", "Close square bracket in QNAME"),
    ("id_openparen(char", "", "Open paren in QNAME"),
    ("id_closeparen)char", "", "Close paren in QNAME"),
    ("id_opencurly{char", "", "Open curly brace in QNAME"),
    ("id_closecurly}char", "", "Close curly brace in QNAME"),
    ("id_lt<char", "", "Less-than in QNAME"),
    ("id_gt>char", "", "Greater-than in QNAME"),
    ("id_question?char", "", "Question mark in QNAME"),
    ("id_exclaim!char", "", "Exclamation mark in QNAME"),
    ("id_tilde~char", "", "Tilde in QNAME"),
    ("id_caret^char", "", "Caret in QNAME"),
    ("id_ampersand&char", "", "Ampersand in QNAME"),
    ("id_percent%char", "", "Percent in QNAME"),
    ("id_dollar$char", "", "Dollar in QNAME"),
    ("id_singlequote'char", "", "Single quote in QNAME"),
    ("id_doublequote\"char", "", "Double quote in QNAME"),
    ("id_backtick`char", "", "Backtick in QNAME"),

    # --- Length axis ---
    ("a", "", "1-char QNAME"),
    ("01234567-89ab-cdef-0123-456789abcdef", "", "36-char UUID (typical Dorado)"),
    ("A" * 100, "", "100-char QNAME"),
    ("B" * 254, "", "254-char QNAME (SAM spec limit)"),

    # --- Numeric/edge-form ---
    ("1234567890", "", "Numeric-only QNAME"),
    ("00000123", "", "Leading-zeros QNAME"),
    ("-id_dash_lead", "", "Leading-dash QNAME (CLI hazard)"),
    ("-42", "", "Negative-looking QNAME"),

    # --- UUID forms ---
    ("01234567-89AB-CDEF-0123-456789ABCDEF", "", "Uppercase UUID"),

    # --- Duplicate QNAME (second occurrence, different QNAME suffix) ---
    ("dup_id", "", "First occurrence of duplicate QNAME"),
    ("dup_id", " duplicate_occurrence", "Second occurrence of duplicate QNAME (different comment)"),

    # --- QNAME with @ symbol ---
    ("id_with_at@symbol", "", "QNAME containing literal @"),
]

# Adversarial/potentially crashing cases — kept separate
ADVERSARIAL_CASES = [
    # Length over SAM spec limit
    ("C" * 256, "", "256-char QNAME (over SAM 254-char limit)"),
    ("D" * 1024, "", "1024-char QNAME (way over limit)"),
    # NUL byte — most FASTQ parsers will reject or truncate
    ("id_nul\x00byte", "", "QNAME with NUL byte (\\x00)"),
]


def write_fastq_record(fh, qname, comment, seq, qual):
    """Write one FASTQ record. header = @<qname>[<comment>]"""
    header = "@" + qname + comment
    fh.write(f"{header}\n{seq}\n+\n{qual}\n")


def main():
    safe_path = os.path.join(OUTDIR, "synthetic_safe.fastq")
    adv_path = os.path.join(OUTDIR, "synthetic_adversarial.fastq")
    cases_path = os.path.join(OUTDIR, "test_cases.tsv")

    with open(safe_path, "w") as safe_fh, open(adv_path, "w") as adv_fh, open(cases_path, "w") as tsv_fh:
        tsv_fh.write("case_num\tread_id\tcomment_raw\tdescription\tfastq_file\n")
        n = 0
        for qname, comment, desc in SAFE_CASES:
            n += 1
            write_fastq_record(safe_fh, qname, comment, SEQ, QUAL)
            # Replace tabs in comment repr for TSV display
            comment_repr = comment.replace("\t", "\\t")
            tsv_fh.write(f"{n}\t{qname}\t{comment_repr!r}\t{desc}\tsafe\n")
        for qname, comment, desc in ADVERSARIAL_CASES:
            n += 1
            write_fastq_record(adv_fh, qname, comment, SEQ, QUAL)
            comment_repr = comment.replace("\t", "\\t")
            # Truncate display of long qnames
            display_qname = qname if len(qname) <= 40 else f"{qname[:20]}...[{len(qname)} chars]"
            tsv_fh.write(f"{n}\t{display_qname}\t{comment_repr!r}\t{desc}\tadversarial\n")

    print(f"Wrote {len(SAFE_CASES)} safe records → {safe_path}")
    print(f"Wrote {len(ADVERSARIAL_CASES)} adversarial records → {adv_path}")
    print(f"Wrote test case index → {cases_path}")

if __name__ == "__main__":
    main()
