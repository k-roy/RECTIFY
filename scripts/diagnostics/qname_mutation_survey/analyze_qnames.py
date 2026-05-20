#!/usr/bin/env python3
"""
Analyze QNAME mutations between synthetic FASTQ input and aligner BAM output.

For each aligner BAM, compare the input FASTQ QNAME (full header after @)
to the BAM QNAME (samtools view field 1). Classify the mutation.

Usage:
    python3 analyze_qnames.py

Writes:
    results.tsv  -- characterization table
"""
import os
import re
import subprocess
import sys

OUTDIR = "/Users/kevinroy/work/rectify/scripts/diagnostics/qname_mutation_survey"
FASTQ_SAFE = os.path.join(OUTDIR, "synthetic_safe.fastq")
FASTQ_ADV = os.path.join(OUTDIR, "synthetic_adversarial.fastq")

# Aligner BAMs to analyze (must already exist)
ALIGNER_BAMS = {
    "minimap2_no_y": os.path.join(OUTDIR, "minimap2_no_y.bam"),
    "minimap2_with_y": os.path.join(OUTDIR, "minimap2_with_y.bam"),
    "gapmm2": os.path.join(OUTDIR, "gapmm2.bam"),
    # H2 bams (added later)
}


def read_fastq_records(fastq_path):
    """Return list of (header_full, qname_first_tok, seq, qual)."""
    records = []
    try:
        with open(fastq_path) as f:
            while True:
                h = f.readline()
                if not h:
                    break
                s = f.readline()
                p = f.readline()
                q = f.readline()
                header_full = h[1:].rstrip('\n')
                # FASTQ: QNAME = everything up to first whitespace (space or tab)
                toks = re.split(r'[ \t]', header_full, maxsplit=1)
                qname_tok = toks[0]
                comment = toks[1] if len(toks) > 1 else ""
                records.append({
                    "header_full": header_full,
                    "input_qname": qname_tok,
                    "comment": comment,
                })
    except FileNotFoundError:
        pass
    return records


def read_bam_qnames(bam_path):
    """Return dict {qname: [aux_tags_string, ...]} from samtools view."""
    result = {}
    if not os.path.exists(bam_path):
        return result
    proc = subprocess.run(
        ["samtools", "view", bam_path],
        capture_output=True, text=True
    )
    for line in proc.stdout.splitlines():
        fields = line.split('\t')
        if not fields:
            continue
        qname = fields[0]
        # aux tags start at field 11 (0-indexed)
        aux = '\t'.join(fields[11:]) if len(fields) > 11 else ""
        result.setdefault(qname, []).append(aux)
    return result


def classify_mutation(input_qname, output_qname):
    """Classify the mutation type."""
    if output_qname is None:
        return "DROPPED"
    if output_qname == input_qname:
        return "IDENTITY"
    # Check if output is a prefix of input (truncated)
    if input_qname.startswith(output_qname):
        return "TRUNCATED_WHITESPACE"
    # Check if output is exactly 254 chars (length truncation)
    if len(output_qname) == 254 and input_qname.startswith(output_qname):
        return "TRUNCATED_LENGTH"
    # Check if a suffix was added
    if output_qname.startswith(input_qname):
        return "SUFFIX_ADDED"
    # Check for underscore encoding of spaces
    if output_qname == input_qname.replace(' ', '_'):
        return "MUTATED_CHARS"
    # Generic mutation
    return "OTHER"


def get_preserved_tags(bam_aux_str, comment):
    """Which FASTQ comment tags appear as BAM aux tags."""
    if not comment or not bam_aux_str:
        return ""
    # Look for SAM aux tag pattern XX:T:value in aux
    # Extract tag keys from comment
    preserved = []
    for m in re.finditer(r'([A-Za-z][A-Za-z0-9]):[A-Z]:(\S+)', comment):
        tag_key = m.group(1)
        tag_val = m.group(2)
        # Check if tag key appears in bam aux
        if re.search(rf'\b{re.escape(tag_key)}:[A-Z]:', bam_aux_str):
            preserved.append(tag_key)
    return ','.join(preserved) if preserved else ""


def main():
    # Load test cases from test_cases.tsv
    cases_path = os.path.join(OUTDIR, "test_cases.tsv")
    test_cases = []
    with open(cases_path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.rstrip('\n').split('\t')
            test_cases.append({
                "case_num": parts[0],
                "read_id": parts[1],
                "comment_raw": parts[2],
                "description": parts[3],
                "fastq_file": parts[4],
            })

    # Load FASTQ records
    safe_records = read_fastq_records(FASTQ_SAFE)
    adv_records = read_fastq_records(FASTQ_ADV)
    all_records = safe_records + adv_records

    # Verify count
    assert len(all_records) == len(test_cases), \
        f"Record count mismatch: {len(all_records)} records, {len(test_cases)} test cases"

    # Results TSV
    results_path = os.path.join(OUTDIR, "results.tsv")
    rows = []
    header = "aligner\ttest_case\tinput_qname\toutput_qname\tmutation_class\tpreserved_tags\tnotes"

    for aligner_name, bam_path in ALIGNER_BAMS.items():
        if not os.path.exists(bam_path):
            # Record as unavailable
            for i, rec in enumerate(safe_records):
                case = test_cases[i]
                rows.append({
                    "aligner": aligner_name,
                    "test_case": case["description"],
                    "input_qname": rec["input_qname"][:80],
                    "output_qname": "N/A",
                    "mutation_class": "CRASHED",
                    "preserved_tags": "",
                    "notes": f"BAM not found: {bam_path}",
                })
            continue

        bam_qnames = read_bam_qnames(bam_path)

        for i, rec in enumerate(safe_records):
            case = test_cases[i]
            input_q = rec["input_qname"]
            comment = rec["comment"]

            # Find matching BAM record(s)
            if input_q in bam_qnames:
                output_q = input_q
                aux_str = ' '.join(bam_qnames[input_q])
            else:
                # Try to find by other means (e.g., truncated)
                # Check if any BAM QNAME is a prefix match of the input
                found = None
                for bq in bam_qnames:
                    if input_q.startswith(bq) or bq.startswith(input_q):
                        found = bq
                        break
                if found:
                    output_q = found
                    aux_str = ' '.join(bam_qnames[found])
                else:
                    output_q = None
                    aux_str = ""

            mutation = classify_mutation(input_q, output_q)
            preserved = get_preserved_tags(aux_str, comment)

            # Build notes
            notes = ""
            if mutation == "TRUNCATED_WHITESPACE":
                notes = f"Stripped after first whitespace; input had comment: {comment[:40]!r}"
            elif mutation == "MUTATED_CHARS":
                notes = "Spaces converted to underscores"
            elif mutation == "DROPPED":
                notes = "Read absent from BAM output"
            elif mutation == "SUFFIX_ADDED":
                suffix = output_q[len(input_q):]
                notes = f"Suffix appended: {suffix!r}"
            elif mutation == "OTHER":
                notes = f"input_len={len(input_q)} output_len={len(output_q) if output_q else 0}"

            display_in = input_q if len(input_q) <= 60 else f"{input_q[:30]}...[{len(input_q)}]"
            display_out = output_q if (output_q and len(output_q) <= 60) else (f"{output_q[:30]}...[{len(output_q)}]" if output_q else "DROPPED")

            rows.append({
                "aligner": aligner_name,
                "test_case": case["description"],
                "input_qname": display_in,
                "output_qname": display_out,
                "mutation_class": mutation,
                "preserved_tags": preserved,
                "notes": notes,
            })

    # Write results
    with open(results_path, 'w') as f:
        f.write(header + '\n')
        for row in rows:
            line = '\t'.join([
                row["aligner"], row["test_case"], row["input_qname"],
                row["output_qname"], row["mutation_class"],
                row["preserved_tags"], row["notes"],
            ])
            f.write(line + '\n')
    print(f"Wrote {len(rows)} rows to {results_path}")
    return rows


if __name__ == "__main__":
    main()
