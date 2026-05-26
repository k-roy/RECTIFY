# Protocol flag: `--dT-primed-cDNA`

**Nanopore direct RNA-seq (DRS):** sequences the entire molecule from the
adapter through the poly-A tail into the RNA body. The poly-A tail IS in
the read. This is the default mode — no protocol flag needed.

**dT-primed cDNA (QuantSeq, oligo-dT cDNA):** the oligo-dT primer
hybridizes at the start of the poly-A tract; sequencing starts at the
first non-A base. The poly-A tail is NOT in the read. Pass
`--dT-primed-cDNA` to enable Module AG (AG mispriming detection), a
cDNA-synthesis artifact caused by reverse-transcriptase slippage at
internal AG runs. This artifact does not occur in DRS and is disabled by
default.

`--polya-sequenced` is a deprecated alias for `--dT-primed-cDNA` and will
be removed in a future version.

## Module activation by protocol

| Module | DRS (default) | dT-primed cDNA (`--dT-primed-cDNA`) | Disabled by `--short-read`? |
|---|:---:|:---:|:---:|
| A-tract walk-back (2E) | ✓ | ✓ | **No** |
| Poly(A) trimming (2B) | ✓ | ✓ | Yes |
| Indel correction (2C) | ✓ (if MD tags present) | ✓ (if MD tags present) | Yes |
| Soft-clip rescue (2G) | ✓ | ✓ | No |
| AG mispriming | ✗ | ✓ (long- and short-read) | No |
| Variant-aware rescue | ✓ | ✓ | No |

**`--short-read` disables only** poly(A) trimming and indel correction,
which require a sequenced poly(A) tail.

**A-tract walk-back (2E) is enabled for short reads** — CSP internal
priming in QuantSeq REV causes poly-A complement T's that, after reverse
complementation for alignment, appear as A's and slide into downstream
genomic A-runs, shifting the reported 3' end rightward.
`find_atract_boundaries()` uses the reference genome, not the read
sequence, so it corrects this regardless of read length. Specifically:
reads ending in `AAAAAAAAAT` (A-run terminated by a genomic T) are
walked back to the last non-A position upstream of the A-run.

**AG mispriming — why DRS is exempt:** DRS (Oxford Nanopore direct RNA)
has no priming step — the RNA is sequenced directly. Oligo-dT
mispriming cannot occur. The module is therefore off by default and
enabled only with `--dT-primed-cDNA`.

**AG mispriming — why short reads need it most:** For QuantSeq/Illumina
dT-primed data, the oligo-dT primer can misprime onto genomic A/G-rich
regions, producing false 3' end sites. The reads stop at the cleavage
site, so the A/G-rich genomic context is **downstream of the reported 3'
end** — absent from the read itself. The module detects mispriming by
calling `get_downstream_sequence()` on the genome, not by inspecting
read sequence. This genomic-context approach works identically for short
and long reads, and short-read QuantSeq data is the primary use case.

Module 2G (`rescue_softclip_at_homopolymer`) runs for all protocols. It
detects 3' soft-clips ≥ 3 bp adjacent to a genomic homopolymer and
extends the 3' end outward. It takes priority over poly-A walk-back
(Module 2E) to prevent opposite-direction corrections from cancelling.
