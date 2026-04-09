# Workshop: Running RECTIFY Tests with GitHub Copilot Agent

RECTIFY is pre-installed on the Hoffman2 cluster. Your VS Code should already
be connected to Hoffman2 via the Remote-SSH extension. The Copilot agent runs
commands directly on the cluster in your integrated terminal.

---

## The Test Command

The full validation command is:

```bash
rectify test --full --keep --output-dir ~/rectify_test --threads 4
```

- `--full` — runs the complete pipeline from FASTQ: align → correct (requires minimap2, which is available on the server)
- `--keep` — retains the output directory so you can inspect the results
- `--output-dir ~/rectify_test` — saves outputs to your home directory
- `--threads 4` — use 4 CPU threads

A passing run ends with:

```
==============================================================
  Result  : 42/42 passed
  Time    : ~60–90s
  Output  : ~/rectify_test
==============================================================
```

---

## Using the Copilot Agent

1. Open the Copilot Chat panel (`Ctrl+Alt+I` / `Cmd+Alt+I`)
2. Switch the mode dropdown at the top of the input from **Ask** to **Agent**
3. Paste this prompt:

> Run `rectify test --full --keep --output-dir ~/rectify_test --threads 4`
> in the terminal. Report whether all 42 tests passed. If any failed, show me
> the read name, expected corrected position, and actual corrected position.

---

## Exploring the Output

Once the test completes, the output directory contains:

```
~/rectify_test/
├── corrected_3ends.tsv       # per-read correction results
├── corrected_3ends_stats.tsv # processing statistics
├── rectified.bam             # hard-clipped corrected BAM
├── rectified.bam.bai
├── rectified_softclip.bam    # soft-clipped BAM (poly-A bases visible in IGV)
├── rectified_softclip.bam.bai
├── minimap2.bam              # raw aligner output
├── minimap2.bam.bai
├── mapPacBio.bam
└── gapmm2.bam
```

Try these follow-up prompts:

> Show me the `read_id`, `strand`, `raw_3prime`, `corrected_3prime`, and
> `correction_category` columns from `~/rectify_test/corrected_3ends.tsv`,
> grouped by correction category.

> Using samtools, show the read count in each BAM file in `~/rectify_test/`.

---

## Correction Categories in the Output

| Category | Description |
|----------|-------------|
| `cat1` | Poly-A walkback — 3' end walked back through aligned A-runs |
| `cat2` | Soft-clip rescue — soft-clipped bases extend the 3' end |
| `cat3` | Splice rescue — 5' end corrected across an unaligned intron |
| `cat4` | False junction removal — spurious N-op absorbed, poly-A walkback applied |

Each category has both plus- and minus-strand examples in the validation set.
