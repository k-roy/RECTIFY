"""`rectify run-all --netseq` — the NET-seq mode (planning 835 WP2 / WP3, unit 856).

Three layers:

* **Command builders.** The planning-829 §4 absolute STAR floors are the expensive part of this
  mode — they cost 42 %% of one library's "unique" alignments to learn — and STAR is not on CI, so
  they are asserted flag by flag on the pure builder. Same for the ABSENCE of ``-u`` in the
  cutadapt command: trimming the randomer would delete six genuine nucleotides from the 58-60 %%
  randomer-free class of a mixed library.
* **Namespace.** ``rectify netseq``'s Namespace is built from its OWN parser, so a new netseq flag
  cannot be silently dropped (the 832 G-1 / 834 §6.4 / 836 CP0b defect shape). The mode's own
  overrides — read5p, corrected tracks, Pol III included, counts not RPM — are asserted.
* **End to end from a BAM**, on a synthetic genome with one intron: the mode runs, writes the
  raw/corrected/deconv tracks and both JSONs, and rescues a read whose 3' end sits 5 nt into
  exon 2. The FASTQ→BAM stages need STAR + cutadapt and are skipped when they are absent.
"""
from __future__ import annotations

import argparse
import json
import random
import shutil
from pathlib import Path

import pysam
import pytest

from rectify.core.commands.run.netseq_pipeline import (
    CHURCHMAN_3P_LINKER,
    NETSEQ_STAR_FLOORS,
    build_combined_genome,
    build_cutadapt_cmd,
    build_netseq_namespace,
    build_netseq_star_cmd,
    build_star_index_cmd,
    parse_idxstats_spikein,
    run_netseq_pipeline,
)

# ══════════════════════════════════════════════════════════════════════════════
# 1 · Command builders
# ══════════════════════════════════════════════════════════════════════════════

#: planning 829 §4, verbatim. A change to any of these must be a deliberate edit HERE too.
_REQUIRED_STAR_FLOORS = {
    '--outFilterMatchNmin': '20',
    '--outFilterMismatchNoverLmax': '0.06',
    '--outFilterMismatchNmax': '3',
    '--alignSJoverhangMin': '12',
    '--outFilterIntronMotifs': 'RemoveNoncanonical',
    '--outFilterMultimapNmax': '1',
    '--alignEndsType': 'Local',
    '--alignIntronMax': '2000',
}


@pytest.mark.parametrize("flag,value", sorted(_REQUIRED_STAR_FLOORS.items()))
def test_star_command_carries_every_829_floor(flag, value):
    cmd = build_netseq_star_cmd('r.fq.gz', Path('/idx'), Path('/out/s.'), threads=4, gzipped=True)
    assert flag in cmd, f"{flag} missing — the 829 §4 absolute floors are the point of this mode"
    assert cmd[cmd.index(flag) + 1] == value


def test_star_floor_table_and_command_agree():
    cmd = build_netseq_star_cmd('r.fq', Path('/idx'), Path('/out/s.'))
    for flag, values in NETSEQ_STAR_FLOORS.items():
        i = cmd.index(flag)
        assert tuple(cmd[i + 1:i + 1 + len(values)]) == values


def test_star_command_writes_md_tags():
    """MD is load-bearing: the terminal-oligo(A) X-mismatch correction parses it for M ops."""
    cmd = build_netseq_star_cmd('r.fq', Path('/idx'), Path('/out/s.'))
    attrs = cmd[cmd.index('--outSAMattributes') + 1:]
    assert 'MD' in attrs[:6] and 'NM' in attrs[:6]


def test_star_command_zcat_only_for_gzipped_input():
    assert '--readFilesCommand' in build_netseq_star_cmd('r.fq.gz', Path('/i'), Path('/o/'),
                                                         gzipped=True)
    assert '--readFilesCommand' not in build_netseq_star_cmd('r.fq', Path('/i'), Path('/o/'),
                                                             gzipped=False)


def test_cutadapt_never_trims_the_randomer():
    """🔴 The one thing this trimmer must NOT do (planning 829 §4: 58-60 %% of reads have none)."""
    cmd = build_cutadapt_cmd(Path('r.fq.gz'), Path('t.fq.gz'))
    assert '-u' not in cmd and '--cut' not in cmd
    assert not any(str(c).startswith('--front') or str(c) == '-g' for c in cmd)


def test_cutadapt_defaults_match_the_validated_829_pipeline():
    cmd = build_cutadapt_cmd(Path('r.fq.gz'), Path('t.fq.gz'), threads=8)
    assert cmd[cmd.index('-a') + 1] == CHURCHMAN_3P_LINKER == 'ATCTCGTATGCCGTCTTCTGCTTG'
    assert '--nextseq-trim=20' in cmd
    assert cmd[cmd.index('-m') + 1] == '18'
    assert cmd[cmd.index('-j') + 1] == '8'


def test_cutadapt_quality_cutoff_replaces_nextseq_trim():
    cmd = build_cutadapt_cmd(Path('r.fq'), None, quality_cutoff=20)
    assert cmd[cmd.index('-q') + 1] == '20'
    assert not any(str(c).startswith('--nextseq-trim') for c in cmd)
    assert cmd[cmd.index('-o') + 1] == '-'   # stdout, for streaming into STAR


def test_star_index_cmd_uses_the_small_genome_sa_size():
    cmd = build_star_index_cmd(Path('g.fa'), Path('/idx'), threads=4)
    assert cmd[cmd.index('--genomeSAindexNbases') + 1] == '11'
    assert '--runMode' in cmd and cmd[cmd.index('--runMode') + 1] == 'genomeGenerate'
    assert '--sjdbGTFfile' not in cmd   # annotation-free index (829); junctions found de novo


# ══════════════════════════════════════════════════════════════════════════════
# 2 · Spike-in accounting
# ══════════════════════════════════════════════════════════════════════════════

def test_build_combined_genome_prefixes_only_the_spikein(tmp_path):
    (tmp_path / 'g.fa').write_text(">chrI a yeast chromosome\nACGT\n>chrII\nTTTT\n")
    (tmp_path / 'sp.fa').write_text(">I Schizosaccharomyces_pombe\nGGGG\n>II\nCCCC\n>III\nAAAA\n")
    counts = build_combined_genome(tmp_path / 'g.fa', tmp_path / 'sp.fa', tmp_path / 'c.fa')
    assert counts == {'primary_contigs': 2, 'spikein_contigs': 3}
    names = [ln[1:].split()[0] for ln in (tmp_path / 'c.fa').read_text().splitlines()
             if ln.startswith('>')]
    assert names == ['chrI', 'chrII', 'Sp_I', 'Sp_II', 'Sp_III']
    # the primary genome is byte-identical, so no coordinate can move
    assert '>chrI a yeast chromosome' in (tmp_path / 'c.fa').read_text()


def test_build_combined_genome_refuses_a_prefix_collision(tmp_path):
    (tmp_path / 'g.fa').write_text(">Sp_weird\nACGT\n")
    (tmp_path / 'sp.fa').write_text(">I\nGGGG\n")
    with pytest.raises(ValueError, match='already starts with the spike-in prefix'):
        build_combined_genome(tmp_path / 'g.fa', tmp_path / 'sp.fa', tmp_path / 'c.fa')


def test_parse_idxstats_spikein():
    text = ("chrI\t230218\t1000\t0\n"
            "chrII\t813184\t2000\t0\n"
            "Sp_I\t5579133\t300\t0\n"
            "Sp_II\t4539804\t0\t0\n"
            "*\t0\t0\t500\n")
    out = parse_idxstats_spikein(text)
    assert out['primary_mapped'] == 3000
    assert out['spikein_mapped'] == 300
    assert out['spikein_fraction'] == pytest.approx(300 / 3300)
    assert out['primary_contigs'] == ['chrI', 'chrII']       # the view filter, order preserved
    assert out['spikein_contigs'] == ['Sp_I', 'Sp_II']


def test_parse_idxstats_no_spikein_is_a_zero_fraction_not_a_crash():
    out = parse_idxstats_spikein("chrI\t100\t7\t0\n*\t0\t0\t1\n")
    assert out['spikein_mapped'] == 0 and out['spikein_fraction'] == 0.0
    assert out['spikein_contigs'] == []


# ══════════════════════════════════════════════════════════════════════════════
# 3 · The netseq Namespace
# ══════════════════════════════════════════════════════════════════════════════

def _run_args(**kw):
    base = dict(organism=None, netseq=True)
    base.update(kw)
    return argparse.Namespace(**base)


def test_netseq_namespace_carries_the_mode_defaults(tmp_path):
    ns = build_netseq_namespace(_run_args(), tmp_path / 's.bam', tmp_path / 'g.fa',
                                tmp_path / 'a.gff', tmp_path / 'out')
    assert ns.rna3p_at == 'read5p'          # Churchman: gene strand = INVERSE of BAM strand
    assert ns.track_position == 'corrected'  # planning 836 rescue/tail output reaches the tracks
    assert ns.include_pol3 is True           # else every snRNA is dropped, incl. the U5 QC peak
    assert ns.no_rpm_normalize is True       # counts: NNLS does not conserve mass
    assert ns.include_rdna is False          # rDNA still excluded
    assert ns.umi_length == 0 and ns.dedup is False
    assert ns.walkback_requires_clip_a is True
    assert ns.gff == tmp_path / 'a.gff'      # the rescue pool has a source


def test_netseq_namespace_inherits_every_subparser_default(tmp_path):
    """The anti-hand-building guarantee: defaults the mode never mentions still arrive."""
    ns = build_netseq_namespace(_run_args(), tmp_path / 's.bam', None, None, tmp_path / 'o')
    for dest in ('rescue_min_k', 'rescue_min_k_with_remainder', 'rescue_max_intronic',
                 'pool_include_trna', 'min_atract_length', 'no_deconvolution', 'min_mapq',
                 'output_format', 'umi_source', 'exclude_mito', 'pol3_flanking'):
        assert hasattr(ns, dest), f"{dest} dropped — Namespace must come from add_netseq_parser"
    assert ns.rescue_min_k == 1 and ns.rescue_min_k_with_remainder == 4


@pytest.mark.parametrize("run_dest,netseq_dest,value", [
    ('netseq_umi_length', 'umi_length', 6),
    ('netseq_rescue_min_k', 'rescue_min_k', 5),
    ('netseq_rescue_min_k_with_remainder', 'rescue_min_k_with_remainder', 7),
    ('netseq_rescue_max_intronic', 'rescue_max_intronic', 4),
    ('netseq_track_position', 'track_position', 'raw'),
    ('netseq_min_mapq', 'min_mapq', 10),
    ('netseq_max_reads', 'max_reads', 1234),
    ('netseq_output_format', 'output_format', ['bedgraph', 'bigwig']),
    ('netseq_pool_include_trna', 'pool_include_trna', True),
    ('netseq_no_deconvolution', 'no_deconvolution', True),
    ('netseq_exclude_mito', 'exclude_mito', True),
])
def test_netseq_flags_are_forwarded(tmp_path, run_dest, netseq_dest, value):
    ns = build_netseq_namespace(_run_args(**{run_dest: value}), tmp_path / 's.bam', None, None,
                                tmp_path / 'o')
    assert getattr(ns, netseq_dest) == value


def test_netseq_inverted_flags(tmp_path):
    ns = build_netseq_namespace(_run_args(netseq_exclude_pol3=True, netseq_rpm=True,
                                          netseq_walkback_unconditional=True),
                                tmp_path / 's.bam', None, None, tmp_path / 'o')
    assert ns.include_pol3 is False
    assert ns.no_rpm_normalize is False
    assert ns.walkback_requires_clip_a is False


def test_dedup_without_a_umi_length_is_refused(tmp_path):
    with pytest.raises(ValueError, match='requires --netseq-umi-length'):
        build_netseq_namespace(_run_args(netseq_dedup=True), tmp_path / 's.bam', None, None,
                               tmp_path / 'o')


# ══════════════════════════════════════════════════════════════════════════════
# 4 · End to end from a BAM, on a synthetic one-intron genome
# ══════════════════════════════════════════════════════════════════════════════

CHROM = 'chrT'
EXON1, INTRON, EXON2 = 400, 120, 400
INTRON_START, INTRON_END = EXON1, EXON1 + INTRON


def _genome_seq(seed=11):
    rng = random.Random(seed)
    exon1 = ''.join(rng.choice('ACGT') for _ in range(EXON1))
    intron = 'GT' + ''.join(rng.choice('ACGT') for _ in range(INTRON - 4)) + 'AG'
    exon2 = ''.join(rng.choice('ACGT') for _ in range(EXON2))
    return exon1 + intron + exon2


@pytest.fixture(scope='module')
def fixture_dir(tmp_path_factory):
    """A minimal NET-seq input set: genome FASTA + GFF + a BAM of synthetic + gene reads."""
    d = tmp_path_factory.mktemp('netseq_fixture')
    seq = _genome_seq()
    fasta = d / 'mini.fa'
    with open(fasta, 'w') as fh:
        fh.write(f">{CHROM}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    pysam.faidx(str(fasta))

    # GFF: one gene / mRNA / two exons / one intron, 1-based inclusive.
    gff = d / 'mini.gff'
    gff.write_text(
        "##gff-version 3\n"
        f"{CHROM}\ttest\tgene\t1\t{len(seq)}\t.\t+\t.\tID=GENE1\n"
        f"{CHROM}\ttest\tmRNA\t1\t{len(seq)}\t.\t+\t.\tID=GENE1_mRNA;Parent=GENE1\n"
        f"{CHROM}\ttest\texon\t1\t{INTRON_START}\t.\t+\t.\tID=e1;Parent=GENE1_mRNA\n"
        f"{CHROM}\ttest\tintron\t{INTRON_START + 1}\t{INTRON_END}\t.\t+\t."
        f"\tID=i1;Parent=GENE1_mRNA\n"
        f"{CHROM}\ttest\texon\t{INTRON_END + 1}\t{len(seq)}\t.\t+\t.\tID=e2;Parent=GENE1_mRNA\n"
    )

    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": CHROM, "LN": len(seq)}]}
    bam = d / 'mini.bam'
    with pysam.AlignmentFile(str(bam), 'wb', header=header) as out:
        # (a) 40 plain gene-body reads: gene '+' => the read maps REVERSE (Churchman).
        for i in range(40):
            start = 100 + i * 3
            r = pysam.AlignedSegment(out.header)
            r.query_name = f"body{i}"
            r.reference_id, r.reference_start, r.flag, r.mapping_quality = 0, start, 16, 255
            r.query_sequence = seq[start:start + 30]
            r.cigartuples = [(0, 30)]
            r.query_qualities = [30] * 30
            out.write(r)
        # (b) 12 reads whose RNA 3' end is 5 nt into exon 2 but whose overhang STAR soft-clipped:
        #     aligned block ends at exon1_last, the 5 exon-2 bases sit in the RIGHT clip.
        #     Correct call = exon2_first + 4; uncorrected = a FALSE splicing intermediate at the 5'SS.
        for i in range(12):
            r = pysam.AlignedSegment(out.header)
            r.query_name = f"overhang{i}"
            r.reference_id, r.reference_start, r.flag, r.mapping_quality = 0, 340 + i, 16, 255
            aligned = seq[340 + i:INTRON_START]
            clip = seq[INTRON_END:INTRON_END + 5]
            r.query_sequence = aligned + clip
            r.cigartuples = [(0, len(aligned)), (4, 5)]
            r.query_qualities = [30] * len(r.query_sequence)
            out.write(r)
        # (c) 8 genuine splicing intermediates: 3' end exactly at exon1_last, no clip. Must NOT move.
        for i in range(8):
            r = pysam.AlignedSegment(out.header)
            r.query_name = f"intermediate{i}"
            r.reference_id, r.reference_start, r.flag, r.mapping_quality = 0, 355 + i, 16, 255
            r.query_sequence = seq[355 + i:INTRON_START]
            r.cigartuples = [(0, INTRON_START - (355 + i))]
            r.query_qualities = [30] * len(r.query_sequence)
            out.write(r)
    pysam.index(str(bam))
    return {'dir': d, 'genome': fasta, 'gff': gff, 'bam': bam, 'seq': seq}


def _mode_args(**kw):
    base = dict(organism=None, netseq=True, threads=1,
                netseq_output_format=['bedgraph'], netseq_no_deconvolution=False)
    base.update(kw)
    return argparse.Namespace(**base)


@pytest.fixture(scope='module')
def mode_output(fixture_dir, tmp_path_factory):
    out = tmp_path_factory.mktemp('netseq_out')
    rc = run_netseq_pipeline(
        _mode_args(),
        input_path=fixture_dir['bam'], input_type='bam',
        output_dir=out, work_dir=out, sample_id='mini',
        genome_path=fixture_dir['genome'], annotation_path=fixture_dir['gff'], threads=1,
    )
    assert rc == 0
    return out


def test_mode_writes_the_expected_track_set(mode_output):
    for name in ('mini.raw.plus.bedgraph', 'mini.raw.minus.bedgraph',
                 'mini.corrected.plus.bedgraph', 'mini.corrected.minus.bedgraph',
                 'mini.netseq_summary.json', 'mini.runall_netseq.json'):
        p = mode_output / name
        assert p.exists(), f"{name} missing from {sorted(x.name for x in mode_output.iterdir())}"


def test_mode_records_its_own_flags(mode_output):
    payload = json.loads((mode_output / 'mini.runall_netseq.json').read_text())
    flags = payload['netseq_flags']
    assert flags['rna3p_at'] == 'read5p'
    assert flags['track_position'] == 'corrected'
    assert flags['include_pol3'] is True
    assert flags['no_rpm_normalize'] is True
    assert payload['mode'] == 'run-all --netseq'


def test_mode_rescues_the_exon2_overhang_and_keeps_the_intermediates(mode_output):
    """Kevin's (a), on a fixture: the 12 clipped reads move to exon 2, the 8 real ones stay."""
    summary = json.loads((mode_output / 'mini.netseq_summary.json').read_text())
    assert summary['reads_near_donor'] >= 20
    assert summary['rescued'] == 12, summary
    assert summary['exon1_end'] == 8, summary
    assert summary['decoy_rescued'] == 0, "a decoy acceptor must not produce these rescues"


def test_mode_moves_the_5prime_splice_site_signal_to_exon_2(mode_output, fixture_dir):
    """The artifact/signal collision, measured on the tracks themselves."""
    def load(path):
        counts = {}
        for line in path.read_text().splitlines():
            if line.startswith('track') or not line.strip():
                continue
            chrom, start, _end, val = line.split('\t')[:4]
            counts[int(start)] = float(val)
        return counts

    raw = load(mode_output / 'mini.raw.plus.bedgraph')
    corrected = load(mode_output / 'mini.corrected.plus.bedgraph')
    exon1_last, exon2_landing = INTRON_START - 1, INTRON_END + 4
    # BEFORE: all 20 donor-adjacent reads pile onto the 5' splice site.
    assert raw.get(exon1_last, 0) == 20
    assert raw.get(exon2_landing, 0) == 0
    # AFTER: only the 8 genuine intermediates remain there; the 12 land in exon 2.
    assert corrected.get(exon1_last, 0) == 8
    assert corrected.get(exon2_landing, 0) == 12


def test_mode_conserves_mass_between_the_raw_and_corrected_tracks(mode_output):
    """planning 829 §4's own gate: correction MOVES reads, it never creates or drops them."""
    def total(*names):
        s = 0.0
        for name in names:
            for line in (mode_output / name).read_text().splitlines():
                if line.startswith('track') or not line.strip():
                    continue
                s += float(line.split('\t')[3])
        return s
    raw = total('mini.raw.plus.bedgraph', 'mini.raw.minus.bedgraph')
    corr = total('mini.corrected.plus.bedgraph', 'mini.corrected.minus.bedgraph')
    assert raw == corr == 60


# ══════════════════════════════════════════════════════════════════════════════
# 5 · The `--netseq-dir` handoff (planning 834 gap 5 / 835 WP2 (e))
# ══════════════════════════════════════════════════════════════════════════════

def test_netseq_dir_loader_selects_one_strand_and_one_track_kind(tmp_path):
    """Scope item (e), with evidence: the loader must not sum the six BigWigs this mode emits."""
    pyBigWig = pytest.importorskip('pyBigWig')
    from rectify.core.netseq.netseq_refiner import NetseqLoader

    names = ['wt.raw.plus.bw', 'wt.raw.minus.bw', 'wt.corrected.plus.bw',
             'wt.corrected.minus.bw', 'wt.deconv.plus.bw', 'wt.deconv.minus.bw']
    for i, name in enumerate(names, start=1):
        bw = pyBigWig.open(str(tmp_path / name), 'w')
        bw.addHeader([('chrI', 1000)])
        bw.addEntries(['chrI'], [100], ends=[110], values=[float(i)])
        bw.close()

    loader = NetseqLoader()
    loader.load_directory(str(tmp_path), pattern='*.bw')
    assert len(loader.bigwigs) == 6
    # raw.plus alone = 1.0. Summing all six would give 21.0; summing both strands of raw, 3.0.
    sig = loader.get_signal('chrI', 100, 110, '+')
    assert float(sig.max()) == pytest.approx(1.0), (
        "get_signal summed more than the raw plus-strand track (planning 834 §4 / 836 D6)"
    )
    assert float(loader.get_signal('chrI', 100, 110, '-').max()) == pytest.approx(2.0)


# ══════════════════════════════════════════════════════════════════════════════
# 6 · FASTQ → BAM stages (need the binaries; skipped on CI)
# ══════════════════════════════════════════════════════════════════════════════

@pytest.mark.skipif(shutil.which('STAR') is None or shutil.which('cutadapt') is None,
                    reason='needs STAR and cutadapt on PATH')
def test_fastq_entrypoint_builds_an_index_and_runs(fixture_dir, tmp_path):
    seq = fixture_dir['seq']
    fq = tmp_path / 'mini.fastq'
    with open(fq, 'w') as fh:
        for i in range(200):
            start = 100 + (i % 200)
            read = seq[start:start + 35] + CHURCHMAN_3P_LINKER
            fh.write(f"@r{i}\n{read}\n+\n{'I' * len(read)}\n")
    out = tmp_path / 'out'
    rc = run_netseq_pipeline(
        _mode_args(star_index_dir=tmp_path / 'idx'),
        input_path=fq, input_type='fastq', output_dir=out, work_dir=tmp_path / 'work',
        sample_id='mini', genome_path=fixture_dir['genome'],
        annotation_path=fixture_dir['gff'], threads=1,
    )
    assert rc == 0
    assert (tmp_path / 'idx' / 'SAindex').exists()
    assert (tmp_path / 'work' / 'netseq_align' / 'mini.spikein.tsv').exists()
    assert (out / 'mini.corrected.plus.bedgraph').exists()


def test_netseq_skip_trim_is_wired(fixture_dir, tmp_path, monkeypatch):
    """A flag that parses but is never read is the same silent no-op class as 832 G-1."""
    import rectify.core.commands.run.netseq_pipeline as np_mod

    called = {}
    monkeypatch.setattr(np_mod, 'netseq_trim',
                        lambda *a, **k: called.setdefault('trim', True))
    monkeypatch.setattr(np_mod, 'ensure_star_index', lambda *a, **k: (_ for _ in ()).throw(
        RuntimeError('reached align stage')))
    fq = tmp_path / 'r.fastq'
    fq.write_text("@a\nACGT\n+\nIIII\n")
    with pytest.raises(RuntimeError, match='reached align stage'):
        np_mod.run_netseq_pipeline(
            _mode_args(netseq_skip_trim=True), input_path=fq, input_type='fastq',
            output_dir=tmp_path / 'o', work_dir=tmp_path / 'w', sample_id='r',
            genome_path=fixture_dir['genome'], annotation_path=fixture_dir['gff'], threads=1)
    assert 'trim' not in called, "--netseq-skip-trim did not skip the trim stage"
