"""COMPASS short-read aligner command builders (P3).

These aligners need the cluster (binaries + human indices) to RUN, so the
runnable path is verified on Sherlock. What IS locally verifiable — and the
riskiest part of blind, untestable wrappers — is that the assembled command
lines match the COMPASS process_reads_and_align.sh invocations exactly,
including the gz-aware adaptations (RECTIFY feeds gzipped chunks; COMPASS fed
uncompressed renumbered reads).
"""

from rectify.core.align.multi_aligner import (
    _build_star_cmd,
    _build_hisat2_cmd,
    _build_magicblast_cmd,
    _build_gsnap_cmd,
    _compass_index_paths,
    _genome_version,
)


def _joined(cmd):
    return ' '.join(str(c) for c in cmd)


# ── Index path derivation ────────────────────────────────────────────────────

def test_index_paths_match_compass_layout():
    p = _compass_index_paths('/x/genome_references/GRCh38_gencode_v44.fasta', read_length=150)
    assert p.genome_version == 'GRCh38_gencode_v44'
    assert p.star_dir.name == 'STAR_annotated_150_bp_SJDB_index'
    assert p.hisat2_index.as_posix().endswith('HISAT2_annotated_index/GRCh38_gencode_v44')
    assert p.splice_sites.name == 'GRCh38_gencode_v44_splice_sites.txt'
    assert p.blast_index.as_posix().endswith('BLAST/GRCh38_gencode_v44')
    assert p.gsnap_dir.name == 'GSNAP'


def test_genome_version_strips_fasta_exts():
    assert _genome_version('/a/b/GRCh38.fasta') == 'GRCh38'
    assert _genome_version('/a/b/GRCh38.fa.gz') == 'GRCh38'
    assert _genome_version('/a/b/GRCh38.fasta.gz') == 'GRCh38'


# ── STAR ─────────────────────────────────────────────────────────────────────

def test_star_default_cmd():
    cmd = _build_star_cmd('R1.fastq.gz', 'R2.fastq.gz', '/idx/STAR', '/out/s.',
                          threads=8, read_length=150, noncanonical=False)
    s = _joined(cmd)
    assert s.startswith('STAR --runThreadN 8')
    assert '--genomeDir /idx/STAR' in s
    assert '--sjdbOverhang 149' in s              # read_length - 1
    assert '--readFilesIn R1.fastq.gz R2.fastq.gz' in s
    assert '--outFileNamePrefix /out/s.' in s
    assert '--alignEndsType EndToEnd' in s
    assert '--outSAMattributes NH HI NM MD AS nM jM jI XS' in s
    assert '--readFilesCommand zcat' in s         # gz adaptation
    assert '--scoreGapNoncan' not in s            # default mode


def test_star_noncanonical_adds_scoregapnoncan():
    cmd = _build_star_cmd('R1.fastq.gz', 'R2.fastq.gz', '/idx', '/out/s.',
                          threads=4, noncanonical=True)
    assert '--scoreGapNoncan 0' in _joined(cmd)


def test_star_uncompressed_input_no_zcat():
    cmd = _build_star_cmd('R1.fastq', 'R2.fastq', '/idx', '/out/s.', threads=4)
    assert '--readFilesCommand' not in _joined(cmd)


# ── HISAT2 ───────────────────────────────────────────────────────────────────

def test_hisat2_default_cmd():
    cmd = _build_hisat2_cmd('R1.fq.gz', 'R2.fq.gz', '/idx/hs', '/idx/ss.txt',
                            'out.sam', 8, 20, 200000, 'novel.txt', 'sum.txt',
                            noncanonical=False)
    s = _joined(cmd)
    assert s.startswith('hisat2 --known-splicesite-infile /idx/ss.txt')
    assert '--no-softclip' in s
    assert '-x /idx/hs' in s
    assert '-1 R1.fq.gz -2 R2.fq.gz' in s
    assert '-S out.sam' in s
    assert '--min-intronlen 20' in s
    assert '--max-intronlen 200000' in s
    assert '--rna-strandness RF' in s
    assert '--new-summary' in s
    assert '--pen-noncansplice' not in s


def test_hisat2_noncanonical_adds_penalty():
    cmd = _build_hisat2_cmd('R1.fq.gz', 'R2.fq.gz', '/idx/hs', '/idx/ss.txt',
                            'out.sam', 8, 20, 200000, 'novel.txt', 'sum.txt',
                            noncanonical=True)
    assert '--pen-noncansplice 0' in _joined(cmd)


# ── Magic-BLAST ──────────────────────────────────────────────────────────────

def test_magicblast_cmd():
    cmd = _build_magicblast_cmd('R1.fq.gz', 'R2.fq.gz', '/idx/blastdb', 'out.sam', threads=12)
    s = _joined(cmd)
    assert s.startswith('magicblast -query R1.fq.gz -query_mate R2.fq.gz')
    assert '-db /idx/blastdb' in s
    assert '-md_tag -fr -no_query_id_trim' in s
    assert '-infmt fastq' in s
    assert '-num_threads 12' in s
    assert '-max_db_word_count 10' in s
    assert '-out out.sam' in s


# ── GSNAP ────────────────────────────────────────────────────────────────────

def test_gsnap_cmd_gzipped():
    cmd = _build_gsnap_cmd('R1.fq.gz', 'R2.fq.gz', '/idx/GSNAP', 'GRCh38', 'out.sam', 8)
    s = _joined(cmd)
    assert s.startswith('gsnap -D /idx/GSNAP -d GRCh38 --use-splicing=GRCh38')
    assert 'R1.fq.gz R2.fq.gz' in s
    assert '--output-file out.sam' in s
    assert '--nthreads=8' in s
    assert '--ambig-splice-noclip' in s
    assert '--novelsplicing=1' in s
    assert '--add-paired-nomappers' in s
    assert '--sam-extended-cigar' in s
    assert '--format=sam' in s
    assert '--gunzip' in s                        # gz adaptation


def test_gsnap_uncompressed_no_gunzip():
    cmd = _build_gsnap_cmd('R1.fq', 'R2.fq', '/idx/GSNAP', 'GRCh38', 'out.sam', 8)
    assert '--gunzip' not in _joined(cmd)
