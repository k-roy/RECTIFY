"""Paired-end (short-read) chunking: shared RN, same-chunk mates, round-trip.

Covers split_fastq_paired — the pair-aware variant of split_fastq that keeps
both mates of a pair together (shared read_num, same chunk) so proper-pairs
survive round-robin chunking.
"""

import gzip

import pytest

from rectify.core.chunking.sidecar import (
    ReadNumSidecar,
    reconstruct_posthoc_sidecar_from_chunks,
)
from rectify.core.commands.split_command import split_fastq_paired


def _write_pair(r1, r2, n=100, mate_style='slash'):
    """Write n synthetic read pairs.

    mate_style='slash' → @readN/1, @readN/2 (suffix stripped on write).
    mate_style='casava' → @readN 1:N:0:AC, @readN 2:N:0:AC (token already shared).
    """
    with open(r1, 'w') as f1, open(r2, 'w') as f2:
        for i in range(n):
            seq1 = 'ACGT' + 'A' * (i % 5)
            seq2 = 'TGCA' + 'C' * (i % 3)
            if mate_style == 'slash':
                h1, h2 = f'@read{i}/1', f'@read{i}/2'
            else:
                h1, h2 = f'@read{i} 1:N:0:AC', f'@read{i} 2:N:0:AC'
            f1.write(f'{h1}\n{seq1}\n+\n{"I" * len(seq1)}\n')
            f2.write(f'{h2}\n{seq2}\n+\n{"I" * len(seq2)}\n')


def _read_chunk(path):
    """Return list of (header, seq, qual) from a gzipped chunk FASTQ."""
    out = []
    with gzip.open(path, 'rt') as fh:
        while True:
            h = fh.readline().rstrip('\n')
            if not h:
                break
            seq = fh.readline().rstrip('\n')
            fh.readline()
            qual = fh.readline().rstrip('\n')
            out.append((h, seq, qual))
    return out


def _rn_of(header):
    return int(header.split('RN:i:')[1].split('\t')[0])


def _bare_qname(header):
    return header.split('\t')[0][1:]  # strip leading '@', take token before first tab


@pytest.mark.parametrize('mate_style', ['slash', 'casava'])
def test_paired_split_shared_rn_same_chunk(tmp_path, mate_style):
    r1 = tmp_path / 'samp_R1.fastq'
    r2 = tmp_path / 'samp_R2.fastq'
    _write_pair(r1, r2, n=100, mate_style=mate_style)

    n_chunks = 4
    r1_paths, r2_paths = split_fastq_paired(r1, r2, tmp_path / 'chunks', n_chunks)
    assert len(r1_paths) == n_chunks and len(r2_paths) == n_chunks

    sc = ReadNumSidecar.open(tmp_path / 'chunks' / 'samp.read_num_sidecar.parquet')
    assert len(sc) == 100  # one row PER PAIR, not per record

    seen_rn = set()
    for k in range(n_chunks):
        rec1 = _read_chunk(r1_paths[k])
        rec2 = _read_chunk(r2_paths[k])
        # Both mate chunks have the same number of records (mates never separated)
        assert len(rec1) == len(rec2)
        for local_idx, ((h1, s1, q1), (h2, s2, q2)) in enumerate(zip(rec1, rec2)):
            rn1, rn2 = _rn_of(h1), _rn_of(h2)
            # Shared RN across mates
            assert rn1 == rn2
            # Round-robin read_num: read_num = local_idx * n_chunks + chunk_idx
            assert rn1 == local_idx * n_chunks + k
            # Both mates carry the SAME bare QNAME (aligners pair by identical QNAME)
            assert _bare_qname(h1) == _bare_qname(h2) == f'read{rn1}'
            seen_rn.add(rn1)
            # Sidecar fingerprint covers R1
            assert sc.verify(rn1, s1, q1).ok
    assert seen_rn == set(range(100))


def test_paired_round_trip_reconstruct(tmp_path):
    r1 = tmp_path / 'samp_R1.fastq'
    r2 = tmp_path / 'samp_R2.fastq'
    _write_pair(r1, r2, n=100, mate_style='slash')

    n_chunks = 4
    r1_paths, _ = split_fastq_paired(r1, r2, tmp_path / 'chunks', n_chunks)

    forward = ReadNumSidecar.open(tmp_path / 'chunks' / 'samp.read_num_sidecar.parquet')

    # Reconstruct the sidecar purely from the R1 chunk FASTQs (the legacy
    # round-robin formula read_num = local_idx * n_chunks + chunk_idx still holds
    # per pair). qname → read_num mapping must match the forward sidecar.
    recon_path = tmp_path / 'recon.parquet'
    reconstruct_posthoc_sidecar_from_chunks(
        [str(p) for p in r1_paths], recon_path,
        sample_id='samp', n_chunks=n_chunks,
    )
    recon = ReadNumSidecar.open(recon_path)
    assert len(recon) == len(forward) == 100
    for rn in range(100):
        assert recon.lookup(rn).original_qname == forward.lookup(rn).original_qname == f'read{rn}'


def test_paired_mate_desync_raises(tmp_path):
    r1 = tmp_path / 'samp_R1.fastq'
    r2 = tmp_path / 'samp_R2.fastq'
    # Deliberately misalign: R2 read 5 is a different fragment id
    with open(r1, 'w') as f1, open(r2, 'w') as f2:
        for i in range(10):
            f1.write(f'@frag{i}/1\nACGT\n+\nIIII\n')
            j = i if i != 5 else 999
            f2.write(f'@frag{j}/2\nTGCA\n+\nIIII\n')
    with pytest.raises(ValueError, match='Mate qname mismatch'):
        split_fastq_paired(r1, r2, tmp_path / 'chunks', 2)


def test_paired_count_mismatch_raises(tmp_path):
    r1 = tmp_path / 'samp_R1.fastq'
    r2 = tmp_path / 'samp_R2.fastq'
    _write_pair(r1, r2, n=10)
    # Truncate R2 to 8 pairs
    lines = (tmp_path / 'samp_R2.fastq').read_text().splitlines(keepends=True)
    (tmp_path / 'samp_R2.fastq').write_text(''.join(lines[: 8 * 4]))
    with pytest.raises(ValueError, match='record count mismatch'):
        split_fastq_paired(r1, r2, tmp_path / 'chunks', 2)
