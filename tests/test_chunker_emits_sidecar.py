import gzip

from rectify.core.chunking.sidecar import ReadNumSidecar
from rectify.core.commands.split_command import split_fastq


def _write_fastq(path, n=100):
    with open(path, 'w') as fh:
        for i in range(n):
            comment = '\tXA:Z:foo\tXC:i:42' if i % 2 else ''
            fh.write(f'@read{i}{comment}\nACGT\n+\nIIII\n')


def test_chunker_emits_rn_headers_and_sidecar(tmp_path):
    src = tmp_path / 'sample.fastq'
    _write_fastq(src, 100)
    chunks = split_fastq(src, tmp_path / 'chunks', 4, prefix='sample')
    sidecar = tmp_path / 'chunks' / 'sample.read_num_sidecar.parquet'
    sc = ReadNumSidecar.open(sidecar)
    assert len(sc) == 100
    headers = []
    for chunk in chunks:
        with gzip.open(chunk, 'rt') as fh:
            while True:
                h = fh.readline().rstrip('\n')
                if not h:
                    break
                headers.append(h)
                seq = fh.readline().rstrip('\n')
                fh.readline()
                qual = fh.readline().rstrip('\n')
                rn = int(h.split('RN:i:')[1].split('\t')[0])
                assert sc.verify(rn, seq, qual).ok
    assert len(headers) == 100
    assert all('\tRN:i:' in h for h in headers)
    tagged = next(h for h in headers if h.startswith('@read1\t'))
    assert '\tRN:i:1\tXA:Z:foo\tXC:i:42' in tagged
    assert sc.lookup(1).fastq_comment == 'XA:Z:foo\tXC:i:42'
