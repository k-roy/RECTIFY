"""Shared indexed reference access for alignment post-passes and consensus."""
from contextlib import contextmanager
import gzip
import hashlib
import json
import logging
from pathlib import Path
import shutil
import tempfile

logger = logging.getLogger(__name__)


@contextmanager
def open_alignment_reference(genome_path, cache_dir):
    """Open FASTA/BGZF, normalizing ordinary gzip once with bounded memory.

    pysam cannot index ordinary gzip. The alignment command already supported
    decompressing it for consensus; post-passes use this same path. Original
    reference bytes are retained, and a source-bound indexed cache is reused.
    """
    import pysam
    source = Path(genome_path)
    try:
        reference = pysam.FastaFile(str(source))
    except OSError:
        if not str(source).endswith('.gz') or not source.is_file():
            raise
        # Do not disguise a corrupt BAM/reference or arbitrary missing path as
        # a format fallback: gzip.open below must successfully decode the file.
        cache = Path(cache_dir)
        cache.mkdir(parents=True, exist_ok=True)
        identity = dict(path=str(source.resolve()), size=source.stat().st_size,
                        mtime_ns=source.stat().st_mtime_ns)
        key = hashlib.sha256(str(source.resolve()).encode()).hexdigest()[:12]
        target = cache / f'{source.name}.{key}.reference.fa'
        receipt = Path(str(target) + '.source.json')
        try:
            reusable = (json.loads(receipt.read_text()) == identity
                        and target.is_file() and Path(str(target) + '.fai').is_file())
        except (OSError, ValueError):
            reusable = False
        if not reusable:
            logger.info('Normalizing gzip reference for indexed access: %s', source)
            with tempfile.TemporaryDirectory(prefix='.reference-', dir=cache) as tmp:
                raw = Path(tmp) / 'reference.fa'
                with gzip.open(source, 'rb') as fin, raw.open('wb') as fout:
                    shutil.copyfileobj(fin, fout, length=1024 * 1024)
                pysam.faidx(str(raw))
                raw.replace(target)
                Path(str(raw) + '.fai').replace(Path(str(target) + '.fai'))
                marker = Path(tmp) / 'source.json'
                marker.write_text(json.dumps(identity) + '\n')
                marker.replace(receipt)
        reference = pysam.FastaFile(str(target))
    try:
        yield reference
    finally:
        reference.close()
