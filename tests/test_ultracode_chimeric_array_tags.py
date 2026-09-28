"""Typed array metadata survives actual chimeric and fallback BAM products."""
from array import array

import pysam
import pytest

from tests.test_ultracode_chimeric_boundary_ownership import make_family
from tests.test_ultracode_chimeric_hardclips import build, with_hardclips


ARRAYS = {
    'b': [-128, -1, 127], 'B': [0, 128, 255],
    'h': [-32768, -1, 32767], 'H': [0, 32768, 65535],
    'i': [-2147483648, -1, 2147483647], 'I': [0, 2147483648, 4294967295],
    'f': [-1.25, 0.0, 3.5],
}


@pytest.mark.parametrize('subtype', list(ARRAYS))
@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('fallback', [False, True])
def test_array_subtype_roundtrip(subtype, reverse, fallback, tmp_path):
    arms, _, genome, annotation = make_family(reverse, soft=True)
    arms = {name: with_hardclips(read, 5, 9) for name, read in arms.items()}
    for read in arms.values():
        read.set_tag('ZZ', array(subtype, ARRAYS[subtype]))
        read.set_tag('zs', 'retained')
        read.set_tag('zi', -129, value_type='s')
        read.set_tag('NM', 17)
        read.set_tag('AS', 3)
        read.set_tag('SA', 'stale-placement')
        read.set_tag('MD', '1A1')
    if fallback:
        # Different valid retained intervals prohibit cross-frame stitching.
        ops = list(arms['uLTRA'].cigartuples)
        ops[0], ops[-1] = (5, 9), (5, 5)
        arms['uLTRA'].cigartuples = ops
    before = {name: read.to_string() for name, read in arms.items()}
    result, out = build(arms, genome, annotation)
    assert result.is_fallback == fallback
    assert out.get_tag('ZZ').typecode == subtype
    assert out.get_tag('ZZ').tolist() == array(subtype, ARRAYS[subtype]).tolist()
    assert out.get_tag('zs') == 'retained'
    assert out.get_tag('zi', with_value_type=True) == (-129, 's')
    assert all(not out.has_tag(tag) for tag in ('NM', 'AS', 'SA', 'MD'))
    path = tmp_path / 'roundtrip.bam'
    with pysam.AlignmentFile(path, 'wb', header=out.header) as bam:
        bam.write(out)
    with pysam.AlignmentFile(path, 'rb') as bam:
        reread = next(bam)
    assert reread.to_string() == out.to_string()
    assert reread.get_tag('ZZ').typecode == subtype
    assert before == {name: read.to_string() for name, read in arms.items()}
