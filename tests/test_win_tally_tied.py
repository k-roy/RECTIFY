"""The consensus win tally must credit ALL aligners tied for the top junction
score, not just the single tiebreak winner.

On a score tie, select.py picks one ``best_aligner`` via an arbitrary ASCII-name
fallback (added in ba06e4f4 purely for determinism). Per Kevin's directive the
tally (``stats['by_aligner']`` -> ``aligner_wins_*``) must keep every tied top
aligner. The write-selection is unchanged — only the accounting.

These are pure unit tests of ``_credit_tied_aligners`` (the crediting the batch
loop performs), so they don't need BAM fixtures.
"""
from collections import defaultdict

from rectify.core.consensus.consensus import ConsensusResult, _credit_tied_aligners


def _stats():
    return {'by_aligner': defaultdict(int)}


def _result(best, tied):
    return ConsensusResult(
        read_id='r',
        best_aligner=best,
        best_alignment=None,
        aligners_compared=list(tied) or [best],
        tied_aligners=list(tied),
    )


def test_all_tied_top_aligners_credited_not_just_winner():
    # 3-way score tie. select.py's ASCII-name fallback picks 'magicblast'
    # (max name) as best_aligner; the tally must ALSO credit the two losers.
    stats = _stats()
    _credit_tied_aligners(stats, _result(best='magicblast',
                                         tied=['bbmap', 'gsnap', 'magicblast']))
    assert stats['by_aligner']['bbmap'] == 1        # tiebreak LOSER — still credited
    assert stats['by_aligner']['gsnap'] == 1        # tiebreak LOSER — still credited
    assert stats['by_aligner']['magicblast'] == 1   # the tiebreak winner
    assert sum(stats['by_aligner'].values()) == 3   # all three, not one


def test_unique_winner_credits_only_that_aligner():
    stats = _stats()
    _credit_tied_aligners(stats, _result(best='minimap2', tied=['minimap2']))
    assert dict(stats['by_aligner']) == {'minimap2': 1}


def test_empty_tied_falls_back_to_best_aligner():
    # Degenerate empty-alignments result (tied=[], best='none'): must stay
    # byte-identical to the old single-credit behavior.
    stats = _stats()
    _credit_tied_aligners(stats, ConsensusResult(
        read_id='', best_aligner='none', best_alignment=None,
        aligners_compared=[], tied_aligners=[]))
    assert dict(stats['by_aligner']) == {'none': 1}


def test_credits_accumulate_across_reads():
    stats = _stats()
    _credit_tied_aligners(stats, _result('gsnap', ['bbmap', 'gsnap']))  # tie
    _credit_tied_aligners(stats, _result('bbmap', ['bbmap']))           # unique
    assert stats['by_aligner']['bbmap'] == 2
    assert stats['by_aligner']['gsnap'] == 1
