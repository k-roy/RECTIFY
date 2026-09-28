"""ISSUE-054: scored Case3 must obey the realized junction-indel contract."""
from array import array
import copy
import json
from pathlib import Path

import pysam

from rectify.core.align.local_aligner import cigar_str_to_ops
from rectify.core.bam import bam_writer as bw
from rectify.core.bam.bam_processor import correct_read_3prime, _build_pool_chrom_index
from rectify.core.bam.output import CORRECTION_TSV_HEADER, correction_result_to_tsv_row
from rectify.core.splice import splice_aware_5prime as s
from rectify.utils.genome import register_genome_contigs

RC = str.maketrans('ACGTN', 'TGCAN')
TAILS = ('G', 'GC', 'GCA', 'GCAT', '', 'AAA')


def inputs():
    base = 'T' * 20 + 'A' + 'ACGTTGCATGCAGTCCATG' + 'GTAAGT' + 'N' * 92 + 'AG' + 'C' * 100
    genomes = {}
    for reverse in (False, True):
        for hp in (False, True):
            g = base[:36] + 'AAAA' + base[40:] if hp else base
            if reverse:
                g = g.translate(RC)[::-1]
            genomes[('hp' if hp else 'plain') + ('M' if reverse else 'P')] = g
    register_genome_contigs(list(genomes))
    header = pysam.AlignmentHeader.from_references(list(genomes), [len(g) for g in genomes.values()])
    records, annotation, descriptions = [], set(), {}
    for reverse in (False, True):
        for tail in TAILS:
            hp = tail == 'AAA'
            chrom = ('hp' if hp else 'plain') + ('M' if reverse else 'P')
            g = base[:36] + 'AAAA' + base[40:] if hp else base
            q = g[28:40] + tail + g[140:200]
            ops = [(4, 12)] + ([(1, len(tail))] if tail else []) + [(0, 60)]
            start, junction = 140, (chrom, 40, 140)
            if reverse:
                q, ops, start, junction = q.translate(RC)[::-1], ops[::-1], 40, (chrom, 100, 200)
            annotation.add(junction)
            kind = 'hp' if hp else ('clean' if not tail else 'non_hp' + str(len(tail)))
            for encoding in ('literal', 'equals'):
                for hard in (0, 5):
                    r = pysam.AlignedSegment(header)
                    r.query_name = f'{kind}_{"minus" if reverse else "plus"}_{encoding}_H{hard}'
                    r.flag = 16 if reverse else 0
                    r.reference_id = header.get_tid(chrom)
                    r.reference_start, r.mapping_quality = start, 60
                    r.cigartuples = (ops + [(5, hard)] if reverse else [(5, hard)] + ops) if hard else ops
                    r.query_sequence = q
                    if encoding == 'equals':
                        encoded = list(q)
                        for qp, rp in r.get_aligned_pairs(matches_only=True):
                            if qp is not None and rp is not None and q[qp] == genomes[chrom][rp]:
                                encoded[qp] = '='
                        r.query_sequence = ''.join(encoded)
                    r.query_qualities = array('B', [35] * len(q))
                    r.set_tag('ZZ', 'placement_control')
                    r.set_tag('XB', 'unrelated_strand_tag')
                    r.set_tag('ZA', array('h', [1, -2, 300]))
                    records.append(r)
                    descriptions[r.query_name] = dict(kind=kind, strand='-' if reverse else '+', encoding=encoding,
                        original_hard=hard, junction=junction,
                        prohibited=kind.startswith('non_hp') and (reverse or len(tail) > 1))
    return genomes, records, annotation, descriptions


def snap(read, genomes):
    encoded = read.to_string()
    r = copy.deepcopy(read)
    bw._decode_eq_seq_inplace(r, genomes)
    qp = 0
    rp = r.reference_start
    junctions = []
    for op, length in r.cigartuples:
        if op == 3:
            junctions.append((rp, rp + length))
        if op in (0, 2, 3, 7, 8):
            rp += length
    return dict(sam=encoded, cigar=r.cigarstring, start=r.reference_start, end=r.reference_end,
        sequence=r.query_sequence, qualities=list(r.query_qualities or []),
        tags=r.get_tags(), query_map=r.get_aligned_pairs(), junctions=junctions,
        molecule_extent=len(r.query_sequence) + sum(n for op, n in r.cigartuples if op == 5))


def run_roundtrip(directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    genomes, reads, annotation, descriptions = inputs()
    pool_index = _build_pool_chrom_index(annotation)
    results, decisions = [], []
    for read in reads:
        decoded = copy.deepcopy(read)
        bw._decode_eq_seq_inplace(decoded, genomes)
        detail = descriptions[read.query_name]
        raw = s.rescue_3ss_truncation(decoded, genomes, annotation, detail['strand'], annotated_junctions=annotation)
        # Pin the scored Case3 path before the independent terminal-peel
        # fallback gets another opportunity to propose a different block.
        scored = s.rescue_3ss_truncation(decoded, genomes, annotation, detail['strand'],
                                       annotated_junctions=annotation, terminal_peel=False)
        decision = correct_read_3prime(copy.deepcopy(read), genomes, annotated_junctions=annotation,
                                      pool_chrom_index=pool_index)[0]
        ops = cigar_str_to_ops(decision.get('five_prime_exon_cigar', ''))
        j = detail['junction']
        gate = s._junction_adjacent_indel_refusal(ops, detail['strand'], genomes[j[0]], j[1], j[2])
        record = dict(read_id=read.query_name, case=detail, before=snap(read, genomes),
                      raw_rescue=raw, scored_without_terminal_peel=scored,
                      decision=decision, emitted_block_gate=gate)
        results.append(record)
        decisions.append(decision)
        (directory / 'decisions.json').write_text(json.dumps(results, indent=2, default=str) + '\n')
    stock, tsv = directory / 'stock.bam', directory / 'corrections.tsv'
    with pysam.AlignmentFile(str(stock), 'wb', header=reads[0].header) as handle:
        for read in reads:
            handle.write(read)
    with tsv.open('w') as handle:
        handle.write('\t'.join(CORRECTION_TSV_HEADER) + '\n')
        for row in decisions:
            handle.write('\t'.join(correction_result_to_tsv_row(row)) + '\n')
    loaded = bw._load_corrections_from_tsv(str(tsv))
    bw.write_corrected_bam(str(stock), str(tsv), str(directory / 'hard.bam'), genomes)
    bw.write_softclipped_bam(str(stock), str(tsv), str(directory / 'soft.bam'), genomes)
    bw.write_dual_bam(str(stock), str(tsv), str(directory / 'dual_hard.bam'), str(directory / 'dual_soft.bam'), genomes)
    outputs = {}
    for arm in ('hard', 'soft', 'dual_hard', 'dual_soft'):
        with pysam.AlignmentFile(str(directory / (arm + '.bam')), 'rb') as handle:
            outputs[arm] = {read.query_name: snap(read, genomes) for read in handle}
    for record in results:
        name = record['read_id']
        record['loaded_correction'] = loaded[name]
        record['outputs'] = {arm: output[name] for arm, output in outputs.items()}
    (directory / 'results.json').write_text(json.dumps(results, indent=2, default=str) + '\n')
    (directory / 'genome.fa').write_text(''.join('>' + c + '\n' + g + '\n' for c, g in genomes.items()))
    return results


def assert_fixed(results):
    assert len(results) == 48
    # the family must still contain blocks the NOVEL-landing predicate refuses
    assert any(r['case']['prohibited'] and r['emitted_block_gate'] == s.JUNCTION_INDEL_REFUSAL
               for r in results)
    for record in results:
        prohibited = record['case']['prohibited']
        decision = record['decision']
        if prohibited:
            # ISSUE-083 (Kevin 2026-09-21): every junction in this family is ANNOTATED, and at an
            # annotated landing an I/D beside the N is "not a deal breaker" — the block is drawn.
            # The predicate still names this exact shape for a NOVEL landing (its default), which
            # is what `emitted_block_gate` records; that is the half of ISSUE-054 that stands.
            assert record['scored_without_terminal_peel']['rescued'], record['read_id']
            assert decision['five_prime_rescued'], record['read_id']
            assert not decision['five_prime_rescue_refused'], record['read_id']
            # '' when the insertion could be slid off the N through identical bases.
            assert record['emitted_block_gate'] in ('', s.JUNCTION_INDEL_REFUSAL), record['read_id']
        else:
            assert record['raw_rescue']['rescued'], record['read_id']
            assert decision['five_prime_rescued'], record['read_id']
            assert not record['emitted_block_gate'], record['read_id']
        b = record['before']
        for arm, out in record['outputs'].items():
            assert out['sequence'] == b['sequence'] and out['qualities'] == b['qualities']
            assert out['molecule_extent'] == b['molecule_extent']
            tags = dict(out['tags'])
            for tag in ('ZZ', 'XB', 'ZA'):
                assert tags[tag] == dict(b['tags'])[tag]
            assert sorted(out['junctions']) == sorted(tuple(x) for x in decision['junctions'])
            assert out['junctions'] == [tuple(record['case']['junction'][1:])]
        assert record['outputs']['hard']['sam'] == record['outputs']['dual_hard']['sam']
        assert record['outputs']['soft']['sam'] == record['outputs']['dual_soft']['sam']


def test_scored_case3_final_indel_gate_all_writer_modes(tmp_path):
    assert_fixed(run_roundtrip(tmp_path))
