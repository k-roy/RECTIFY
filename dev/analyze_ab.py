import sys

def load(path):
    rows = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            rows[p[idx["read_id"]]] = p
        return rows, idx

B, idx = load(sys.argv[1])
P, _ = load(sys.argv[2])
c5 = idx["five_prime_position"]; cr = idx["five_prime_rescued"]; cj = idx["junctions"]; cg = idx["gene_id"]

new_rescue = lost_rescue = diff_junc = diff_pos_only = same = 0
examples = []
for rid, b in B.items():
    p = P.get(rid)
    if p is None:
        continue
    changed = (b[c5] != p[c5]) or (b[cj] != p[cj])
    if not changed:
        same += 1
        continue
    if b[cr] == "0" and p[cr] == "1":
        new_rescue += 1; tag = "NEW_RESCUE"
    elif b[cr] == "1" and p[cr] == "0":
        lost_rescue += 1; tag = "LOST_RESCUE"
    elif b[cj] != p[cj]:
        diff_junc += 1; tag = "DIFF_JUNCTION"
    else:
        diff_pos_only += 1; tag = "DIFF_POS"
    if len(examples) < 12:
        examples.append((tag, rid, b[cg], b[c5], b[cj], "->", p[c5], p[cj]))

def jlen(jstr):
    # junctions like "35749-37076" (may be empty or multiple comma-sep)
    if not jstr:
        return 0
    best = 0
    for seg in jstr.split(","):
        if "-" in seg:
            a, b2 = seg.split("-")[:2]
            try:
                best = max(best, abs(int(b2) - int(a)))
            except ValueError:
                pass
    return best

print(f"total_common={len(B)} same={same}")
print(f"post_NEW_rescue={new_rescue} post_LOST_rescue={lost_rescue} DIFF_junction={diff_junc} DIFF_pos_only={diff_pos_only}")
# intron-length distribution of the POST junction on changed reads
post_lens = sorted(jlen(P[rid][cj]) for rid, b in B.items()
                   if rid in P and ((b[c5]!=P[rid][c5]) or (b[cj]!=P[rid][cj])))
import statistics
if post_lens:
    print(f"post_junction_len on changed reads: min={post_lens[0]} median={statistics.median(post_lens)} max={post_lens[-1]}")
    print(f"  changed reads with post junction >1500bp: {sum(1 for x in post_lens if x>1500)} / {len(post_lens)}")
    print(f"  changed reads with post junction >2000bp: {sum(1 for x in post_lens if x>2000)} / {len(post_lens)}")
# distribution of changed gene_ids
from collections import Counter
genes = Counter(b[cg] for rid, b in B.items() if rid in P and ((b[c5]!=P[rid][c5]) or (b[cj]!=P[rid][cj])))
print("changed_by_gene:", genes.most_common(10))
print("examples (tag, read, gene, base_5p, base_junc -> post_5p, post_junc):")
for e in examples:
    print("  ", e)
