"""Prototype: vectorized DP for _score_junction, table=None. Verify byte-identical vs the per-k reference.
Strategy: compute ALL t1(k) in one DP and ALL t2(k) in one DP, then min_k(t1+t2). Harness is the arbiter."""
import sys, os, random
sys.path.insert(0, '/Users/kevinroy/work/rectify/.claude/worktrees/agent-a25a2c1e784ad37dc')
import rectify.core.splice.junction_scoring as J
from rectify.core.splice.hp_penalty import _precompute_del_costs

INS, DEL_N, DEL_HP, SUB, HP_MIN = 1.25, 1.0, 0.5, 1.0, 4
INF = float('inf')

def score_hp_anchored_ref(query, ref, del_costs):
    """Faithful copy of the table=None DP in _score_hp_anchored (pure-python branch)."""
    Q, R = len(query), len(ref)
    if Q == 0: return 0.0
    if R == 0: return Q * INS
    prev = [INF]*(R+1); prev[0]=0.0
    for j in range(1, R+1): prev[j]=prev[j-1]+del_costs[j-1]
    for i in range(1, Q+1):
        curr=[INF]*(R+1); curr[0]=i*INS; qi=query[i-1].upper()
        for j in range(1, R+1):
            cost_sub = 0.0 if qi==ref[j-1].upper() else SUB
            curr[j]=min(prev[j-1]+cost_sub, prev[j]+INS, curr[j-1]+del_costs[j-1])
        prev=curr
    return min(prev)

def all_t1(rescue, ref, del_costs):
    """t1(k) = score_hp_anchored(rescue[k:], ref, del_costs) for k=0..L-1, as a vector.
    Build by running the DP once over the FULL rescue but recording, for each query-start k,
    the min-over-ref-end. Trick: process query rows; a DP started at row k = the DP for rescue[k:].
    We compute all starts by a backward pass: G[k] = min over j of D_k[L][j].
    Do it with L independent inits fused: maintain, for each start k, its own prev row is expensive.
    -> Use the reversal: t1(k) over reversed seqs = a query-PREFIX DP. Implement + verify."""
    L=len(rescue); R=len(ref)
    rr=rescue[::-1]; rf=ref[::-1]; dc=del_costs[::-1]  # reversed del costs align to reversed ref
    # forward DP over rr vs rf: query global (rr prefix consumed), free ref LEFT prefix (=orig right suffix),
    # answer at ref RIGHT end R (=orig left, the anchor). D[m][j], m=0..L, j=0..R.
    # init: D[0][j]=0 (free ref left prefix skip), D[m][0]=m*INS
    prev=[0.0]*(R+1)
    t1=[None]*L
    # m=0 row done (prev). t1(k) uses query-prefix length m=L-k -> t1(L-0)=... we need t1(k)=D[L-k][R]
    # record after each row
    D_at_R=[prev[R]]  # m=0
    for m in range(1, L+1):
        curr=[INF]*(R+1); curr[0]=m*INS; qm=rr[m-1].upper()
        for j in range(1, R+1):
            cost_sub=0.0 if qm==rf[j-1].upper() else SUB
            curr[j]=min(prev[j-1]+cost_sub, prev[j]+INS, curr[j-1]+dc[j-1])
        prev=curr; D_at_R.append(prev[R])
    for k in range(L):
        t1[k]=D_at_R[L-k]
    return t1

# ---- verify all_t1 against per-k reference ----
random.seed(3); mism=0; first=None
for t in range(3000):
    L=random.randint(1,30); R=random.randint(1,50)
    rescue="".join(random.choice("ACGT") for _ in range(L))
    ref="".join(random.choice("ACGT") for _ in range(R))
    # inject HP runs
    if random.random()<0.5:
        p=random.randint(0,max(0,R-8)); ref=ref[:p]+random.choice("ACGT")*random.randint(4,8)+ref[p+0:]
        ref=ref[:R+8]; R=len(ref)
    dc=_precompute_del_costs(ref, None, 0, False, None, HP_MIN, DEL_N, DEL_HP)
    ref_v=all_t1(rescue, ref, dc)
    for k in range(L):
        got=ref_v[k]; want=score_hp_anchored_ref(rescue[k:], ref, dc)
        if abs(got-want)>1e-12:
            mism+=1
            if first is None: first=(rescue,ref,k,got,want)
            break
print(f"all_t1 vs per-k reference: mismatches={mism}/3000")
if first: print("first mismatch:", first)
else: print("all_t1 BYTE-IDENTICAL")

# ============ FULL _score_junction replacement (table=None) via all_t1 twice ============
def score_junction_fast(query, q_split, intron_end, genome_seq):
    _MAX_RESCUE=30
    rescue=query[q_split:q_split+_MAX_RESCUE].upper()
    L=len(rescue)
    if L==0: return 0.0
    gs=len(genome_seq)
    _BUF=max(L+20,50)
    ref_exon2=genome_seq[intron_end:min(gs,intron_end+_BUF)].upper()
    ref_intron_end=genome_seq[max(0,intron_end-_BUF):intron_end].upper()
    ref_intron_end_rev=ref_intron_end[::-1]
    dc_fwd=_precompute_del_costs(ref_exon2, genome_seq, intron_end, False, None)
    dc_rev=_precompute_del_costs(ref_intron_end_rev, genome_seq, intron_end-1, True, None)
    t1=all_t1(rescue, ref_exon2, dc_fwd)                       # t1[k]=score(rescue[k:], ref_exon2)
    V2=all_t1(rescue[::-1], ref_intron_end_rev, dc_rev)        # V2[m]=score(rev(rescue)[m:], ...)
    best=INF
    for k in range(L):
        t2 = 0.0 if k==0 else V2[L-k]
        s=t1[k]+t2
        if s<best: best=s
        if best==0.0: break
    return best

# ---- verify full score_junction_fast vs the REAL J._score_junction (table=None) ----
import random as _r; _r.seed(11); mism=0; first=None; N=8000
for t in range(N):
    gl=_r.randint(40,140); genome="".join(_r.choice("ACGT") for _ in range(gl))
    # inject HP run(s)
    for _ in range(_r.randint(0,2)):
        p=_r.randint(0,max(0,gl-12)); genome=genome[:p]+_r.choice("ACGT")*_r.randint(4,12)+genome[p:]
    gl=len(genome)
    ie=_r.randint(20,gl-20)
    rl=_r.randint(1,30); st=max(0,ie-_r.randint(0,rl))
    q=genome[st:st+rl]
    q="".join(c if _r.random()>0.15 else _r.choice("ACGT") for c in q)  # errors
    want=J._score_junction(q,0,0,ie,genome,hp_pen=0.25,W=15,max_slide=0,penalty_table=None)[0]
    got=score_junction_fast(q,0,ie,genome)
    if abs(got-want)>1e-12:
        mism+=1
        if first is None: first=(q,ie,got,want)
print(f"\nscore_junction_fast vs REAL _score_junction (table=None): {N-mism}/{N} EXACT; mismatches={mism}")
if first: print("first mismatch:", first)
else: print("FULL BYTE-IDENTICAL on the table-free path (MECH1 handled by k in [0,L))")

# ---- speedup: old per-k _score_junction vs fast, on a representative batch ----
import time as _t
_r.seed(5); batch=[]
for _ in range(400):
    gl=_r.randint(60,140); genome="".join(_r.choice("ACGT") for _ in range(gl))
    p=_r.randint(0,gl-12); genome=genome[:p]+_r.choice("ACGT")*_r.randint(4,12)+genome[p:]; gl=len(genome)
    ie=_r.randint(20,gl-20); q=genome[max(0,ie-12):ie+18]
    q="".join(c if _r.random()>0.15 else _r.choice("ACGT") for c in q)
    batch.append((q,ie,genome))
t0=_t.time()
for q,ie,g in batch: J._score_junction(q,0,0,ie,g,hp_pen=0.25,W=15,max_slide=0,penalty_table=None)
told=_t.time()-t0
t0=_t.time()
for q,ie,g in batch: score_junction_fast(q,0,ie,g)
tnew=_t.time()-t0
print(f"\nSPEEDUP (400 candidate scorings, table=None): old {told*1000:.0f}ms ({told/len(batch)*1000:.2f}ms/call) "
      f"-> fast {tnew*1000:.0f}ms ({tnew/len(batch)*1000:.2f}ms/call) = {told/tnew:.1f}x")
print(f"DP-pass count: old ~2L(<=60)/candidate -> fast 2/candidate")
