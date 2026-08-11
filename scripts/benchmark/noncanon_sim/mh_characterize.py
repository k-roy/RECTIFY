"""Move-aware microhomology: for an acceptor shift ne->je (k=je-ne), the drift is enabled iff the
k-mer absorbed into the intron (genome[ne:je]) near-matches the new exon2 k-mer (genome[je:je+k]).
That is EXACTLY the spike-in's installed repeat genome[A:A+k]~genome[A+k:A+2k]. Measure fraction matched."""
import os, csv, statistics as st
def load_fa(p):
    g={}; n=None
    for l in open(p):
        if l.startswith(">"): n=l[1:].split()[0]; g[n]=[]
        else: g[n].append(l.strip())
    return {k:"".join(v).upper() for k,v in g.items()}
def move_mh_frac(seq, ne, je):
    """fraction of the |k|-mer that repeats at the drift distance (the drift-enabling microhomology)."""
    k=je-ne
    if k>0:   # acceptor shifted DOWNSTREAM: absorbed genome[ne:je] vs new-exon genome[je:je+k]
        a=seq[ne:je]; b=seq[je:je+k]
    else:     # shifted UPSTREAM: symmetric
        k=-k; a=seq[je:ne]; b=seq[je-k:je]
    if len(a)<k or len(b)<k or k==0: return 0.0, k
    return sum(1 for x,y in zip(a,b) if x==y)/k, k

g=load_fa("fab_sweep/sim_ref.fa")
print("=== FAB: microhomology of the DRIFT move (true acceptor -> drift acceptor) ===")
fab=[]
for r in csv.DictReader(open("fab_sweep/fab_contig_truth.tsv"), delimiter="\t"):
    c=r["chrom"]; A=int(r["acceptor"]); Ad=int(r["drift_acceptor"]); k=int(r["drift_dist"])
    if c in g:
        f,kk=move_mh_frac(g[c], A, Ad); fab.append((k,f))
        print(f"  {c}: move {A}->{Ad} (k={k})  microhomology_frac={f:.2f}")
gr=load_fa("mix_r3b_out/sim_ref.fa")
print("\n=== R3 (genuine non-canonical): microhomology of the REAL move (canonical decoy -> true non-canonical) ===")
r3=[]; seen=set()
# panel_truth has decoy_acceptor (canonical) + true_acceptor (non-canonical)
for r in csv.DictReader(open("mix_r3b_out/panel_truth.tsv"), delimiter="\t"):
    if r.get("motif_rung")!="R3": continue
    c=r["chrom"]; A=int(r["true_acceptor"]); dec=r.get("decoy_acceptor","")
    if not dec or c not in gr or (c,A) in seen: continue
    seen.add((c,A)); f,kk=move_mh_frac(gr[c], int(dec), A); r3.append(f)
    print(f"  {c}: move {dec}->{A} (k={kk})  microhomology_frac={f:.2f}")
if fab and r3:
    fabf=[f for k,f in fab]
    print(f"\n★ SEPARATION: FAB-drift mh_frac mean={st.mean(fabf):.2f} (min {min(fabf):.2f}) vs R3-real mean={st.mean(r3):.2f} (max {max(r3):.2f})")
    print("=> a microhomology-frac threshold between them separates fabricated drift from real discovery.")
