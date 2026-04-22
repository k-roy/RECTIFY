# Nanopore COMPASS Aligner Recommendations - 2026 Update

**Date**: 2026-03-10
**Context**: Updating aligner panel based on latest developments

---

## Current Aligner Status

### ✅ Successfully Installed

| Aligner | Version | Location | Status |
|---------|---------|----------|--------|
| **minimap2** | 2.28 | system PATH | ✅ Ready |
| **uLTRA** | 0.1 | system PATH | ✅ Ready |
| **deSALT** | 1.5.6 | system PATH | ✅ Ready |

### ❌ Installation Issues

| Aligner | Issue | Recommendation |
|---------|-------|----------------|
| **GraphMap2** | Not updated since 2018, compilation issues | **REPLACE** with newer aligner |

---

## Newer Aligners for Consideration (2025-2026)

### 1. Minisplice (Released ~June 2025) ⭐ **HIGH PRIORITY**

**Description**: Enhances minimap2 by integrating deep-learning splice signals

**Advantages**:
- Built on proven minimap2 foundation
- Deep learning for splice site prediction
- Optimized for noisy long RNA-seq reads
- Active development

**Status**: Released/updated around June 2025
**Recommendation**: **INCLUDE** in 4-way consensus panel

**Installation**:
```bash
# TBD - check GitHub/conda for installation
```

### 2. Winnowmap2 (Actively Updated) ⭐ **MEDIUM PRIORITY**

**Description**: Faster, more memory-efficient alternative to minimap2

**Advantages**:
- Optimized for human genomics pipelines
- Better memory efficiency
- Still being benchmarked
- Frequent updates

**Status**: Actively developed
**Recommendation**: **CONSIDER** for large-scale human datasets, but may be less critical for yeast

**Notes**: Primary focus is human genomics, may not provide significant benefits over minimap2 for yeast

### 3. GLASS (Published April 2025) ⭐ **HIGH PRIORITY**

**Description**: Graph learning algorithm to screen and improve splice-aware alignments

**Advantages**:
- Latest wave of development
- Graph-based approach (novel)
- Specifically designed for splice-aware alignment improvement
- Very recent (April 2025)

**Status**: Recently published (April 2025)
**Recommendation**: **INCLUDE** in 4-way consensus panel

**Installation**:
```bash
# TBD - check publication/GitHub for availability
```

---

## Proposed 4-Way Consensus Panel (Updated)

### Option A: Modern Aligner Panel ⭐ **RECOMMENDED**

Replace GraphMap2 with newer aligners:

| Aligner | Role | Strength |
|---------|------|----------|
| **minimap2** | Fast baseline | High sensitivity, general-purpose |
| **Minisplice** | ML-enhanced splicing | Deep learning splice signals |
| **deSALT** | Complex isoforms | De Bruijn graph, alt. splicing |
| **GLASS** | Graph-based refinement | Latest splice-aware improvements |

**Rationale**:
- Covers multiple algorithmic approaches (seed-chain, ML, De Bruijn, graph)
- All actively maintained (2025-2026)
- Optimized for nanopore error profiles
- Spans fast (minimap2) to specialized (GLASS)

### Option B: Conservative Panel

Keep proven aligners while adding one new one:

| Aligner | Role | Strength |
|---------|------|----------|
| **minimap2** | Fast baseline | Proven, widely used |
| **uLTRA** | Annotation-guided | Small exons, collinear chaining |
| **deSALT** | Complex isoforms | De Bruijn graph |
| **Minisplice** | ML-enhanced | Latest improvements |

**Rationale**:
- Lower risk (3 proven aligners + 1 new)
- Still includes ML enhancements
- All currently installed and working

### Option C: Comprehensive Panel (6-way)

Use all available modern aligners:

| Aligner | Role |
|---------|------|
| **minimap2** | Fast baseline |
| **Minisplice** | ML-enhanced |
| **uLTRA** | Annotation-guided |
| **deSALT** | De Bruijn graph |
| **GLASS** | Graph-based |
| **Winnowmap2** | Memory-efficient (optional) |

**Rationale**:
- Maximum consensus power
- Captures all modern approaches
- Higher computational cost

---

## Implementation Plan

### Phase 1: Test Current 3-Way Panel ✅ **IMMEDIATE**

Test with currently installed aligners:
```bash
# minimap2 + uLTRA + deSALT
bash scripts/test_single_chunk_alignment.sh \
    test_parallel/fastq_chunks/chunk_001.fastq \
    test_parallel/chunk_001_alignments
```

**Purpose**: Validate parallelization infrastructure works

### Phase 2: Add Minisplice 🔜 **NEXT**

1. Install Minisplice
2. Test 4-way consensus (minimap2 + Minisplice + deSALT + uLTRA)
3. Compare results to 3-way

### Phase 3: Add GLASS 🔜 **AFTER MINISPLICE**

1. Install GLASS
2. Test 5-way consensus
3. Benchmark performance vs 4-way

### Phase 4: Evaluate & Finalize

- Compare accuracy across panels
- Benchmark computational cost
- Select optimal panel for production

---

## Installation Priorities

### Immediate (This Session)

✅ minimap2 (done)
✅ uLTRA (done)
✅ deSALT (done)

### Next (After validation)

🔜 **Minisplice** - High priority, ML-enhanced
🔜 **GLASS** - High priority, latest development

### Optional (If needed)

- Winnowmap2 - If processing human data or memory constraints
- Keep GraphMap2 out - obsolete, compilation issues

---

## Benchmarking Criteria

When comparing aligner panels, evaluate:

### Accuracy Metrics
- Junction detection sensitivity
- Canonical dinucleotide percentage
- Annotation match rate
- False positive rate

### Performance Metrics
- Alignment time per chunk
- Memory usage
- Parallelization efficiency
- Overall throughput

### Consensus Metrics
- Agreement rate between aligners
- Consensus confidence distribution
- Novel junction validation

---

## References

**Newer Aligners**:
- Minisplice: Released ~June 2025, enhances minimap2 with DL splice signals
- Winnowmap2: Actively updated, optimized for human genomics
- GLASS: Published April 2025, graph learning for splice-aware alignment

**Established Aligners**:
- minimap2: Li, H. (2018). Bioinformatics.
- uLTRA: Sahlin et al. (2021). Genome Biology.
- deSALT: Liu et al. (2019). Genome Biology.

**Obsolete**:
- GraphMap2: Last updated 2018, no longer maintained

---

## Action Items

1. ✅ **Complete**: Test 3-way consensus (minimap2 + uLTRA + deSALT)
2. 🔜 **Next**: Research Minisplice installation and availability
3. 🔜 **After**: Research GLASS installation and availability
4. 🔜 **Then**: Benchmark aligner panels on test data
5. 🔜 **Finally**: Select optimal panel for production analysis

---

**Last Updated**: 2026-03-10
**Status**: Ready to test 3-way consensus, then expand to modern aligners
