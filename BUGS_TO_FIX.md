# RECTIFY Bugs to Fix

## Bug 1: SettingWithCopyWarning in clustering.py

**Location:** `rectify/core/analyze/clustering.py:473`

**Warning:**
```
SettingWithCopyWarning: A value is trying to be set on a copy of a slice from a DataFrame.
Try using .loc[row_indexer,col_indexer] = value instead
```

**Current code:**
```python
assigned['_effective_count'] = assigned[count_col]
```

**Fix:** Use `.loc[]` or ensure `assigned` is a copy:
```python
assigned = assigned.copy()
assigned['_effective_count'] = assigned[count_col]
```

---

## Bug 2: Sample column naming inconsistency

**Issue:** Default sample column is "sample" but data often uses "replicate"

**Current behavior:** Fails with `ValueError: Missing required columns: ['sample']`

**Workaround:** Use `--sample-column replicate`

**Suggested fix:** Either:
1. Auto-detect sample column from common names: ["sample", "replicate", "sample_id", "sample_name"]
2. Better document the requirement in help text

---

## Enhancement: GFF3 parser improvements

**Fixed in this session:** The `_parse_gtf()` function now handles both GTF and GFF3 formats:
- GFF3: Extracts `ID` → gene_id, `gene` → common name, `Name` → systematic name
- GTF: Extracts `gene_id` and `gene_name`

**Location:** `rectify/core/analyze_command.py:707-770`

---

## Last Updated: 2026-03-18
