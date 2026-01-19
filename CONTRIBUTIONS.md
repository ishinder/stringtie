# Contributions to StringTie3

This fork contains development work for **StringTie3's long-read module** that was integrated into the main [StringTie repository](https://github.com/gpertea/stringtie).

## Publication

> Shinder I, Pertea G, Hu R, Rudnick Z, Pertea M. **StringTie3 Improves Total RNA-seq Assembly by Resolving Nascent and Mature Transcripts.** *Nature Methods* (2025). [Accepted]
>
> Preprint: [bioRxiv 10.1101/2025.05.21.655404](https://www.biorxiv.org/content/10.1101/2025.05.21.655404)

---

## Attribution Note

These contributions were developed in this fork (commits [`23af2ba`](https://github.com/ishinder/stringtie/commit/23af2ba), [`e71e3c6`](https://github.com/ishinder/stringtie/commit/e71e3c6)) and integrated
into the main repository by Geo Pertea in [commit 29cd748](https://github.com/gpertea/stringtie/commit/29cd7485362756c28c5f5ad05193f233310f6b66),
where he also fixed an edge-case bug he identified during testing. The commit message
credits this work as "Ida's patch."

The full development history is preserved in this fork, and the implementation is
described in the Methods section of the StringTie3 paper.

---

## Commit History

### Commits in This Fork (by Ida Shinder)

| Commit | Date | Description | Stats |
|--------|------|-------------|-------|
| [`23af2ba`](https://github.com/ishinder/stringtie/commit/23af2ba) | 2025-09-24 | Feature/tts poly a graph (#3) - Main polyA implementation | +1006, -157 |
| [`e71e3c6`](https://github.com/ishinder/stringtie/commit/e71e3c6) | 2025-09-26 | Poly verbose v2 (#4) - Statistics and logging | Additional refinements |

### Related Commit in Main Repository

| Commit | Author | Description |
|--------|--------|-------------|
| [`29cd748`](https://github.com/gpertea/stringtie/commit/29cd748) | Geo Pertea | "long read infinite loop fix for Ida's patch" |

---

## Lines of Code

Based on `git blame` analysis of **this fork**:

| File | Lines Contributed | Description |
|------|-------------------|-------------|
| `rlink.cpp` | ~944 | Core polyA detection, transcript boundary logic |
| `rlink.h` | ~30 | Data structure definitions |
| `stringtie.cpp` | ~25 | Global statistics variables |
| `bundle.h` | ~14 | Read alignment data structures |
| **Total** | **~1,013** | |

---

## Long-Read Module Refinements

The long-read module refinements described in the StringTie3 paper (Methods: *Long-read mode refinement*) were developed in this fork.

### Overview

**1. Poly(A/T) Tail Detection**
- Scan up to 20 bp inward from alignment ends
- Mark terminal aligned poly(A/T) when window contains ≥5 A/T bases and ≥80% A/T content
- Handle both aligned segments and soft-clipped unaligned tails

**2. Terminal Alignment Artifact Handling**
- Detect erroneous terminal exons from poly(A) tails splice-aligning to A/T-rich genomic loci
- Discard single-exon alignments with terminal aligned poly(A/T) (frequent artifacts)
- For multi-exon: shorten terminal exon to 3 bp to preserve splice junction

**3. Poly(A) Site Anchoring (CPAS)**
- Cluster unaligned tail coordinates into 5 bp bins
- Promote bins with ≥20 supporting reads to cleavage/poly(A) anchors
- Constrain transcript boundaries during graph traversal

---

## Functions Implemented

| Function | Location | Purpose |
|----------|----------|---------|
| `poly_window_meets()` | rlink.cpp:510-545 | Core A/T content checking |
| `check_aligned_polyT_start()` | rlink.cpp:550-565 | Aligned poly-T at read start |
| `check_aligned_polyA_end()` | rlink.cpp:568-582 | Aligned poly-A at read end |
| `check_unaligned_polyA_end()` | rlink.cpp:585-600 | Soft-clipped poly-A detection |
| `check_unaligned_polyT_start()` | rlink.cpp:603-617 | Soft-clipped poly-T detection |
| `check_last_exon_polyA()` | rlink.cpp:619-637 | Last exon A-fraction check |
| `check_first_exon_polyT()` | rlink.cpp:639-656 | First exon T-fraction check |
| `shortenFirstExon()` | rlink.cpp:14118-14139 | Trim first exon to 3 bp |
| `shortenLastExon()` | rlink.cpp:14142-14163 | Trim last exon to 3 bp |
| `cluster_positions_with_counts()` | rlink.cpp:379-426 | CPAS position clustering |
| `add_cpas_trimpoint()` | rlink.cpp:428-447 | Merge CPAS into graph anchors |
| `ok_to_demote()` | rlink.cpp:362-370 | Conservative junction demotion |
| `junc_ratio_threshold()` | rlink.cpp:354-359 | Tiered demotion thresholds |
| `is_splice_between()` | rlink.cpp:450-453 | Check splice gap between nodes |
| `has_lr_witness_two_splices()` | rlink.cpp:483-506 | Long-read validation for splice pairs |

---

## Constants Defined
```cpp
POLY_TAIL_WIN         = 20    // Window size for aligned/unaligned end scan
POLY_TAIL_MIN_COUNT   = 5     // Minimum A/T bases in window
POLY_TAIL_MIN_FRAC    = 0.8   // Minimum A/T fraction (80%)
POLY_TAIL_MIN_CONSEC  = 5     // Minimum consecutive A/T bases
POLY_TAIL_STOP_COUNT  = 8     // Minimum reads for hardstart/hardend
CPAS_POS_BIN          = 5     // Cluster CPAS within 5 bp
CPAS_MIN_SUPPORT      = 20    // Require ≥20 reads per CPAS cluster
```

---

## Data Structures Added

**CTransfrag** (rlink.h):
```cpp
uint16_t poly_start_aligned;
uint16_t poly_end_aligned;
uint16_t poly_start_unaligned;
uint16_t poly_end_unaligned;
```

**CReadAln** (bundle.h):
```cpp
uint16_t aligned_polyT;
uint16_t aligned_polyA;
uint16_t unaligned_polyT;
uint16_t unaligned_polyA;
int sort_tiebreaker;  // Stable sorting in mixed mode
```

**CGraphnode** (rlink.h):
```cpp
bool hardstart:1;
bool hardend:1;
uint16_t polyStartAligned;
uint16_t polyStartUnaligned;
uint16_t polyEndAligned;
uint16_t polyEndUnaligned;
```

**Global Statistics** (stringtie.cpp):
```cpp
g_longread_alignments_total
g_longread_singleton_screen_total
g_longread_singleton_discard_total
g_longread_removed
g_longread_unaligned_tail_total
g_longread_shortened
```

---

## Development History

To explore the development history and implementation details, clone this fork:
```bash
git clone https://github.com/ishinder/stringtie.git
cd stringtie

# View commit history
git log --oneline --author="Ida"

# Count lines by author in key files
git blame rlink.cpp | grep -i "ida\|shinder" | wc -l

# View specific contributions
git blame -L354,430 rlink.cpp    # Constants and helper functions
git blame -L510,660 rlink.cpp    # PolyA detection functions
git blame -L14118,14165 rlink.cpp # shortenExon functions

# Show commit statistics
git show 23af2ba --stat
```

---

## Syncing This Fork

This fork's `master` branch may be behind the main repository. The contributions described above have been integrated upstream. For the latest StringTie3:
```bash
git clone https://github.com/gpertea/stringtie
cd stringtie
make release
```

---

## Contact

**Ida Shinder**
- Email: [ishinde1@jhmi.edu](mailto:ishinde1@jhmi.edu)
- Main StringTie repository: https://github.com/gpertea/stringtie
