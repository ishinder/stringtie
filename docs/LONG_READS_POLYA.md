# Long reads: polyA/T tail evidence

How StringTie uses polyA/T tails on long-read alignments (`-L`, and the long-read
side of `--mix`) to drop artifacts, assign strand, and anchor transcript ends.

> **Verified against:** `origin/master` at `4bed517` (2026-11-14). Line numbers
> refer to that commit. Function names are more stable than line numbers; if they
> disagree, trust the function and re-check this doc.
>
> **Status labels:** **[code]** read directly in the code at the cited location.
> **[tested]** confirmed by a run, with the command or data named.
> **[hypothesis]** reasoned but not tested. Hypotheses live only in
> [Open questions](#open-questions-and-hypotheses).

## Terms

| Term | Meaning in the code |
| --- | --- |
| aligned polyA end | The last aligned bases of the read (before the right soft clip) are A-rich: `check_aligned_polyA_end()` |
| aligned polyT start | The first aligned bases (after the left soft clip) are T-rich: `check_aligned_polyT_start()` |
| unaligned polyA end | The right **soft clip** starts with an A-run: `check_unaligned_polyA_end()` |
| unaligned polyT start | The left **soft clip** ends with a T-run: `check_unaligned_polyT_start()` |
| CPAS | Cleavage/polyadenylation site (transcript 3' end) |
| hardstart / hardend | `CGraphnode` flags meaning "a transcript may really start/end here"; they block end trimming and change transfrag compatibility |

Reads are in genome orientation. A plus-strand transcript's tail shows as A's at
the right end; a minus-strand transcript's tail shows as T's at the left end.

**Interpretation built into the code [code]:** an *unaligned* (soft-clipped) tail
is treated as a real polyA tail. An *aligned* A/T run with no soft-clipped tail is
treated as a genome-encoded A-run, i.e. a likely internal-priming artifact.

## Detection thresholds [code]

`rlink.cpp:400-406`

| Constant | Value | Use |
| --- | --- | --- |
| `POLY_TAIL_WIN` | 20 | bases scanned from the alignment boundary |
| `POLY_TAIL_MIN_CONSEC` | 5 | a run of 5 A/T starting at the boundary is enough on its own |
| `POLY_TAIL_MIN_COUNT` / `POLY_TAIL_MIN_FRAC` | 5 / 0.8 | otherwise: a prefix from the boundary with ≥5 matches and ≥80% A/T |
| `POLY_TAIL_STOP_COUNT` | 8 | unaligned-tail reads needed to mark hardstart/hardend |
| `CPAS_POS_BIN` | 5 | max gap (bp) when clustering tail positions |
| `CPAS_MIN_SUPPORT` | 20 | min weighted reads for a CPAS cluster |

The inline comments on `CPAS_POS_BIN` ("within 25 bp") and `CPAS_MIN_SUPPORT`
(">= 3 reads") do not match the values. Trust the values.

Bases are read straight from the BAM record with `readbase()` (`67b2508`), so only
upper-case `A`/`T` occur; the lower-case comparisons are harmless leftovers.

## Pipeline

### 1. `processRead()`: producer thread, once per alignment (`rlink.cpp:710`)

Applies to long reads only (`longreads`, or `brec.uval` in mixed mode). Secondary
and supplementary long-read alignments are ignored up front.

1. **Short-alignment screen** (`rlink.cpp:744-751`) [code]. An alignment with
   ≤2 exons that does not overlap a guide is **discarded** if it has an *aligned*
   polyT start or an *aligned* polyA end. The `&& !check_unaligned_...` clauses are
   commented out, so the read is dropped even when it also has a soft-clipped tail.
2. **Terminal artifact-exon trimming** (`rlink.cpp:754-798`) [code]. For ≥2 exons
   with no guide overlap: the last exon is removed if ≥80% of its read bases are A
   (`check_last_exon_polyA`), the first if ≥80% are T (`check_first_exon_polyT`).
   If both fire on a 2-exon read, the read is dropped. Otherwise the alignment
   strand decides: `strand==1` keeps the first exon, `strand==-1` keeps the last,
   `strand==0` removes both. The matching `juncsdel`/`alndata.juncs` entries are
   removed too.
3. **Evidence flags on a new read entry** (`rlink.cpp:865-879`) [code]. The four
   flags are set to 0/1. A trimmed exon is converted to *unaligned* evidence
   (`rmFirst` sets `unaligned_polyT=1`, `aligned_polyT=0`; `rmLast` likewise for A).
4. **Strand from tails** (`rlink.cpp:881-887`) [code]. Only when the alignment
   strand is 0: plus if any polyA evidence (aligned **or** unaligned) and no polyT;
   minus if the reverse; still 0 if both or neither.
5. **Duplicate alignments** (`rlink.cpp:928-941`) [code]. When the alignment
   matches a stored read, its four flags are added to the stored counts, which
   saturate at 65535. The stored strand is not re-evaluated.

### 2. `infer_transcripts()`: worker thread, before graph building (`rlink.cpp:17912-17935`)

**Internal-priming shortening** [code]. A read with an aligned tail but no
unaligned tail has that terminal exon shrunk to 3 bp (`shortenLastExon` /
`shortenFirstExon`, `rlink.cpp:14450-14495`). The tail end is chosen by strand:
+ → last exon, − → first exon, 0 → whichever ends have the evidence.

### 3. `build_graphs()`: CPAS anchors (`rlink.cpp:14553-14607`) [code]

Runs for `longreads || mixedMode`. Each long read adds its end position to
`raw_plus`, once per unaligned-polyA count, and its start position to
`raw_minus`, once per unaligned-polyT count. Strand-1 reads contribute only to
plus, strand −1 only to minus, strand-0 reads to both. `cluster_positions_with_counts()`
does single-linkage clustering (gap ≤5 bp). Each cluster with total weight ≥20 becomes a
`CTrimPoint` in `tstartend`: a sink link for plus, a source link for minus. It is
merged into an existing point within 5 bp. Only **unaligned** evidence feeds CPAS.

### 4. Transfrags and graph nodes [code]

- Read counts are summed into `CTransfrag::poly_*` fields as transfrags are built
  (`rlink.cpp:5010-5017`).
- `process_transfrags()` (`rlink.cpp:5923`) sums them onto the first/last graph
  node. On the minus graph (`s==0`) it sets `hardstart` when
  `polyStartUnaligned ≥ 8`; on the plus graph (`s==1`) it sets `hardend` when
  `polyEndUnaligned ≥ 8`.
- Effects of hardstart/hardend: the end-trimming in the coverage-drop logic is
  skipped (`rlink.cpp:4852-4940`), and transfrag compatibility checks treat the
  ends as fixed (`rlink.cpp:5760-5850`).
- `parse_trflong()` (`rlink.cpp:10411`; trimming at 10480-10530): when a transfrag
  has ≥8 unaligned tail reads, an extended path is cut back to that transfrag's
  last node (+) or first node (−).

### 5. Reporting

With `-L`/`--mix`, each bundle's log ends with a cumulative "Poly-tail stats"
line (`stringtie.cpp:1961`).

## Why strand from polyA/T can be wrong

Code facts that bear on this [code]:

- Strand inference (step 1.4) counts **aligned** A/T runs as well as soft-clipped
  ones. An A-rich genomic end (the internal-priming signal of step 2) can
  therefore set a strand. Step 2 only shortens the exon; it does not undo the
  strand.
- The tail-inferred strand reaches the stored read, never its junctions.
  `processRead()` creates them with the alignment's strand (`jstrand=strand`, still 0).
  `count_good_junctions()` fixes unstranded long reads from their junctions or coverage,
  but only for reads with `rd.strand==0` (`rlink.cpp:17748`). A tail-stranded read is
  therefore skipped by that rescue, even when its junctions would give a different
  strand. (`rlink.cpp:17699` copies `rd.strand` onto a junction only in the
  deletion-restore path.)
- A read whose strand came from its tail is stored with strand ±1, but the next
  identical alignment arrives with strand 0. `samereadmatch()` requires equal
  strands, so the two are **never pooled**. The comment above `keystart()`
  (`rlink.cpp:270-276`) acknowledges this. Each such alignment becomes its own
  `readlist` entry.

- After that, `build_graphs()` resets the read to strand 0 when its junctions are
  not stranded (`if(!keepstrand) rd.strand=0;`). A spliced read whose strand
  came only from its tail can therefore end up contributing no stranded
  transcript. A junction still unstranded at that point is treated as bad and
  the read is cut there. Geo's `polyA20.sam` reproducer (#515) stores 20
  entries, not 1, and produces no transcript.

### Measured: how often the tail contradicts the aligner [tested]

Tagged, spliced reads; "wrong" means the tail-implied strand disagrees with
the `ts` tag. Where CHESS introns can decide on chr19, they side with the tag
17,689 times and with the tail once.

| Data | Any tail | Soft-clipped tail only |
| --- | ---: | ---: |
| ONT cDNA, GM12878 chr19 (realigned) | 7.4% | 3.5% |
| ONT cDNA, neuron Day 0 | 2.4% | 1.2% |
| ONT direct RNA, WTC11 | 1.0% | 0.5% |
| PacBio, WTC11 | 0.02% | ~0% |

Narrowing the rule on GM12878 chr19:

| Condition | Wrong |
| --- | ---: |
| Soft-clipped tail at one end | 3.5% |
| … and no aligned A/T run at the opposite end | 3.1% |
| … and no aligned A/T run at either end | **1.6%** (best simple rule) |

A longer or stricter soft-clip requirement never went below about 1.7% and threw
away most correct tails. The current test is loose: 5 A/T in the first 6 bases
passes.

### Why ONT cDNA tails mislead

- **[tested] Where the wrong tails are.** Wrong tails are mostly a T-rich
  stretch at the *left* end of a plus-strand read. 32% of wrong soft-clipped tails
  lie within 50 bp of an annotated transcript start, against 0.02% of correct
  ones. RPS11 alone accounts for 15%; UBA52, RPL18A and RPS16 also cluster. In
  the largest cluster, 6,450 of 8,192 reads have the T-run aligned, not clipped.
- **[hypothesis, literature-backed] Second-strand mis-priming.** Leftover oligo-dT
  anneals to the internal A-run that a transcript T-run leaves in first-strand
  cDNA, so second-strand synthesis starts mid-transcript. The read is
  5'-truncated and begins at a genomic T-run, which the code reads as a
  minus-strand polyA end. See Sessegolo et al. 2019, *Sci Rep*,
  doi:10.1038/s41598-019-51470-9: the effect is detectable from about 9 T's,
  absent in direct RNA, and proposed rather than proven there. T-rich
  ribosomal-protein starts may contribute without any mis-priming; which cause
  dominates is unresolved.
- **[code, protocol text] PCR-cDNA is worse than direct cDNA.** The PCR-cDNA
  kit (SQK-PCB109) carries leftover oligo-dT into PCR with no clean-up; the
  direct cDNA kit (SQK-DCS109) cleans up before a one-cycle second strand. That
  fits 7.4% vs 2.4%, but the GM12878 library prep was not confirmed.
- **The opposite error: first-strand internal priming.** An A-run inside the
  transcript gives the right strand but the wrong 3' end. Steps 1.1 and 2 exist
  for this case.

### Strand-rule experiment (#515) [tested]

GM12878 ONT cDNA, all of chr19, realigned with `minimap2 -ax splice` (about
1M reads). The `ts` tag was removed from a random 50% of tagged reads (3 draws);
the run used `-L`. Untagged reads were scored against the removed tag, and
recovery counted the 3,303 intron chains found with full tags. The experiment
used code from before #513.

| Rule | Wrong strand | No strand | Chains recovered |
| --- | ---: | ---: | ---: |
| A: current code | 1,432 | 92,647 | 1,977 |
| G: junctions take the tail strand (Geo's `lr-polya-strand-junctions`) | 5,814 | 5,715 | 2,448 |
| D: spliced reads never take tail strand | 17 | 5,916 | 2,236 |
| E: tail strand used last | 19 | 5,901 | 2,241 |
| E + per-read junction counts | 18 | 5,692 | **2,743** |
| F: per-junction vote | 32 | 5,714 | 2,728 |

- The vote's gain came from per-read counts (see the over-count note in
  [STRAND_AND_JUNCTIONS.md](STRAND_AND_JUNCTIONS.md)), not from voting.
- On well-tagged real data (neuron ONT cDNA, 99.6% of spliced reads tagged), every
  rule gives the same sensitivity and precision. The rule only matters when tags
  are missing.
- **Proposed direction, not implemented:**
  - strand order: aligner tag → stranded junctions from other reads → coverage →
    tail last, using only a soft-clipped tail with no aligned A/T run at either end;
  - per-read junction counts;
  - no vote.

### Discarding short polyA alignments (step 1.1): design history [tested]

- Discarding regardless of a soft-clipped tail is deliberate. The
  `!check_unaligned_…` terms were commented out on purpose, to suppress fragment
  transcripts with different starts. The screen was later extended from 1 to 2 exons.
- **Re-tested.** Restoring the commented-out `!check_unaligned_…` terms on neuron
  Day 0 (whole library):
  - transcripts: +1,089 (+76 matching CHESS);
  - transcript precision: 48.1 → 47.2%;
  - locus precision: 73.4 → 71.3%;
  - sensitivity: unchanged (16.27 → 16.32%).

  On IGL, most screened reads with clipped tails are 3' fragments of V-J-C
  molecules whose V region didn't align. This supports keeping the terms out.
- With `-G`, reads overlapping a guide skip the screen and the trimming
  (`!ovlpguide`), so `-G` runs keep more reads. Neuron Day 5 rep 1: 2,907,199 vs
  3,258,748 reads kept.

## Open questions and hypotheses

Untested unless marked otherwise. Each line says what would settle it.

- **[hypothesis] Stats undercount shortening.** `g_longread_shortened` is
  incremented only when a *first* exon is shortened (`rlink.cpp:17920`,
  `17933`). It is also incremented on worker threads without a lock, so with
  `-p >1` the count may be racy. This affects the log line only.
- **[hypothesis] Some transfrags may be skipped by node aggregation.**
  `process_transfrags()` skips transfrags whose first node is the source **or**
  whose last node is the sink (`tf->nodes[0]==0 || tf->nodes.Last()==gno-1`), but
  the comment says "only source/sink". If long-read transfrags carry source/sink
  nodes at this point, their tail evidence never reaches the graph nodes.
  *Settle by:* logging how many long-read transfrags are skipped on a test bundle.
  The "Reset node-level counters first" comment has no matching code; the counters
  start at 0 in the `CGraphnode` constructor (`rlink.h:437`), so this is probably
  a stale comment only.
- **[hypothesis] Unpooled tail-stranded reads.** Unpooled duplicates (see above)
  should not change the CPAS or transfrag totals, because those sums run over all
  entries. They do increase `readlist` size and per-bundle work for deep
  polyA-rich loci. *Settle by:* comparing `readlist.Count()` with and without
  strand inference on a deep ONT locus.

- **[tested on IGL; impact unmeasured] Smaller defects:**
  - Evidence counters are attached to every piece of a split read, including
    pieces without the read end: about 30k of 77k uses on IGL.
  - The exon A/T fraction (`check_first_exon_polyT` / `check_last_exon_polyA`)
    assumes no indels inside the exon; 38 decisions change on IGL.
  - After whole-exon removal, the forced "unaligned" evidence sits on a splice
    site, not a read end.
  - One soft-clipped copy blocks shortening for a whole pooled entry. On IGL
    only 466 of 7,331 blocked entries were blocked by a *different* copy.
- **[code] Dead or misleading bits with no output effect:**
  - Aligned-evidence counters on nodes and transfrags are never read.
  - The `CReadAln` copy constructor does not copy the polyA counters; no caller
    of it was found.
  - The 2-exon read dropped when both ends are A/T-rich is not counted in the stats.
- **Not yet tested:**
  - a full ablation with the screen off;
  - a gffcompare class breakdown of the 1,089 extra transcripts;
  - a strand-aware screen;
  - keeping a clipped-tail read when the genome beyond the aligned run is not A/T-rich;
  - accuracy of first/last-exon trimming;
  - `--mix` with trimmed long reads.
- **IGL is a misleading test locus for strand questions.** At IGLC2/IGLC3 the
  minimap2 tags are wrong; reads use an acceptor 7 bp upstream with a CT-AC
  motif. An early "tail beats tag" conclusion drawn from IGL was reversed on
  chr19. Use chr19 or whole libraries.
