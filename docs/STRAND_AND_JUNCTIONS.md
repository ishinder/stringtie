# Long reads: strand, junctions, and deletions next to introns

How a read and its junctions get a strand, how deletions beside an intron are
handled (#514), and the order and coverage effects found while studying them.
Mostly long-read (`-L`) and mixed (`--mix`) code paths.

> **Verified against:** `origin/master` at `4bed517` (2026-11-14). Line numbers
> refer to that commit; trust function names if they drift.
> Status labels as in [LONG_READS_POLYA.md](LONG_READS_POLYA.md): **[code]**,
> **[tested]**, **[hypothesis]**. Hypotheses are collected in
> [Open questions](#open-questions-and-hypotheses).
>
> **Not in 4bed517:** Geo's branches `origin/lr-junction-deletions` (#514) and
> `origin/lr-polya-strand-junctions` (#515), and the local stable-sort commit
> `1166acd`. Statements about them are marked as such.
>
> **Experimental numbers:** "GM12878 chr19, half the tags removed" results come from
> instrumented builds based on pre-#513 master unless marked "Geo's branch".
> They compare rules against each other; re-measure before quoting them as
> absolute effects on current master.

## Terms

| Term | Meaning in the code |
| --- | --- |
| tentative junction | A `CJunction` with `strand==0` whose coordinates were widened over an adjacent deletion in `processRead()`. Told apart from a truly unstranded junction only because its ends do not line up with the read's exon ends. |
| `juncsdel` | `GSamRecord::juncsdel` (`gclib/GSam.h:45`): meant to hold one `GSeg` per intron, `.start` = bases deleted just before the intron, `.end` = just after. Exists only on the incoming record, not on `CReadAln`. |
| DEL AWARE | The code that widens junctions over deletions (`processRead()`) and later keeps or restores them (`count_good_junctions()`). Long/mixed reads only. Added in May 2021 (`6567202`, `28a9a2e`, `122f896`). |
| guide snap | With `-G`, moving a read's splice site onto a guide splice site within `-E` bp. |
| neighbour step | In `build_graphs()`, demoting a weak junction in favour of a better one within `-E` bp. Runs with or without `-G`. |

## 1. How a read gets its strand [code]

On the producer thread, before `processRead()`:

1. **`XS` tag**, else the minimap2 **`ts`** tag. `ts` is read-relative, so it is
   flipped for reverse-strand alignments: `CAlnTags::spliceStrand()`
   (`bundle.h:276-282`), called at `stringtie.cpp:861`. Since `d1ac335` all tags
   come from one pass (`alntags.get()`); `GSamRecord::spliceStrand()` is no
   longer on this path. There is no splice-motif check.
2. **Library type** (`--fr`/`--rf`), only if no tag (`stringtie.cpp:862-877`).
3. **Guide overlap** (`-G` only), only if still unknown:
   `BundleData::evalReadAln()` sets the strand when all similar-size guide
   exon overlaps agree (`bundle.cpp:235-237`).
4. **PolyA/T tail** (long reads, strand still 0): see
   [LONG_READS_POLYA.md](LONG_READS_POLYA.md), Pipeline step 1.4.

Junctions are created in `processRead()` with the read's arrival strand
(`jstrand=strand`, `rlink.cpp:900`), so a tail-inferred strand never reaches
them; see LONG_READS_POLYA.md, "Why strand from polyA/T can be wrong".

Later, on the worker:

- `count_good_junctions()`: the guide snap can change a read's strand
  (`rlink.cpp:17648-17649`, `17668-17669`), and the unstranded pass (§2) gives
  strand to spliced reads still at 0.
- `build_graphs()`: a junction still at strand 0 is a bad junction. The read is
  cut there, and if none of its junctions is stranded the read loses its strand
  (`rlink.cpp:15390-15395`). The pieces count only as coverage.

### Alignment pitfalls [tested]

- **minimap2 `-uf` stamps `ts:A:+` on every alignment**, spliced or not. On an
  unstranded library, strand then equals read orientation, which is wrong for
  about half the reads, and polyA/T strand inference never fires. The public
  GM12878 nanopore cDNA BAMs (minimap2 2.5 `-ax splice -uf -k14`) are like this
  (IGL: all 143,660 primary alignments `ts:A:+`; 58,646 forward, 85,014
  reverse). **Realign without `-uf` before testing anything about strand.**
- **pbmm2 (`--preset ISOSEQ`, v1.18.0) writes no `ts` or `XS`.** Every PacBio
  read aligned this way reaches StringTie at strand 0 and depends entirely on
  the unstranded pass.
- **IGL is a poor locus for strand questions.** At IGLC2/IGLC3 the minimap2 tags
  are wrong (reads use an acceptor 7 bp upstream, motif CT-AC), almost all reads
  are on +, and V-J joins have no splice motif. Use chr19 or whole libraries.

## 2. The unstranded pass in `count_good_junctions()`

### Mechanism [code]

`rlink.cpp:17797-17864`. Reads are collected at `17748` (spliced long/mixed
reads with `rd.strand==0`) and handled one at a time, in `readlist` order:

1. If any of the read's junctions is stranded, or a stranded junction exists at
   the same coordinates, take that strand (`17810-17834`).
2. Otherwise compare + and − coverage over the read's exons (`17835-17840`).
   This is done only for `longreads`; in `--mix`, unstranded long reads get a
   strand from junctions only (`pos`/`neg` stay 0).
3. For each of the read's strand-0 junctions: if a stranded copy exists, move
   the read to it and add counters (`17848-17859`). Otherwise **relabel the
   shared record in place**: `else rd.juncs[i-1]->strand = rd.strand;` (`17860`).

### The first read decides a shared junction [code] [tested]

Because of step 3's in-place relabel, the first read that reaches a shared
strand-0 junction sets its strand for every later read, which then find it as a
stranded junction in step 1. There is no vote.

- [tested] On 9 tonsil IGK reads, two tie orders gave the transcript on + in one
  run and on − in the other (confirmed in a debugger).
- [tested] GM12878 chr19, half the tags removed: 479 unstranded junctions had
  reads split across strands and 11,874 had reads left unstranded; with the
  strand order proposed for #515 (LONG_READS_POLYA.md, "Strand-rule experiment"),
  32 and 0. Proposed principle: a junction should never exist on both strands.

A short comment about this line was posted on #515.

### Unstranded-pass over-count [code] [tested]

When a read moves to an existing stranded junction (`17848-17857`), six counters
of the strand-0 record are added as **whole totals, once per read**: `nreads`,
`nreads_good`, `nm`, `mm`, `leftsupport`, `rightsupport`. Every read that shares
the record adds the full totals again, and the old record is not reduced.

- 6 stranded + 4 unstranded reads on one junction: `nreads` = 22, not 10.
- 30 tagged + 30 untagged: 930, not 60.
- Inflated support makes `good_junc()` drop weaker junctions that share the
  splice site, through `jd.nreads*10<ERROR_PERC*jd.leftsupport` (and
  `rightsupport`) at `rlink.cpp:14394`. In one chr19 draw, 1,148 junction
  records were lost this way.

**Fix tested, not implemented.** Move each read's own share: `read_count` for
`nreads`, and the fraction `read_count/nreads` of the old record for the other
five, subtracting the same amounts from the old record.

- GM12878 chr19, half the tags removed (3 draws): intron chains recovered rose
  from 2,241 to **2,743 of 3,303**, the largest single effect in the #515
  experiment. A "move once, then zero the old record" variant did as well
  (2,741 vs 2,738); the fractional form was preferred.
- It must go in with the self-confirmation fix (§3, defect 6): that bug inflates
  `nreads` and drives the fraction toward 0 (5 junctions lost in one run).

### Unstranded junctions cut reads [code] [tested]

A spliced read with no strand evidence keeps strand-0 junctions, which
`build_graphs()` treats as bad. Ten identical untagged one-intron reads give two
single-exon transcripts with strand `.`, coverage 10 each. Design direction
(not implemented): do not read GT/AG from the genome; report such reads as a
spliced `.` transcript only if they are not used in any stranded transcript.

When the same junction exists on + and −, `build_graphs()` removes the strand of
the copy with fewer `nreads` **and** fewer `nreads_good` (`rlink.cpp:14779-14786`).
Counts inflated by the over-count above feed this comparison.

## 3. Deletions next to introns (#514)

### Pipeline [code]

**A. CIGAR parsing.** `GSamRecord::setupCoordinates()` (`gclib/GSam.cpp:214-262`)
fills `juncsdel`. Only StringTie's bundled gclib has it; a fix touches no other
project.

**B. `processRead()`.** For a **new** entry (long/mixed reads), an intron with a
non-zero `juncsdel` is stored as a widened junction with strand 0, unless the
aligned intron is guide-matched (`rlink.cpp:906-912`). The read's exons stay as
aligned. Duplicate pooling compares `juncsdel` with offsets recomputed from the
stored entry (`deljuncmatch()`, `rlink.cpp:239-249`); see
[PROCESSREAD.md](PROCESSREAD.md), "Duplicate pooling".

**C. `count_good_junctions()`, `/* DEL AWARE*/`** (`rlink.cpp:17673-17705`). For
a strand-0 junction whose ends differ from the read's exon ends, look up
`CJunction(start,end,rd.strand)` at the widened coordinates (`17677`). If found,
move the read there and trim its exon ends. Otherwise restore the as-aligned
junction, reusing an existing record or rewriting the tentative one.

Support rules, verified with synthetic reads: one clean read with the wider
intron is enough; reads with the same deletion do not support each other; the
wider intron wins even when outnumbered (`122f896` set this order deliberately).

### Defects present in 4bed517

1. **Two records for one intron when the deletion follows it (`N D M`)** [code]
   [tested]. The `D` case (`GSam.cpp:235-242`) adds a record but does not clear
   `intron`, so the next `M` adds another (`227-231`). `100M1000N3D97M` gets
   `(0,3),(0,0)`. `deljuncmatch()` fails its count check, so identical copies
   never pool; in multi-intron reads later introns read a shifted record.
   Since `a856df5` the index keys these reads apart, so they no longer slow the
   lookup, but they are still stored once per copy (PROCESSREAD.md, "Pooling
   behaviors"). Scale: PacBio PBMC IGK 50,184 of 80,500 entries never pooled;
   HG002 ONT IGL 64% of reads. Most of these are 1-bp deletions at V-J
   rearrangement gaps, not real splices (PBMC PacBio: 22.5% of V-J gaps vs
   26 of 304,994 real introns; GM12878 ONT IGL: 10.2% vs 1.3%). Also, two `D`
   ops in a row keep only the last length (`del=`, `GSam.cpp:236`).
2. **Wrong side for a deletion before the intron (`D N M`)** [code] [tested].
   `GSeg(uint s,uint e)` swaps its arguments when `s>e` (`gclib/GBase.h:381-384`).
   `GSeg(prevdel,0)` at `GSam.cpp:228` therefore stores `97M3D1000N100M` as
   `(0,3)`, widening on the acceptor side. GTF output shows `1-100,1104-1200`
   instead of `1-100,1101-1200`. In effect, "widen if supported" has only ever
   worked for deletions after an intron. **REVISED:** an early session note said
   `D N` was handled correctly; it gets one record, on the wrong side.
3. **Guide-matched intron never pools** [code] [tested]. With `-G`, the widening
   is skipped when `alndata.juncs[i-1]->guide_match` (`rlink.cpp:907`), so the
   stored offsets recompute as 0 while the incoming `juncsdel` is not 0.
   Identical reads never pool, and pooling depends on arrival order. The guide
   is matched against the *aligned* intron; if the annotation has the wider
   one, the read is widened as usual.
4. **DEL AWARE adds whole totals** [code] [tested]. Moving to an existing record
   adds `rd.juncs[i-1]->nreads` (and `mm`, `nm`) at `17681-17683` and
   `17693-17695`, the record's total, not the read's count. Geo's `count.sam`
   gives 12 where 6 is correct; his branch fixes it (`85b71ae`).
5. **`-G` bypasses the support check** [code] [tested]. DEL AWARE is the
   `else` of `if(modified)` (`17632`). In any `-G` bundle where some junction was
   snapped to a guide (`modified=true`, `17524`), the `modified` branch moves
   exon ends onto a widened strand-0 junction without support, and the read
   loses its strand (`17668-17669`). 20 reads `97M3D1000N100M300N100M` give one
   transcript without `-G`, none with a guide that corrects only the second
   junction. Present since `122f896` (2021). Not on Geo's list.
6. **An unstranded read confirms its own tentative junction** [code] [tested].
   When `rd.strand==0`, the lookup key at `17677` is the read's own tentative
   record (same start, end, strand 0). The widened junction is always accepted,
   even with 10 clean reads on the as-aligned intron, and `17681` adds the
   record to itself, doubling `nreads` per entry (6 entries: 6 → 384 on master;
   6 → 12 on Geo's branch; 6 is correct). Output GTF can look normal; the
   damage is in internal counts. Chr19, half the tags removed (Geo's branch):
   41,114 self-accepts, where the as-aligned intron was annotated 31,609 times
   and the widened one 20 times. Not on Geo's list.
   Reproducer with Geo's `mkreads2.py` from #515 (`name:count:CIGAR:strand`,
   `.` = no tag):
   ```sh
   python3 mkreads2.py selfjunc.sam d1:1:100M1000N3D97M:. d2:1:100M1000N3D96M:. \
     d3:1:100M1000N3D95M:. d4:1:100M1000N3D94M:. d5:1:100M1000N3D93M:. \
     d6:1:100M1000N3D92M:. c:10:100M1000N100M:.
   stringtie -L -m 50 -o out.gtf selfjunc.sam
   ```
7. **The restore path rewrites a shared record in place** [code] [tested]. When
   no as-aligned record exists, `17698-17703` changes the tentative record's
   start, end and strand. Other reads may share that record, and `junction` is
   not re-sorted until `17795`, while `junction.Found()` binary searches
   continue. On chr19 (half tags removed), about 2,000 lookups missed and about
   800 keys were duplicated.

### Geo's branch `lr-junction-deletions` (#514, not merged)

Geo lists four cases: wrong side (2 above), two records (1), `-G` never pools
(3), rescue over-count (4). His fixes 1-3 write one record per intron and, at
guide-matched junctions, only require stored offsets of 0 inside
`deljuncmatch()`. He reports NA12878 chr1 precision 40.5 → 41.5, sensitivity
unchanged.

- [tested] Tip `48a5e68` agrees with an independent prototype on 46 of 51
  synthetic CIGAR cases and all 16 pool/no-pool pairs; it changes `mix_reads`
  and `mix_reads_guided` in `run_tests.sh` (coverage only, up to 3.4%).
- [tested] Differences: `100M1000N5D1000N100M` (N D N) gets one record for two
  introns, and `processRead()` reads `juncsdel[1]` out of range (debug build
  aborts; release build widens the wrong intron). `3D2D` keeps only 2. Neither
  shape occurs in about 4.4M real spliced alignments (minimap2 and pbmm2, IGK,
  IGL and chr19 sets). `N I N` does occur (about 1,000 in tonsil) and is
  handled correctly by all builds.
- It does not fix defects 5-7. Replies to #514 and #515 are drafted, not posted.

### Prototype fixes for 6 and 7 (V4, V4n) [tested]

Tested on top of `48a5e68`; not in any branch or PR.

- **V4:** add `bool aligned` to `CJunction` (false by default). Set it in
  `processRead()` when the junction is not widened, and in DEL AWARE's
  "new junction" case. In DEL AWARE, accept the read's own record (`nj==dj`)
  only if `dj->aligned`, and never add counts when `nj==dj`. Results: 6
  deletion reads + 10 clean as-aligned reads → as-aligned intron, count 16;
  6 + 10 clean wider reads → wider intron, count 16. Geo's reproducers and the
  repo tests are unchanged.
- **V4n:** in the restore path, `AddIfNew` a new record instead of rewriting the
  shared one.
- **REFUTED (V1):** "never accept the read's own record". The read's tentative
  record is the same object as a clean unstranded wider read's record, so this
  corrupts records that clean reads rely on.
- **REFUTED:** `intron=false` in the `D` case as a one-line fix for defect 1. It
  false-pools a D-before read with a D-after read (294/325/144 wrong merges on
  three real inputs), leaves defect 2, and breaks hanging-intron handling that
  also reads `intron` since `62551bb` (`GSam.cpp:251`, `287`).

### "One clean read" loses the junction [code] [tested]

1 × `97M1003N100M` plus 20 × `97M3D1000N100M`: Geo's branch outputs two
single-exon pieces. Master outputs the spliced transcript only by accident
(defect 2 sends the reads to the as-aligned intron); the mirror case fails on
master too.

- Cause: coverage is added (`17617`) before DEL AWARE moves the exon ends, so
  the deleted bases stay counted. The `higherr` coverage check in
  `build_graphs()` (`tolerance=1-ERROR_PERC`, `rlink.cpp:14931`; `rightcov >
  tolerance*leftcov` at `14953`) then sees 21 → 20 and marks the junction bad.
- The junction is lost whenever clean reads are 10% or fewer of the reads on it.
  Fixing counters alone does not help.
- A 6-line fix (subtract the read's coverage with `cov_edge_add(..., -read_count)`
  before moving exon ends, long reads only) fixes all synthetic cases and passes
  the regression sets, but changed **0** structures on chr19, IGL and PacBio
  IGK; the pattern appeared once in about 1.3M alignments. Not adopted.
  Do not "move" coverage instead; see §5.

## 4. Read order and ties

### The `-L` sort has no tie-break [code] [tested]

`infer_transcripts()` sorts long reads with a plain `bundle->readlist.Sort()`
(`rlink.cpp:17940-17941`): `GList` quicksort on (start, end). The stable path
with `sort_tiebreaker`, `coordStableCmp` (`rlink.cpp:67`) and `pair_idx` remap
(`17943-17985`) runs for `mixedMode` only. The sort was added with the polyA/T
feature (Sept 2025), because trimming and shortening break start order; it is
still needed with the hash index.

- One entry fewer moved 13,965 of 80,798 tied positions on IGL. This explains
  Geo's "one collapsed duplicate changes 36 of 70 transcripts".
- Paired short reads given `-L` keep stale `pair_idx`; don't run them with `-L`.
- Local commit `1166acd` (not upstream, not in a PR) routes `-L` through the
  mixed block; accuracy-neutral on neuron ONT Day 0/3/5, for reproducibility.

### Tie order changes output even with a stable sort [tested]

A stable sort only removes dependence on *unrelated* entries; ties still follow
BAM order. Three sites depend on the order of tied reads:

- **A. End trimming against partial coverage.** `update_abundance()`
  (`rlink.cpp:4835`; trim tests `4860-4957`) compares node coverage that is
  still accumulating while reads are walked, so the first of two tied reads is
  trimmed and the second is not. Class: bug.
- **B. `longtrCmp` ties.** `transfrag.Sort(longtrCmp)` (`6112`; comparator
  `1196-1210`) returns 0 for equal abundance, node count and pattern size, and
  the greedy keep loop picks the representative by position (`6140-6320`).
  Class: arbitrary choice; changing `longtrCmp` was judged not worth it.
- **C. Shared junction records in `count_good_junctions()`.** The first-read
  relabel (`17860`, §2) and the DEL AWARE in-place rewrite (`17698-17703`).
  Class: bug.
- Plus float sums in read order (node coverage reaches 1.26e9 on tonsil).

With all of these neutralised, 8 tie orders gave byte-identical output on every
input. Without: tonsil `-L` gave 8 distinct outputs from 8 orders (18
transcripts in only some; coverage swings up to 96%); GM12878 IGL `-L` 8/8
distinct, up to 68%.

For equivalence checks use `-p 1` (STRG IDs and line order vary with `-p >1`),
and do not mix Linux and macOS: post-sort tie order differs between them.

## 5. Stale coverage [code] [tested]

`add_read_to_cov()` runs at the top of each read's iteration in
`count_good_junctions()` (`rlink.cpp:17617`), before the read is changed:

| Later change | What goes stale | Wrong output shown? |
| --- | --- | --- |
| DEL AWARE exon-end moves (`17685-17686`) | exon ends | yes (§3, one clean read) |
| `-G` guide-snap exon-end moves (`17632-17670`) | exon ends | yes |
| `-G` snap strand change | strand | stale state only |
| unstranded-pass strand (§2) | strand | stale state only |
| `build_graphs()` neighbour exon-end rewrite | exon ends | stale state only |
| `build_graphs()` strand resets | strand | likely low impact |

`-G` example: 20 identical one-intron reads with a guide donor 3 bp away give
two single-exon pieces with `-G`, one transcript without. On chr19 `-G`, 204
junctions were dropped only because of stale coverage, none annotated.

**REFUTED fix: move the coverage update after the read changes.** It fixes the
synthetic cases but hurts with `-G`: IGL precision 69.6 → 65.3, PacBio IGK
51.3 → 40.8, chr19 `-G` 15 structures removed and 19 added, and the repo test
`long_reads_guided` breaks. Stale coverage accidentally rejects bad guide snaps
across short exons the aligner missed. At BSG, about 1,284 reads snap 19 bp
across an annotated 25-bp exon; with the move, a false exon-skip isoform appears
at coverage 1,234.

## 6. The `-E` splice-site window [code] [tested]

One option drives two steps; `-E 0` disables both, and they cannot be set
separately.

- **Guide snap** (`-G`, long/mixed): `count_good_junctions()`
  `rlink.cpp:17440-17520`, reach |d| ≤ E. The read's exon end moves onto the
  guide site. Opposite-strand junctions are searched within `2*sserror+1`
  (`17533`, `17569`).
- **Neighbour step** (always, inside `if(higherr)`, `build_graphs()`
  `14876-15147`): window `juncsupport` = E for `-L`, `sserror/DROP` = 2E for
  `--mix` (`14879-14881`), `-a` for short reads. Branches for a less reliable
  neighbour are also gated by `< sserror` (`14988`, `15017`, `15082`, `15111`).
  Demotion requires `ok_to_demote()` (`rlink.cpp:409-425`): a junction with ≤25
  reads is demoted if the neighbour has at least as many; above that the ratio
  must be ≤0.5, ≤0.25 or ≤0.10 (>200 reads). [tested] For short reads, one
  intron differed between `-E 0` and `-E 25`.

**Neighbour reassignment is broken (bug since v2.1.6)** [code] [tested]. In the
`build_graphs()` read loop the exon ends are moved to the neighbour
(`15261-15262`), but the read is not moved onto the neighbour junction; its own
junction loses its strand (`15210-15211`) and the read is cut. The reassignment
code (`15263` on) ran 0 times in 368 traced runs. Commit `6567202` (2021-05-27,
first tag v2.1.6) commented out the earlier `jd.mm>=0` condition (`15218`) and
enabled `if(jd.mm>=0)`; for a stranded neighbour, `jk`/`ek` are also reset to
-1 (`15232`, `15253`), which skips the search. Restoring it moves about 46,500
chr19 reads at `-E 25`; every metric moves by ≤0.4 points.

**Error profile [tested].** ONT near-miss sites are motif-quantized (±1 bp is
nearly absent; modes at 4, 7 and 9 bp). SIRV: misplaced sites lie within 11 bp
(PacBio 8). Real sites 3 bp apart are common (44.5% of such acceptor pairs have
Illumina support), so distance and read ratio cannot separate them from errors.
**REFUTED:** a one-sided (shift-side) rule; real pairs lie mostly on the
"error" side (56 vs 31).

**`--junc-bed` matters.** The StringTie 3 paper's long-read BAMs were aligned
with minimap2 `--junc-bed`; 97.7-99.5% of read introns then match the list
exactly, and `-E 0` → 25 gains 1-3 precision points for 1,100-1,400 chains.
Without a junction list (LRGASP WTC11, unguided `-L`, correct chains / chain
precision %):

| `-E` | ONT dRNA | ONT cDNA | PacBio |
| --- | --- | --- | --- |
| 0 | 17,002 / 47.7 | 17,992 / 43.9 | 22,364 / 55.2 |
| 5 | 16,478 / 58.3 | 17,603 / 52.0 | 21,669 / 56.1 |
| 25 (default) | 16,033 / 61.9 | 17,125 / 54.8 | 20,996 / 56.8 |

Going 0 → 25 loses about 1,800 annotated introns per library, 95-97%
Illumina-supported. A constant 0.25 ratio in `ok_to_demote()` at `-E 25` beat
`-E 5` on every platform (dRNA chain Sn/Pr 12.5/58.8 vs default 12.0/61.9); it
needs a code change. **Guided:** every base of guide window costs held-out
nearby introns (dRNA: 498 of 1,927 recovered at 0, 193 at 25) and adds false
exon skips (22 → 164); with the full guide, 25 gives 200-650 more correct
chains than 0. **Outcome:** the default stays at `-E 25`.

## 7. Guided mode drops well-supported novel isoforms [code] [tested]

Two rules in `printResults()` interact:

1. A non-guide prediction overlapping a guide-matched one is switched off if
   its coverage is below `ERROR_PERC` (0.1, `rlink.h:14`) of the guide's
   (`rlink.cpp:22233-22266`).
2. The later duplicate merge of same-intron-chain predictions (loop at `22351`,
   which does not skip switched-off records) copies the higher-coverage
   prediction's ends and coverage into the other record but keeps that
   record's `flag`. The flag is restored only when the winner is a guide:
   `if(pred[n]->t_eq) pred[m]->flag=true;` (`22478`, `22508`).

So a novel isoform with coverage 42,554 merged into a switched-off copy with
coverage 210 disappears. On a chr19 locus (48,450 reads) it has coverage 47,045
without `-G` and is absent with a guide holding only the neighbour isoform.
Both rules came with `2bdb26f` (2021-10-25, first tag v2.2.0); no option
affects them. On chr19 the merge path fired 37 times (coverage up to 43,767).
Small synthetic cases do not trigger it. No fix written.

## Open questions and hypotheses

Each says what would settle it.

- **Should one clean read be enough to widen a junction?** It rarely applies (14
  of about 840k reads on chr19). A design question for the maintainers.
- **First-read relabel (§2):** fix by grouping or a vote, in a separate change.
  *Settle by:* the tie-order matrix (§4) plus the #515 scoring on chr19.
- **PR for sort + tie sites A and C:** not started; needs a maintainer decision.
- **[hypothesis] Inflated counters make the wrong copy win** in the both-strands
  step (`14779-14786`). *Settle by:* a synthetic case with the over-count.
- **[hypothesis] Guide snap skips acceptors.** On chr19 at `-E 25`, 3,740 acceptor
  ends were not snapped although a guide site was in the window. *Settle by:*
  tracing one such read through `17486-17514`.
- **Snap exon-skip check (in progress).** Coverage fix plus refusing a snap
  whose junction is unannotated and jumps a complete annotated exon ≤50 bp (and
  cutting the read there). Chr19 full guide `-E 25`: false skips 20 → 4,
  precision 86.1 → 86.6; BSG fixed. *Settle by:* the unfinished
  whole-transcriptome dRNA run.
- **Most guided losses of held-out introns are elsewhere** (109 of 129 nearby,
  71 of 89 random): the isoform is never assembled. *Settle by:* tracing one
  such locus through path selection.
- **Other `-E` changes not adopted:** separate guide and neighbour windows; the
  0.25 ratio rule; restoring neighbour reassignment.
- **[hypothesis] Insertions next to junctions** could be hints the way
  deletions are. A design question for the maintainers.
- **[hypothesis] V-J deletions are rearrangement leftovers** (no splice motif),
  not splice errors. Counts are measured; the explanation is not.
