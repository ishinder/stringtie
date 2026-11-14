# Read intake: `processRead()` and duplicate pooling

What happens to each alignment between the BAM reader and a queued bundle, and
how identical alignments are pooled into one `CReadAln`.

> **Verified against:** `origin/master` at `4bed517` (2026-11-14). Line numbers
> refer to that commit; trust function names if they drift.
> Status labels as in [LONG_READS_POLYA.md](LONG_READS_POLYA.md): **[code]**,
> **[tested]**, **[hypothesis]**.

## Where it runs [code]

`main()` reads alignments in coordinate order on a **single producer thread**.
For each alignment (`stringtie.cpp:1080-1166`) it:

1. Starts a new bundle when the read begins more than `runoffdist` past
   `currentend` (`stringtie.cpp:945`), or the chromosome changes. The finished
   bundle is queued to a worker.
2. Extends `currentend` with the read and with any guides it overlaps.
3. Calls `BundleData::evalReadAln()`, which returns `ovlpguide` (`stringtie.cpp:1156`).
4. Calls `processRead(currentstart, currentend, *bundle, hashread, alndata, ovlpguide)`
   (`stringtie.cpp:1166`). In `-e` mode only reads overlapping a guide get here.

Everything in `processRead()` is serial. Per-read cost there limits throughput
no matter what `-p` is.

## Steps inside `processRead()` (`rlink.cpp:710`) [code]

1. Long reads drop secondary and supplementary alignments.
2. Long-read polyA/T screening and terminal-exon trimming. This can return early
   or edit `brec.exons`; see [LONG_READS_POLYA.md](LONG_READS_POLYA.md) §1.
3. **Duplicate lookup** (not in merge mode); see below.
4. Merge mode instead sets `tinfo->qpass` from the coverage/FPKM/TPM filters.
   Rejection is deferred so repeated structures can rescue a filtered input
   ([MERGE_MODE.md](../MERGE_MODE.md)).
5. **New read:** a `CReadAln` is built. Its junctions are added to the bundle with
   `junction.AddIfNew()`, deletion-aware for long/mixed reads. It is appended to
   `readlist`.
   **Duplicate:** keep the smaller `nh` and add the polyA/T counts.
6. Add the read weight (`YC` tag, or the `YK` unitig coverage, divided by `nh`
   unless `nomulti` is set) to `read_count`. Update per-junction `nreads`, `nm`
   (mismatch-flagged support) and `mm` (long-anchor support).
7. Mate pairing through `hashread`, keyed by read name, mate start and `HI`.
8. Index the new read for later duplicate lookups (`rlink.cpp:1058-1065`).

## Duplicate pooling (issue #511) [code]

**What counts as "the same read":** `samereadmatch()` (`rlink.cpp:261`)
requires all of these to be equal:

- strand
- long/short flag
- unitig flag (with `-U`)
- exon coordinates (`exonmatch`)
- for long/mixed reads, the deletion offsets at each junction (`deljuncmatch`)

**How candidates are found:** stored reads are indexed by a 64-bit
`wyhash64` key over exactly those fields:

- `BundleData::last_read_with_key` maps a key to the newest `readlist` index with that key.
- `BundleData::prev_read_same_key[i]` links to the previous read with the same key (`bundle.h:453-458`).
- `processRead()` looks up `alnkey(brec,…)` and follows the chain newest first,
  confirming each candidate with `samereadmatch()`. The chain is longer than one
  step only after a hash collision.
- New reads are indexed at the **end** of `processRead()` with `readkey(rd)`,
  once their strand and unitig flag are final.

`alnkey()` and `readkey()` must stay in step with `samereadmatch()` and
`deljuncmatch()`: any pair the matcher accepts must get equal keys. The matcher
reads the two sides from different data (`brec.juncsdel` versus offsets derived
from `rd->juncs`), which is why there are two key functions.

**Debug check:** build with `-DREADIDX_CHECK`. Every lookup is then repeated as
a full backward scan, and `GError` is raised if the index and the scan disagree
(`rlink.cpp:816-824`).

**Before the index:** `processRead()` scanned `readlist` backwards over every
stored read with the same start. When very many distinct alignments share a
start (IGK/IGL loci in B-cell or plasma-cell long-read data), that is quadratic
and runs on the producer thread. History:

- PR #513 (`b022e91`, merged as `6a8ca72`) indexed by an exon-only hash. On a
  tonsil IGK locus the bundle went from 1,437 s to 17 s.
- PR #516 (`a856df5`, Geo) switched to keys over every compared field. With an
  exon-only key, reads that share exons but can never match pile into one chain:
  - a deletion right after an intron
  - with `-G`, a deletion at a guide junction
  - a strand taken from the polyA/T tail

### #511: what was measured

- **[tested] Reproduced on 3.0.4 with synthetic data.** All reads share one start
  and have distinct first-exon/intron lengths; run with `-L -c 2 -v -p 4`:

  | Reads | Time |
  | --- | ---: |
  | 200k, distinct starts (control) | 0.54 s |
  | 50k, same start | 2.68 s |
  | 100k, same start | 10.57 s |
  | 200k, same start | 47.8 s |

  In the 200k run, about 40 s pass before the bundle's "begins processing"
  line. macOS `sample` put 95% of main-thread samples in the scan loop. The
  symptom matches the report: one core at 100%, idle workers, samtools blocked
  in `pipe_write`.
- **[tested] The same-start scan is in every release from v1.1.0 to v3.0.3.**
  Check with `git show $t:rlink.cpp | grep -F 'readlist[n]->start==brec.start'`.
  Use `-F`, because an unescaped `[` silently matches nothing.
- **[code] The 2.2.2 → 2.2.3 changes are only in `count_good_junctions()`.**
  That runs after "begins processing", so it does not explain the reporter's
  longer pre-bundle gap on 2.2.3.
- **[tested] Hotspots are the immunoglobulin loci:** IGK (chr2:88.5-89.3 Mb) and
  IGL (chr22:21.9-23.1 Mb). Abundant Ig transcripts with varied V-J joins and
  ragged long-read ends give many distinct structures.
- **[tested] No command-line option avoids the scan in 3.0.4.**
- **[tested] An earlier "hybrid" fix gave byte-identical output.** It did a linear
  scan up to 64 entries, then indexed the same-start run.
  - `run_tests.sh` 11/11.
  - 28/28 identical GTFs.
  - A cross-check build agreed on every read, about 6.8M reads.
  - The 200k pile went from 46.5 s to 9.5 s; short-read timings were unchanged.

  It was dropped for the simpler bundle-wide index, which is allowed to change
  output (next section).

### Pooling behaviors to keep in mind

These were found with a Python model of the original scan. The model reproduced
StringTie's per-bundle distinct counts exactly. Each behavior was then checked
on a synthetic case.

- **[tested] Strand compares against the stored, possibly inferred, strand.**
  - A later *tagged* read can pool into an entry whose strand came from a polyA/T tail.
  - A later *untagged* copy cannot pool into it.
  - An untagged entry's strand is fixed by its first copy; later copies add tail
    counts but never re-infer strand.
- **[tested] The old scan missed duplicates across a trimmed read.** It only
  searched the contiguous same-start run at the end of `readlist`. A long read
  whose polyT first exon is removed is appended with a later start, which ends
  the run, so identical reads on either side of it stopped pooling. In a 3k
  synthetic case that gave 1,545 entries instead of 154. The bundle-wide index
  fixes this; it is one of the intended output changes.
- **[tested] A deletion right after an intron (`N D M`) never pools, even with an
  identical copy.** In `GSamRecord::setupCoordinates()` (gclib `GSam.cpp`) the
  `D` case does not reset the intron state, so an extra `juncsdel` entry is
  added, and `deljuncmatch()`'s count check always fails. The comment above
  `keystart()` (`rlink.cpp:270-276`) lists this case as "can never match"; it is
  hashed apart, not fixed. A separate issue was proposed but not filed.
- **[tested] Pooling is not a pure optimization.** Pooled entries share fields:
  - the minimum `nh` is kept;
  - polyA/T counts are summed. Shortening (LONG_READS_POLYA.md §2) needs zero
    unaligned evidence, so one soft-clipped copy blocks shortening for every copy.

  More pooling can therefore change output beyond merged counts.
- **[code] Order-dependent code to re-check when changing pooling:**
  - `hashread` mate pairing;
  - the mixed-mode `sort_tiebreaker` / `coordStableCmp`; a local commit (`1166acd`, not upstream) keeps arrival order for ties;
  - the `readlist.Sort()` for long reads in `infer_transcripts()`.

### Real-data test set for #511 [tested]

Public GM12878 nanopore BAMs are range-requestable, so `samtools view URL region`
can slice them:

- `https://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-cDNA-1D.pass.dedup.fastq.hg38.minimap2.sorted.bam`
- the DirectRNA equivalent

Caveats:

- **The tags are wrong.** They were aligned with `minimap2 -ax splice -uf`, which
  stamps `ts:A:+` on every read of this **unstranded** cDNA library. Realign
  without `-uf` to get junction-derived strand tags.
- **GM12878 is lambda-heavy.** IGL is a good stress test (about 200k alignments in
  the slice); IGK is not (about 5.5k). A kappa-heavy public dataset was not found.
- On 16 input/option combinations, unpatched and patched 3.0.4 gave identical
  output. Inferring strand from polyA/T did not break the index.
- CHESS 3.1.3 has 1,508 lines with doubled quotes in the chr2+chr22 subset;
  StringTie still loads every transcript.

## Open questions and hypotheses

- **The `N D M` deletion bug** (above): file it separately. Fixing it changes output.
- **Should strand and shortening use pooled evidence?** Should strand
  inference use all pooled copies rather than the first, and should "any
  unaligned evidence vetoes shortening" become evidence-weighted? No decision yet.
- **[hypothesis] Possibly long hash chains on real data.** Under the hybrid fix's
  exon-only key, chains reached 446 steps (mean 7.5) on realigned IGL. The
  suspected cause, same exons with different deletions next to junctions, is
  exactly what `a856df5`'s full key separates. Not re-measured on `a856df5`.
- **Benchmark not run.** A server benchmark on the StringTie 3 paper datasets
  was prepared but never run.
