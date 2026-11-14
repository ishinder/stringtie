# Guides, `-e`, `-A` and abundance filters

How reference transcripts (guides, `-G`) are tracked through a run, what the
`-A` gene table and the `-e` mode actually report, how `-M` filters
multi-mapped loci, and where opposite-strand junctions and transcripts are
dropped. Merge internals are in [MERGE_MODE.md](../MERGE_MODE.md).

> **Verified against:** `origin/master` at `4bed517` (2026-11-14). Line numbers
> refer to that commit. Function names are more stable than line numbers; if they
> disagree, trust the function and re-check this doc.
>
> **Status labels:** **[code]** read directly in the code at the cited location.
> **[tested]** confirmed by a run, with the commit, command or data named. Most
> runs behind this doc used 3.0.3 or `e61eeb9`-era binaries; that is stated
> where it applies. **[hypothesis]** reasoned but not tested. Hypotheses live only
> in [Open questions](#open-questions-and-hypotheses).

## Guide status

Each guide carries a status that only moves upward [code]. The enum is at
`bundle.h:28-31`; accessors are `getGuideStatus`/`setGuideStatus`
(`bundle.cpp:10`, `bundle.cpp:13`).

| Value | Name | Set where |
| --- | --- | --- |
| 0 | `GBST_UNSET` | never pulled into a bundle |
| 1 | `GBST_IN_BUNDLE` | `BundleData::getReady()`, for every kept guide, just before the bundle is queued (`bundle.h:480`) |
| 2 | `GBST_ALL_INTR_COV` | every intron has at least one spliced read (`rlink.cpp:14682`) |
| 3 | `GBST_STORED` | matched by a kept prediction, so it will be printed (e.g. `rlink.cpp:10170`, `13085`, `15866`, `21289`) |

A single-exon guide at status ≥2 becomes STORED in `build_graphs()` only when
its coverage is non-zero (`rlink.cpp:15859-15866`, from `8653184`) [code].

**Bundles pull in guides transitively** [code]. When a read extends the bundle,
the `do { … } while (cend_changed)` loop adds every guide that starts before the
new end and extends the end to the guide's end (`stringtie.cpp:1135-1149`). The
inline comment "no transitive overlap check here" refers to the single guide
test, not the loop. So overlapping transcripts of one gene always share a
bundle. A gene's guides are split between bundled and never-bundled only when
one `gene_id` labels non-overlapping loci on the same chromosome.

## The final guide pass: `writeUnbundledGuides()`

`stringtie.cpp:2149`. One pass over every guide on every chromosome
(`GRefData`, `bundle.h:38`; `crefd.rnas` are the guides). "Unbundled" is a
historical name: since `8653184` it also handles guides that were in a bundle.
It runs after all bundles are printed [code]:

- normal mode, with `-G`: `stringtie.cpp:1295`, passing `g_out` (non-NULL only with `-A`);
- merge mode with `-G`: `stringtie.cpp:1206`, **without** `gout`, so the `-A` branch never runs.

One shared skip gates two jobs [code]:

1. Skip any guide with status ≥ `GBST_ALL_INTR_COV` (`stringtie.cpp:2173`).
2. Otherwise, for `-e` or `--merge`, write a 0-coverage transcript record
   (`stringtie.cpp:2180-2219`). In default mode nothing is written to `-o`.
3. Otherwise, with `-A`, collect a 0.0 gene row into `geneabs`, **only if the
   guide's status is < `GBST_IN_BUNDLE`** (`stringtie.cpp:2226`). Rows are
   written per chromosome by `writeUnbundledGenes()` (`stringtie.cpp:2134`).

In merge mode with `--dedup-guides`, a branch before the skip suppresses
structurally identical guides (`stringtie.cpp:2158-2168`); see MERGE_MODE.md.

### #506: duplicate 0.0 gene rows (fixed)

`printResults()` writes a real `-A` row for every reference gene with a guide in
a processed bundle (`refgene`, `rlink.cpp:22131`; rows at `22985-22996`). A
bundle is printed when `pred.Count()>0 || ((eonly || geneabundance) &&
keepguides.Count()>0)` (`stringtie.cpp:1921`). The only early return that skips
gene rows is `-R` (`rlink.cpp:21918-22101`) [code].

Root cause [code + git history]:

- `8653184` changed the shared skip from "any non-zero status" (2.2.3's
  `td.in_bundle`) to ≥ STORED.
- `e7be600` (#493) changed it to ≥ ALL_INTR_COV and made `printMergeResults()`
  mark matched guides STORED, which fixed `--merge -G` printing every reference
  twice.
- A guide left at status 1 (in a bundle, not fully covered, not stored) then
  passed the skip. Its gene, already written by `printResults()`, got a second
  all-zero row with the same coordinates. #506 (also #509) is the `-A` sibling
  of #493.

The fix (`e61eeb9`, PR #510) is the `<GBST_IN_BUNDLE` test on the gene branch
only (`stringtie.cpp:2226`). It restores 2.2.3's gene-row rule and leaves the
transcript skip unchanged.

- **[tested on `e61eeb9` vs. its parent]** 0 duplicates on all repros; only
  0.0 rows were removed and no gene ID disappeared; `-o`, `-C` and `-b` were
  identical; `-p 4` matched `-p 1`; 15 regression outputs across `-e`,
  assembly and `--merge -G` were identical. Before the fix, `long_reads.bam -L
  -G human-chr19_P.gff` had 17 duplicate IDs and the SARS-CoV-2 repro
  (`--rf -G genome.gtf -A`, nf-core `sarscov2` `test.paired_end.sorted.bam`)
  duplicated `orf1ab`.
- **[code]** `-e` was never affected. In `printResults()` under `-e`
  (`rlink.cpp:22148-22157`), every non-nascent guide below STORED is printed
  with 0 coverage; status-1 guides become STORED and status-2 guides stay at 2.
  Either way the final pass skips them.
- **[code at 4bed517]** `run_tests.sh` guards this: `check_gene_abund`
  (`run_tests.sh:37-49`, called at `:117`) fails on any repeated gene ID and on
  an unexpected count of non-`STRG` genes.

**Rejected alternative.** Pre-collecting every gene with any in-bundle guide and
suppressing its zero rows was rejected. It drops the zero row of read-less
copies of a multi-locus gene ID, and it does not make rows consistent across
samples either. Both versions gave the same output on the test data
[tested on `e61eeb9`-era build].

### Gene keys and multi-locus gene IDs

- **[code]** The gene key is `gene_id`, falling back to `gene_name`
  (`stringtie.cpp:2224-2225`; `rlink.cpp:20263-20264`). In `printResults()` the
  key also carries the strand (`gid+=pred[n]->strand`, `rlink.cpp:20266`).
- **[code]** In the final pass, all never-bundled guides of one gene ID on one
  chromosome merge into **one** 0.0 row from the minimum start to the maximum
  end, with no overlap check (`stringtie.cpp:2232-2233`). A gene on both
  strands only prints a warning (`stringtie.cpp:2229-2231`). The 2.2.3 source
  has the same rule.
- **[tested on 3.0.3 and `e61eeb9`]** With gene IDs edited in
  `human-chr19_P.gff`, a three-copy ID with no reads gave one 33 Mb 0.0 row; a
  gene with one expressed and one read-less copy gave 2 rows.
- **[tested on a UCSC download]** In hg38 `refGene`, the gene symbol is the
  `gene_id`. 420 of 28,250 IDs on the main chromosomes span more than one
  non-overlapping locus on the same chromosome and strand (e.g. U2AF1 on
  chr21; tandem arrays such as CT47A1 and RNA5S1 have about 10 copies).
- **[code]** gclib's GFF reader joins records that share a transcript ID on
  the same chromosome and strand if they start within `GFF_MAX_LOCUS` = 7 Mb
  (`gclib/gff.cpp:7`, `1503`); exons farther apart raise "duplicate GFF ID"
  (`gff.cpp:1601`). On 3.0.3, a reused ID at two loci 2 Mb apart became one
  2 Mb chimeric transcript [tested on 3.0.3, merge input]. This affects `-G`
  guides too.

## `-e`, `-A` and default mode

`-e` needs `-G`; without it StringTie exits with "invalid -e usage, GFF
reference not given" (`stringtie.cpp:1691-1692`). With `--merge`, `-e` is turned
off and sets `includecov=true` instead (`stringtie.cpp:1687-1689`) [code].

| | default (`-G`, no `-e`) | `-e -G` | no `-G` |
| --- | --- | --- | --- |
| Reads used | all | only reads overlapping a guide (`stringtie.cpp:1158`) | all |
| Non-guide junctions | normal filters | all rejected in `good_junc()` (`rlink.cpp:14304-14307`) | normal filters |
| `-o` GTF | assembled transcripts only; no 0-cov guide records | every guide, 0-cov ones included | assembled only |
| `-A` reference gene rows | yes; 0.0 for genes never in a bundle | yes; 0.0 likewise | none |
| `-A` `STRG.N` rows | one per assembled locus, gene name `-`, even if it matches the annotation (`rlink.cpp:22999-23038`) | none (`if(!eonly)`, `rlink.cpp:22999`) | yes |

Row sources are [code]; the `STRG` behavior and the double counting below were
also seen in runs [tested on 3.0.3].

**How gene values are computed** (`printResults()`, `rlink.cpp:22940-22996`)
[code]:

- **Default mode:** gene coverage is strand-resolved read coverage (`bpcov`)
  over the merged reference exons, divided by their length (`rlink.cpp:22986`).
  The FPKM/TPM source (`covsum`) is the sum of coverages of predictions that
  match a guide of the gene (`t_eq`, `rlink.cpp:20260-20271`) plus leftover
  coverage on the reference exons that no matching prediction explains
  (`rlink.cpp:22979-22981`).
- **`-e`:** the read scan is skipped (`if(!eonly)`, `rlink.cpp:22949`). Gene
  coverage is the length-weighted transcript coverage over the merged exon
  length, and `covsum` is the sum of transcript coverages, so gene TPM is the
  sum of transcript TPMs. In `mix_guides.gff`, gene38337 has 12 transcripts and
  gene TPM 999999.998899 [tested on 3.0.3].
- FPKM and TPM for both transcripts and genes are computed in one place from
  the temp file: `tcov*1e9/Frag_Len` and `tcov*1e6/Cov_Sum`
  (`stringtie.cpp:1265-1291`).
- In default mode a reference-matching locus is counted twice in the table:
  once in the reference gene row and once in its `STRG` row (e.g. `STRG.1`
  equals `uc060qou.2` on the test data) [tested on 3.0.3].

`human-chr19_P.gff` uses each transcript ID as its gene ID, so it is a weak test
of gene grouping. `mix_guides.gff` has one gene with 12 isoforms.

## `-M`: multi-mapping fraction (#500, #324)

Multi-mapped alignments add `1/NH` of their count to the read and to `bpcov`
unless `-u` is given (`rlink.cpp:974`; `bpcov` uses `read_count`,
`rlink.cpp:148-152`) [code]. The `-M` cutoff (`mcov`, default 1,
`stringtie.cpp:173`, `1639-1641`) then differs by bundle strand [code]:

| Bundle | Kept if | Where |
| --- | --- | --- |
| stranded | `multi/cov <= mcov` **and** length ≥ `-m`, **or** `nolap>0` | `rlink.cpp:16063-16065` |
| neutral (no stranded reads) | `multi/cov <= mcov*(1-ERROR_PERC)`, or `guide_ovlp`, or `-R` | `rlink.cpp:15847` |

- `ERROR_PERC` is 0.1 (`rlink.h:14`). A locus covered only by NH>1,
  unspliced, unstranded reads therefore needs `-M` ≥ 1/0.9 ≈ 1.111. On the
  #324 BAM (6,247 alignments, no junctions, all multi-mapped; attached to the
  issue), `-M 1` and `-M 1.11` gave 0 transcripts, `-M 1.112` and `-M 1.12` gave
  1 [tested on 3.0.3]. This dates from the 2015 source.
- `nolap` counts same-strand guides with status ≥ ALL_INTR_COV that overlap the
  bundle (`rlink.cpp:16037-16047`).
- `guide_ovlp` only considers **single-exon** guides (the loop skips
  `exons.Count()>1`) and only tests overlap with the bundle's first node
  (`rlink.cpp:15840-15845`).
- With `-G`, a guide over the locus rescued it at the default `-M`, for both a
  one-exon and a two-exon guide, and the output took the guide's strand. A
  guide elsewhere on the contig did not rescue it [tested on 3.0.3]. The
  two-exon case cannot come from `guide_ovlp`; see Open questions.
- The help text says only "fraction of bundle allowed to be covered by
  multi-hit reads" (`stringtie.cpp:67`, `README.md:190`). A clarification of
  the 0.9 factor was drafted but never applied. #500 and #324 were closed with
  no code change.

## Opposite-strand junctions and transcripts, short reads (#498)

Places that drop an opposite-strand junction [code; reporter's data was not
available]:

1. **Same junction on both strands** (`build_graphs()`, `rlink.cpp:14778-14784`):
   a strand loses it only when it has fewer reads **and** fewer good reads.
2. **`--mix` only**: a junction whose reads all have mismatches loses its strand
   if an opposite-strand junction with some clean reads lies within `sserror`
   (25 bp, `stringtie.cpp:167`) at both ends (`rlink.cpp:14744-14765`).
3. **Splice-site consensus**, only with `-S`/`--rseq` (genome loaded at
   `stringtie.cpp:1516-1520`, copied per bundle at `968-973`): a non-guide
   junction with a non-consensus site and fewer than `DROP/ERROR_PERC` = 5 good
   reads loses its strand (`rlink.cpp:14773`, `14831`).
4. **`good_junc()`** (`rlink.cpp:14299`): besides the `-e` rule and the
   `junctionthr` minimum, it drops a junction whose reads are a tiny fraction of
   nearby coverage when that coverage does not drop across the splice site.

Transcript-level filter in `print_predcluster()` (`rlink.cpp:19132`;
loop at `19482-19555`) [code]. Predictions are visited by coverage, so `n1`
has coverage ≥ `n2`. For each overlapping pair, a prediction that matches a
guide (`t_eq`) is never removed (`rlink.cpp:19493`). Several checks run before
the strand test and apply to opposite-strand pairs too: very low coverage next
to a reference-matching `n1`, short terminal exons, and `retainedintron()`. When
the strands differ (`rlink.cpp:19514-19553`):

- `n2` is removed if it is single-exon, or if `n1` is a long-read prediction
  and `n2` is a short-read prediction with coverage < 10;
- otherwise the **higher**-coverage `n1` is removed if it is single-exon,
  matches no guide, and its coverage is < `n2` + `singlethr`;
- otherwise `n2` is removed if its structure is included in `n1`
  (`included_pred()`, `rlink.cpp:18232`);
- otherwise, for a short-read `n2` outside `--mix`, `n2` is removed if on
  every exon it holds no more than half of the total coverage.

STAR `--outSAMstrandField intronMotif` assigns a junction's strand from its
motif, so one coordinate pair cannot carry both strands. Opposite-strand
junctions that overlap forward ones must sit a few bp apart. The reply drafted
for #498 asked for `SJ.out.tab` lines and suggested EASTR to remove spurious
spliced alignments; whether it was posted is not recorded.

## Other issue answers

- **#501 (endogenous vs. exogenous transcripts).** StringTie never compares
  read or transcript sequences to tell isoforms or alleles apart. It uses
  sequence only for the `-S` splice-site check above and for long-read polyA
  detection (see LONG_READS_POLYA.md). The posted advice: add the construct's
  full transcribed sequence as an extra contig and align with minimap2
  `-ax splice`; vector UTRs help assignment. A one-base wild-type/mutant
  difference must be counted per base (e.g. `mpileup`). The issue was left
  open.
- **#499 (HISAT2 `--dta-cufflinks`).** Per the HISAT2 manual (not run), it
  implies `--dta`, limits novel splice sites to GT/AG, GC/AG and AT/AC, and
  adds `XS:A` to every spliced alignment, so its BAMs work with StringTie. A
  reply was drafted; posting is not recorded.

## Open questions and hypotheses

Untested unless marked otherwise. Each line says what would settle it.

- **[hypothesis] How a two-exon guide rescues a multi-mapped locus.**
  `evalReadAln()` returns the overlapping guide's strand in `xstrand`, which is
  copied to `alndata.strand` (`stringtie.cpp:1156-1162`). Reads then build a
  stranded bundle, where `multi/cov <= 1` passes at `-M 1`; `nolap` may also
  apply. *Settle by:* rerunning the #324 BAM with a two-exon guide on 4bed517
  and logging which branch keeps the bundle.
- **[hypothesis] `-M` never filters stranded loci at the default.** If
  `multi ≤ cov` always holds, `multi/cov <= 1` always passes. All-multi spliced
  loci can still fall below `-c`/`-s` because of the 1/NH weighting.
  *Settle by:* an all-multi spliced test locus at `-M 1` and `-M 0.5`.
- **[hypothesis] Novel isoforms and default-mode gene rows.** From the code, a
  novel isoform's reads on shared reference exons reach the gene through the
  leftover term; reads on its novel exons do not. No test assembled a novel
  isoform. *Settle by:* a guide set with one isoform removed, compared with
  the full set.
- **[hypothesis] #498 cause.** The reporter's high-coverage reverse transcripts
  are probably removed by the `print_predcluster()` strand branch. Screenshots
  (about 70 bp/pixel) showed reverse junction ends at forward exon boundaries
  (~92.1, 119.9, 141.5 kb), which suggests misalignment near exon edges or
  repeats. *Settle by:* the reporter's BAM region and `SJ.out.tab`.
- **Merged-span 0.0 rows.** The single row spanning every read-less copy of a
  multi-locus gene ID (`stringtie.cpp:2232-2233`) is long-standing and was
  never filed. *Settle by:* deciding whether rows should split per
  non-overlapping locus, then filing an issue.
- **[hypothesis] Opposite-strand single-exon collapse in `--merge`.** The
  exact-duplicate loop in `build_merge()` (`rlink.cpp:16737-16763`) compares
  junction pointers, or exact start/end for single-exon inputs, and never
  checks strand. Two stranded single-exon inputs with identical coordinates on
  opposite strands may collapse into the first one's strand. (The 3.0.3 case
  with merely overlapping single exons no longer applies; unstranded inputs
  bypass this loop, `rlink.cpp:16661-16666`.) *Settle by:* a two-record GTF
  run through `--merge`.
- **[hypothesis] Negative `-g`.** `bundledist` is `uint` (`stringtie.cpp:171`)
  and `-g` is assigned with no sign check (`stringtie.cpp:1607`). Line 1608
  then raises `runoffdist` to the wrapped value, and `stringtie.cpp:945` casts
  it back to `int`. On 3.0.3, `-g -50` kept nearby merge single exons apart.
  The current effect on bundling is unknown. *Settle by:* rejecting negative
  values, or a run with `-g -50` on test data.
- **[code] `-e -R`.** The message says `-e` is ignored, but `eonly` is never
  reset (`stringtie.cpp:1683-1685`), and it is printed only with `-v`.
  *Settle by:* deciding whether to clear `eonly` there.
