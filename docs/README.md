# Developer notes

Subsystem notes for people and coding agents working on StringTie. Each note
says which commit it was checked against and labels every claim:

- **[code]**: read in the source at the cited location
- **[tested]**: confirmed by a run; the data and command are named
- **[hypothesis]**: reasoned, not tested; these live only in each note's
  "Open questions" section

If the code and a note disagree, the code wins. Fix the note, or flag the
claim as stale.

| Note | Covers |
| --- | --- |
| [PROCESSREAD.md](PROCESSREAD.md) | Read intake on the producer thread, duplicate pooling and its index (#511) |
| [LONG_READS_POLYA.md](LONG_READS_POLYA.md) | PolyA/T tail detection, artifact screening, end anchoring, why tail-derived strand misleads on ONT cDNA (#515) |
| [STRAND_AND_JUNCTIONS.md](STRAND_AND_JUNCTIONS.md) | Read and junction strand, deletions next to introns (#514), read-order effects, `-E` |
| [GUIDES_AND_ABUNDANCE.md](GUIDES_AND_ABUNDANCE.md) | Guide status, `-e`/`-A` output, `-M` (#506, #500) |
| [BUILD_AND_TEST.md](BUILD_AND_TEST.md) | Test suites, test data, platform build problems |
| [../MERGE_MODE.md](../MERGE_MODE.md) | `--merge` workflow (maintained separately) |

The `worklogs` branch holds dated development logs, for example the polyA vs
RiboZero assembly observations. These notes record conclusions; worklogs
record how they were reached.

## Keeping notes current

- When a change invalidates a claim, update the note in the same commit.
- When re-verifying a note, update its "Verified against" commit.
- A hypothesis becomes a fact only with a cited test or code location.
