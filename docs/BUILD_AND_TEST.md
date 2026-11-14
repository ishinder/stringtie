# Building and testing

> **Verified against:** `origin/master` at `4bed517` (2026-11-14).
> Status labels: **[code]** read in the repo, **[tested]** confirmed by a run,
> **[hypothesis]** untested.

## Test suites [code]

`make test` (`Makefile:181`) runs `./run_tests.sh` and then `./test_merge_policy.sh`.
`make test-merge` (`Makefile:184`) runs only the merge suite.

**`run_tests.sh`** downloads `gpertea/stringtie-testdata` at the tag pinned in
the script (`run_tests.sh:31`). It runs 11 output comparisons:

- short reads; short reads with super-reads; short reads `-G`
- `-L`; `-L -G`
- `--mix`; `--mix -G`
- `-N -G`
- `--nasc -G`
- `--ptf`, short and long reads

`check_gene_abund()` fails if a gene ID repeats in an `-A` table (#506, `4b040d5`).
The script exits non-zero on any mismatch.

**`test_merge_policy.sh`** builds small GTFs with `mkgtf` and runs `--merge` through
`fx NAME args…`. Its helpers are `has_chain`, `no_chain`, `tx_count` and
`known_fail`. A `known_fail` is a documented defect: it passes while the command
fails, and stops the suite once a fix makes it succeed, so the check gets
promoted to an assertion. It supports `KEEP_TMP=1` and `STRINGTIE_BIN=`.

**Test data:** the `gpertea/stringtie-testdata` repo (about 1.7 MB of BAMs plus the
11 MB `human-chr19_P.gff`). Its `e_test/` fixture (`missed_regions.bam`,
`sel_refs.gtf`) is unused and could back an `-e` zero-coverage test.
`tests/README.md` is stale: it points at the old `test_data` tarball URL
(`tests/README.md:5`) and lists 9 tests.

## Known problems

- **[tested] Fresh builds fail on Apple Silicon.** `htslib/Makefile:38` puts
  `-msse4.1 -mssse3 -mpopcnt` in the global `CFLAGS` (since `3875d93`, the htslib
  1.21 upgrade). arm64 clang rejects these with `unsupported option '-msse4.1'`.
  Checkouts that already have a `libhts.a` don't hit it.
  A fix tested in a scratch clone (clean build, 11/11 tests) adds the
  flags only on x86:
  ```make
  CFLAGS   = -g -Wall -O2 -fvisibility=hidden -D HAVE_LIBDEFLATE -D HAVE_LIBLZMA
  ifneq (,$(filter x86_64 amd64 i386 i686,$(shell uname -m)))
  CFLAGS  += -msse4.1 -mssse3 -mpopcnt
  endif
  ```
  An alternative tested earlier (#504 work) removes the global flags entirely.
  htscodecs then compiles only its runtime-dispatched SSE4 file with those
  flags. That needs `htscodecs.mk` to be cleared on clean: `build_lib.sh clean`
  and `make clean` leave it in place, and a stale copy fails with
  `'_mm_mullo_epi32' requires target feature 'sse4.1'`. That stale file is the
  likely, undocumented, reason the global flags were added.
- **[tested] The merge suite fails on case-insensitive filesystems (default
  macOS).** Outputs like `$t/G3.gtf` (`fx G3 …`, `test_merge_policy.sh:327`)
  overwrite inputs like `$t/g3.gtf`. 16 name pairs collide (af/Af, g3/G3,
  h3/H3, …). Renaming the lowercase inputs to `*_in.gtf` fixed it in a scratch
  copy.
- **[code, from #475] The output comparison is byte-exact.** Long-read and `--mix`
  tests can differ in the last digit of `cov`/FPKM/TPM across platforms. A
  numeric-tolerance compare has been proposed, not written.
- **[tested] With `-p >1`, transcript IDs can differ between runs.** Compare
  structures, not IDs, when testing multi-threaded output.

## Time-savers

- `make release` takes roughly 6-15 s and `make test` under 2 s on an M-series Mac.
- The htslib build downloads libdeflate, bzip2 and xz into `htslib/xlibs` at build
  time. CI needs network access, or a cache of that directory.
  `prep_offline_source.sh` bundles them for offline builds.
- Top-level `make clean` does not clean htslib; use `make clean-htslib` or `clean-all`.
- Rebuild before trusting results; a stale `./stringtie` binary has misled
  past debugging.
- `--merge` creates its temp dir next to the `-o` path, so `-o /dev/null` fails.
- `gh issue view N --comments` sometimes prints nothing; use
  `gh issue view N --json number,title,body,comments`.

## CI (status as of 2026-11-14)

Geo agreed to GitHub Actions that run on pull requests and by hand only, never on
pushes to master, with heavy tests kept local. A draft workflow exists but has
not been submitted.
