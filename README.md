# E. coli NMR Metabolomics Pipeline

A pipeline for comparing two groups of 1D ¹H NMR spectra from Bruker TopSpin,
built for an *E. coli* metabolomics study but written to run on any two-group
NMR comparison. It picks peaks directly off the raw spectra, matches each
sample's peak list against a compound reference library, and produces the
statistics and figures for a control-vs-treatment comparison.

Peak detection uses a tunable shape-based detector rather than a stock
quantification library, specifically so each candidate compound is scored
by how much of its own characteristic peak pattern actually shows up in
that sample, not just whether one peak happens to land nearby. See "Known
limitations" below for the reasoning behind that scoring.

## Requirements

**R** (tested on R 4.x):
- [ASICS](https://bioconductor.org/packages/release/bioc/html/ASICS.html)
  (Bioconductor). Only `importSpectraBruker()`/`createSpectra()` are used
  from it, to read raw Bruker spectra into R; nothing in this pipeline
  calls ASICS's own quantification (`ASICS()`) or its stock reference
  library. It is still a required dependency for that import step alone,
  since there's no separate lightweight package for reading Bruker-format
  spectra in R.
- ggplot2
- ggrepel (optional; figures fall back to base `geom_text` if it isn't
  installed, with more overlapping labels)

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ASICS")
install.packages(c("ggplot2", "ggrepel"))
```

**Python** (tested on 3.10+):
- pandas
- numpy
- matplotlib (only needed for `analyze_standard_dilution_linearity.py`'s
  plot; see "Try the pipeline first" below)

```bash
pip install pandas numpy matplotlib
```

**Python, only if you use the optional HSQC add-on** (see "Optional: 2D
[1H,13C]-HSQC support" below). Nothing else in this pipeline needs these:
- nmrglue (reads the raw 2D Bruker spectrum)
- scipy (used by the 2D peak picker)
- matplotlib (used only to save each sample's HSQC contour plot)

```bash
pip install nmrglue scipy matplotlib
```

## Installing nmrpick and nmrmatch

This pipeline's reusable logic also ships as two installable packages, one
per language, so a new project can `library(nmrpick)` / `import nmrmatch`
instead of copying loose scripts around. Nothing under `src/` breaks if
you skip this: the R runner scripts (`run_*.R`) now call `library(nmrpick)`
instead of `source()`-ing the function files directly, and the Python
scripts that got absorbed into `nmrmatch` (`peak_matcher.py`,
`match_r_feature_table.py`, `compute_fold_induction.py`, and the four
reference-library builders) are now thin shims that import from it, so
`python match_r_feature_table.py ...` etc. still work exactly as before.
Installing the packages is what makes `library(nmrpick)` and the shims
resolve, and it's the path a new project building on this pipeline should
use rather than copying `src/` wholesale.

**`nmrpick`** (R) -- baseline correction, peak detection, lineshape
refinement, cross-sample alignment, and the group-comparison statistics
(`compound_stats_table()`, `wilcoxon_exact_floor()`, etc.):

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("shrinidhi1997/Ecoli-nmr-metabolomics-pipeline", subdir = "nmrpick")
```

**`nmrmatch`** (Python) -- reference-library matching, confidence
scoring, fold-change, and the reference-library-building utilities:

```bash
pip install "git+https://github.com/shrinidhi1997/Ecoli-nmr-metabolomics-pipeline.git#subdirectory=nmrmatch"
```

or, from a local clone, editable for development:

```bash
pip install -e ./nmrmatch
```

Installing `nmrmatch` also installs its command-line tools
(`nmrmatch-match-feature-table`, `nmrmatch-fold-change`, and others --
see `nmrmatch/README.md` for the full list), which the commands below use
in place of `python match_r_feature_table.py` / `python
compute_fold_induction.py`; either form does the same thing.

Repository layout, with the two packages alongside everything else:

```
.
├── nmrpick/              R package (baseline correction, peak detection, stats)
├── nmrmatch/             Python package (matching, scoring, fold-change, reflib)
├── src/                  Driver/example/benchmark/plotting scripts + CLI shims;
│                         the two packages above are built from files in here
├── reference_library/    The two shipped reference-peak CSVs + their license
├── example_data/         Small, static, committed input for the demo pipeline
├── example_results/      Small, static, committed output for the demo pipeline
└── output/               Your own run's output -- regenerated, .gitignore'd
```

## Try the pipeline first: a simulated standard dilution series

You don't need your own Bruker data to see what this pipeline actually
produces. `example_data/` and `example_results/` ship a fully SIMULATED
(not real instrument, not real *E. coli*) standard dilution series for
L-alanine at five concentrations (0.5, 1, 2, 5, 10 mM), run end to end
through the real, unmodified detection and matching code in this src/
folder, for both the 1D and the optional 2D HSQC path. This is also a
template: once you have your own real standards at known concentrations,
the same three scripts below work on them unchanged (the 1D one needs a
new `data_dir`/`sample_group_map`-style setup for your real Bruker
folders in place of the synthetic spectra; the matcher and linearity
scripts are already generic).

**1D**, from `src/`:

```r
source("simulate_alanine_standard_dilution_1d.R")
```

```bash
nmrmatch-match-feature-table \
    --feature-table ../example_results/1d/peak_list_alanine_dilution.csv \
    --reference ../reference_library/reference_peaks_combined_ecoli_only.csv \
    --tolerance 0.008 --top-n 500 --min-confidence 0.5 \
    --out ../example_results/1d/compound_matches_alanine_dilution.csv
# (equivalently: python match_r_feature_table.py with the same flags, once
# nmrmatch is installed -- see "Installing nmrpick and nmrmatch" above)

python analyze_standard_dilution_linearity.py \
    --matches ../example_results/1d/compound_matches_alanine_dilution.csv \
    --concentrations ../example_data/1d/concentrations.csv \
    --compound-id HMDB0000161 \
    --out ../example_results/1d/alanine_dilution_linearity.csv
```

What this run actually found (not a guess -- see the printed output when
you run it yourself): L-alanine (`HMDB0000161`, already in the shipped
E. coli reference library) comes out the top match at every one of the 5
concentrations, and `matched_intensity_mean` vs. concentration fits a
line with R² = 0.999. The honest caveat underneath that: only
`find_shape_based_peaks()` is used here, not
`find_peaks_second_derivative()` -- the second detector returned several
hundred spurious peaks on this synthetic noise, a known fragility
already documented in its own docstring in `baseline_peak_picking.R`, not
a new problem introduced by this demo. And at the two lowest
concentrations (0.5, 1 mM), only 4 of alanine's 6 real reference peaks
clear the default detection threshold; all 6 are found from 2 mM up.
See `simulate_alanine_standard_dilution_1d.R`'s own header comment for
the full explanation of both.

**2D HSQC**, from `src/`:

```bash
python simulate_alanine_standard_dilution_hsqc.py
```

This builds its own tiny, clearly-labeled demo-only HSQC reference
library (`example_data/hsqc/reference_hsqc_demo_alanine.csv`), since (as
below) no real one ships by default. Both of alanine's cross-peaks were
matched at every concentration in this run; the 13C shift values used are
typical published values for L-alanine, not a live BMRB pull -- see that
script's own header for specifics.

Everything under `example_data/` and `example_results/` is intentionally
committed and static (about 57 MB total, no single file over 13 MB, well
inside GitHub's limits) -- unlike `output/`, which
is regenerated output and excluded by `.gitignore` (see "What not to
commit" below).

## Input data layout

`ASICS::importSpectraBruker()` expects, directly under one folder, one real
subfolder per sample, each containing a single numbered experiment folder
(the Bruker "expno") with that experiment's own `pdata/<procno>` inside it:

```
<data_dir>/
  sample_01/
    1/              <- expno
      pdata/
        1/          <- procno; the actual processed spectrum files live here
  sample_02/
    1/
      pdata/
        1/
```

If your raw experiment folders already sit directly in that shape, point
`data_dir` (in `config_control_vs_sample.R`) straight at their parent folder and
skip the rest of this section.

### Staging real Bruker data without copying it

More often, your real experiment folders are scattered across other project
folders, each ending at its own `pdata/<procno>`. Rather than copying that
data, stage it with a junction (Windows) or symlink (macOS/Linux) per
sample, pointing at the folder *containing* the `pdata` folder (i.e. the
expno folder), not at `pdata` itself. `importSpectraBruker()` picks the
numerically-largest subfolder name it finds and expects it to be a real
experiment number; junctioning straight at `pdata` breaks that, since none
of `pdata`'s own subfolders are numbers.

**Windows**, one sample:

```bat
mkdir "<data_dir>\sample_01"
mklink /J "<data_dir>\sample_01\1" "<path to the real experiment folder>\1"
```

**macOS/Linux**, one sample:

```bash
mkdir -p "<data_dir>/sample_01"
ln -s "<path to the real experiment folder>/1" "<data_dir>/sample_01/1"
```

Repeat per sample. The outer folder name (`sample_01` above) becomes that
sample's column name in every output table and figure, regardless of what
the real folder was called; the inner link's name only needs to be the real
experiment number.

To rebuild the whole staging folder from scratch for a new sample set
(safe, since it holds only junctions/symlinks, never your real data):

```bat
rmdir /S /Q "<data_dir>"
mkdir "<data_dir>"
```

then repeat the `mklink`/`ln -s` step per sample.

## Configuration: the one file to edit

**`src/config_control_vs_sample.R`** is the single file every runner script
sources, and the only file you should need to edit for a new dataset or a
different condition. It has four parts:

1. `data_dir`, the staging folder built above.
2. `sample_group_map`, a named vector mapping each sample's folder name to
   one of exactly two group labels of your choice (`"control"`/`"treatment"`,
   `"wild_type"`/`"mutant"`, anything).
3. Detection-sensitivity knobs for the peak picker (`prominence_mult`,
   `min_height_mult`, and related settings), documented where they're used
   in `baseline_peak_picking.R`.
4. `out_base`, where output CSVs and figures are written.

Two further settings live outside this file, as command-line flags on
`match_r_feature_table.py`/`nmrmatch-match-feature-table` rather than in R:
`--tolerance` (the ppm window for matching a peak to a reference shift) and
`--min-confidence` (the confidence floor for keeping a candidate compound).
Raise `--tolerance` if your instrument shows more shift drift than the
default; see that script's own `--help` for the rest.

## Running the pipeline

From `src/`:

```r
source("run_peak_picking_control_vs_sample.R")
```

This writes `peak_list_control_vs_sample.csv` under `out_base`, and prints the
exact next command to run, filled in with your actual sample names. Then,
from a shell:

```bash
nmrmatch-match-feature-table \
    --feature-table ../output/control_vs_sample/peak_list_control_vs_sample.csv \
    --reference <path to your reference peak library CSV> \
    --tolerance 0.008 \
    --top-n 500 \
    --min-confidence 0.88 \
    --out ../output/control_vs_sample/compound_matches_control_vs_sample.csv
```

(equivalently, `python match_r_feature_table.py` with the same flags, if
you'd rather not install `nmrmatch` -- see "Installing nmrpick and
nmrmatch" above.)

Then, from R, for the full statistics and figures:

```r
source("run_stats_and_plots_control_vs_sample.R")
```

or, for just a plain fold-change CSV with no R session and no figures:

```bash
nmrmatch-fold-change \
    --matches ../output/control_vs_sample/compound_matches_control_vs_sample.csv \
    --group control:sample_01,sample_02,sample_03 \
    --group treatment:sample_04,sample_05,sample_06 \
    --out ../output/control_vs_sample/fold_induction_control_vs_sample.csv
```

### The reference peak library

`--reference` needs a CSV with `compound_id`, `compound_name`, `nucleus`,
and `shift_ppm` columns (`intensity` and `multiplicity` are optional), one
row per characteristic peak per compound. Two are included under
`reference_library/`:

- `reference_peaks_combined.csv` (4,109 compounds, BMRB + HMDB) is the full
  library, not filtered to any one organism.
- `reference_peaks_combined_ecoli_only.csv` (580 compounds) is that same
  library filtered to compounds confirmed as *E. coli* metabolites against
  ECMDB. This is the one the example commands above point `--reference` at.

Both already have an `is_ecoli_metabolite` column from that ECMDB check, so
you can re-filter `reference_peaks_combined.csv` for a different organism's
curated database the same way, once you have an equivalent cross-reference
for it.

The scripts that built these (also in `src/`, and packaged as
`nmrmatch.reflib` -- see "Installing nmrpick and nmrmatch" above) are
included for provenance and for extending the library, not because you
need to re-run them to use this pipeline:

- `fetch_nmr_reference.py` and `build_bmrb_library.py` pull real peak data
  from ECMDB and BMRB's public APIs (`nmrmatch-fetch-reference` /
  `nmrmatch-build-bmrb` once `nmrmatch` is installed).
- `merge_reference_libraries.py` combines a BMRB-built and an HMDB-built
  library into one shared column set (`nmrmatch-merge-libraries`).
- `filter_reference_to_ecoli.py` is the ECMDB cross-reference filter that
  produced the `_ecoli_only` file above from the full one; this is the
  "peak filter" step, run as `python filter_reference_to_ecoli.py
  --reference reference_library/reference_peaks_combined.csv --out
  reference_library/reference_peaks_combined_ecoli_only.csv` (or
  `nmrmatch-filter-ecoli` with the same flags).

One gap worth knowing about: `merge_reference_libraries.py` and
`filter_reference_to_ecoli.py` both also reference an HMDB-side builder
and an ECMDB compound-identity builder, which aren't included here, so the
HMDB half of `reference_peaks_combined.csv` isn't independently
reproducible from just this repository, only the BMRB half is. The data
files themselves are already complete either way. These builder scripts
also weren't re-run in this pass (they call live BMRB/ECMDB APIs this
environment can't reach), so give them a quick test run yourself before
relying on them against a fresh pull.

**These two CSVs carry their own license, separate from this repository's
code.** Part of their content comes from HMDB, which is released under
CC BY-NC 4.0, non-commercial use only, with commercial use requiring
HMDB's permission. See `reference_library/DATA_LICENSE.md` before you
publish or reuse this data, since it is not covered by the MIT license
in this repository's root.

## Optional: 2D [1H,13C]-HSQC support

Everything above is 1D ¹H only. If you also collected a 2D [1H,13C]-HSQC
spectrum for a sample, `hsqc_analysis.py` and `hsqc_matcher.py` add a
separate, optional path for it. Neither file is imported by the 1D
pipeline, so skipping this section entirely leaves the rest of the
pipeline untouched.

What it does, per sample:

1. `hsqc_analysis.py` reads the raw 2D Bruker HSQC spectrum (a `pdata`
   folder containing `2rr`, `procs`, and `proc2s`) and picks cross-peaks,
   excluding the water/solvent stripe on the 1H axis.
2. `hsqc_matcher.py` matches those cross-peaks against an HSQC reference
   library and ranks candidate compounds by confidence, the same
   coverage/evidence/shift-score approach `nmrmatch`'s `score_compound()`
   uses for 1D matching (`score_compound_2d`, imported from the installed
   `nmrmatch` package rather than a same-directory `utils.py`).

```bash
pip install nmrglue scipy matplotlib   # HSQC-only dependencies, see Requirements
python hsqc_analysis.py --sample-id sample_treatment_01 \
    --path /path/to/your/staged/bruker/samples/sample_treatment_01-HSQC/pdata/1 \
    --out-dir ../output
python build_bmrb_library.py --hsqc --out ../output/reference_peaks_hsqc_bmrb_library.csv
# (equivalently: nmrmatch-build-bmrb --hsqc ..., once nmrmatch is installed)
python hsqc_matcher.py \
    --reference ../output/reference_peaks_hsqc_bmrb_library.csv \
    --experimental ../output/sample_treatment_01_hsqc_peaks.csv \
    --out ../output/match_results_hsqc.csv
```

Two things worth knowing before you rely on this:

- No HSQC reference library is shipped in `reference_library/`. Unlike
  the 1D CSVs, you build this one yourself with the `build_bmrb_library.py
  --hsqc` command above, which needs live access to BMRB's API.
- The 2D peak-picking math and the HSQC matching/scoring logic were both
  tested against synthetic data while adding this (known cross-peaks
  correctly recovered, water-region exclusion working, confidence scoring
  correctly favoring a compound with more matched cross-peaks). The raw
  Bruker-file reading in `read_bruker_2d()` was not tested against an
  actual 2D instrument export, since none was available in this
  environment, only unit-tested against mocked nmrglue output. Run it
  against one of your own real HSQC folders and check the output plot and
  peak list before trusting it for the paper.

This repo also has no runner script wiring HSQC into
`config_control_vs_sample.R`'s per-sample loop the way the 1D files are
wired together; run the two scripts above by hand, per sample, for now.

## What each file does

| File | What it does |
|---|---|
| `nmrpick/` | **Installable R package** (see "Installing nmrpick and nmrmatch" above) wrapping `peak_picking_functions.R`, `baseline_peak_picking.R`, `read_bruker_1d_direct.R`, and the group-comparison statistics factored out of `run_stats_and_plots_control_vs_sample.R` (`build_compound_matrix()`, `compound_stats_table()`, `wilcoxon_exact_floor()`). `library(nmrpick)` replaces `source()`-ing these files directly; all runner scripts below now do this. |
| `nmrmatch/` | **Installable Python package** (see "Installing nmrpick and nmrmatch" above) wrapping `utils.py`, `peak_matcher.py`, `compute_fold_induction.py`, `match_r_feature_table.py`, and the four reference-library-builder scripts, plus console-script entry points (`nmrmatch-match-feature-table`, `nmrmatch-fold-change`, etc.). |
| `config_control_vs_sample.R` | The one file to edit for a new dataset: data location, group assignment, detection knobs, output location. |
| `peak_picking_functions.R` | Reduces a raw TopSpin peak-picking export to a smaller set of real peaks, then aligns that list across samples. An older, simpler alternative to `baseline_peak_picking.R`'s detector. Packaged as `nmrpick::pick_strong_peaks()`/`merge_peaks_across_samples()`/`merge_peaks_reference_guided()`/`write_metaboanalyst_csv()`; kept here unchanged as the file the package was built from. |
| `baseline_peak_picking.R` | The main peak detector: adaptive baseline fitting (airPLS), a valley-to-valley prominence-based peak finder, optional Lorentzian lineshape refinement for true peak area, and reference-guided merging of peak lists across samples. Packaged as `nmrpick::fit_airpls_baseline()`/`find_shape_based_peaks()`/`find_peaks_second_derivative()`/`refine_peaks_lorentzian()`; kept here unchanged as the file the package was built from. |
| `run_peak_picking_control_vs_sample.R` | Runs peak picking end to end: imports spectra, fits baselines once per sample, detects peaks with two complementary detectors, merges across samples, and writes the feature table the Python matcher reads. Now calls `library(nmrpick)` instead of `source()`-ing the two files above. |
| `utils.py` | The compound-scoring core: for one candidate compound, matches its reference peaks against a sample's experimental peaks and combines coverage, the number of matched peaks, shift accuracy, and how prominent the matched peaks are in that sample's own spectrum into a single confidence score. Packaged as `nmrmatch.scoring`/`score_compound()`; kept here unchanged as the file the package was built from. |
| `peak_matcher.py` | Thin shim: imports `match_sample()`/`match_all_samples()` from the installed `nmrmatch.matching` package and runs its `main()`. `python peak_matcher.py ...` still works exactly as before; `nmrmatch-match` is the equivalent console script. |
| `match_r_feature_table.py` | Thin shim over `nmrmatch.cli.match_feature_table`: reshapes the wide feature table `run_peak_picking_control_vs_sample.R` writes into the long format the matcher expects, then runs it. This is the script (or `nmrmatch-match-feature-table`) to call after peak picking. |
| `compute_fold_induction.py` | Thin shim over `nmrmatch.fold_change`: turns the matcher's per-sample output into a per-compound fold-change table between your two groups, with no R session and no figures. `nmrmatch-fold-change` is the equivalent console script. |
| `run_stats_and_plots_control_vs_sample.R` | Reads the matcher's output and produces the statistics and figures: per-compound Welch t-test, Wilcoxon rank-sum, BH-FDR correction, fold-change, a volcano plot, a group-averaged heatmap, individual-replicate boxplots, a confidence-distribution histogram, a detected-compound-count bar chart, sample correlation, PCA, and a stacked relative-composition chart. The per-compound stats now run through `nmrpick::build_compound_matrix()`/`compound_stats_table()` rather than inline code. |
| `fetch_nmr_reference.py` | Thin shim over `nmrmatch.reflib.fetch`. Pulls real NMR peak data from ECMDB for the reference library. See "The reference peak library" above. |
| `build_bmrb_library.py` | Thin shim over `nmrmatch.reflib.build_bmrb`. Pulls real NMR peak data from BMRB's metabolomics collection for the reference library. |
| `merge_reference_libraries.py` | Thin shim over `nmrmatch.reflib.merge`. Combines a BMRB-built and an HMDB-built reference library into one shared column set. |
| `filter_reference_to_ecoli.py` | Thin shim over `nmrmatch.reflib.filter_ecoli`. The peak-filter step: keeps only the rows confirmed as real *E. coli* metabolites via an ECMDB cross-reference, dropping everything else. |
| `peak_picking_2d.py` | Optional. Reads a raw 2D [1H,13C]-HSQC Bruker spectrum and picks cross-peaks. Only used by `hsqc_analysis.py`; the 1D pipeline never imports it. Not packaged (2D support is a smaller, separate add-on; see "Optional: 2D HSQC support" above). |
| `hsqc_analysis.py` | Optional. Per-sample HSQC pipeline: reads the 2D spectrum via `peak_picking_2d.py`, picks cross-peaks, saves a peak-list CSV and a contour plot. See "Optional: 2D HSQC support" above. Not packaged. |
| `hsqc_matcher.py` | Optional. Matches `hsqc_analysis.py`'s picked cross-peaks against an HSQC reference library using `score_compound_2d` (imported from the installed `nmrmatch` package), ranking candidates by confidence the same way `peak_matcher.py` does for 1D. The 2D matching logic itself (`match_sample_2d`/`match_all_samples_2d`) isn't packaged, only the scoring function it calls. |
| `simulate_alanine_standard_dilution_1d.R` | Builds the simulated 1D standard dilution series under "Try the pipeline first" above: synthetic L-alanine spectra at 5 concentrations, run through the real baseline/peak-picking code. |
| `analyze_standard_dilution_linearity.py` | Generic: checks whether a matched compound's intensity scales linearly with known concentration across a dilution series. Ships with the simulated demo above, but reusable as-is on real standards data. |
| `simulate_alanine_standard_dilution_hsqc.py` | The 2D HSQC counterpart to `simulate_alanine_standard_dilution_1d.R`: synthetic cross-peaks at 5 concentrations, run through the real `pick_peaks_2d()`/`hsqc_matcher.py` code, plus a small demo-only HSQC reference library. |

Note: the benchmark comparing this pipeline's own peak detection against
ASICS's and BATMAN's stock quantification is run with separate standalone
scripts that are not part of this repository.

## Where to change things for a different condition or dataset

Almost everything is `config_control_vs_sample.R`: point `data_dir` at your
staged samples, list your samples and their two groups in
`sample_group_map`, adjust the detection knobs if your spectra are cleaner
or noisier than the dataset this pipeline was tuned on, and set `out_base`
to wherever you want results written. No runner script should need editing.

Outside that file, the only dataset-specific choices are command-line flags
on `match_r_feature_table.py`/`nmrmatch-match-feature-table` (`--tolerance`
for shift drift, `--min-confidence` for how permissive to be about weak
matches) and the reference library CSV path itself, which you supply.

## Known limitations

- **Statistical power at small sample sizes.** With `n` replicates per
  group, the exact two-sided Wilcoxon rank-sum test can reach at best
  `p = 2 / choose(2n, n)`, regardless of how extreme the actual difference
  is. At `n=3` that floor is 0.1, above the conventional 0.05 threshold, so
  no compound can clear raw Wilcoxon significance at all at that sample
  size. The Welch t-test used elsewhere in this pipeline doesn't have that
  floor, but is far more sensitive to a single noisy replicate with few
  samples per group. Treat any single-compound p-value from a small study
  as a lead for a follow-up experiment, not a standalone result, and always
  read the BH-FDR-adjusted q-value column over the raw p-value.
- **Peak height is not peak area.** By default (`use_lorentzian_refinement
  <- FALSE`), the detected "intensity" is the baseline-corrected peak
  height, which is fast but not proportional to molar amount when linewidth
  varies between samples (shimming, viscosity, temperature). Set
  `use_lorentzian_refinement <- TRUE` for real Lorentzian-fitted peak area
  instead, at the cost of a per-peak curve fit.
- **Confidence scoring rewards evidence, not just ratio.** A compound's
  match confidence combines how large a fraction of its reference peaks
  matched (`coverage`), the absolute number of peaks that matched
  (`evidence_score`), how close those matches were in shift (`shift_score`),
  and how prominent the matched peaks are in the sample's own spectrum
  (`peak_prominence`). Coverage alone cannot tell a compound matched on 2 of
  its 2 reference peaks apart from one matched on 27 of its 28, even though
  the second is far stronger evidence; see `score_compound()`'s docstring
  (`utils.py`, packaged as `nmrmatch.scoring`) for the full rationale and
  a worked example.
- **Exact Wilcoxon's small-n floor applies here too.** `nmrpick::
  wilcoxon_exact_floor(n1, n2)` makes the first bullet above a callable
  fact rather than something to recompute by hand; `compound_stats_table()`
  prints it for the group sizes actually used in a given run.

## What not to commit

This repository is code plus the two reference-library CSVs under
`reference_library/` (1-8 MB each, well inside GitHub's limits). Keep the
following out of version control (the included `.gitignore` already
excludes the generated ones):

- **Raw NMR data.** Bruker experiment folders are usually large and are
  often subject to a separate data-deposition requirement for publication
  (e.g. a repository like [MetaboLights](https://www.ebi.ac.uk/metabolights/)
  or [OSF](https://osf.io/)); link to that deposition from this repository
  instead of committing the data itself.
- **A larger or different reference library.** If you regenerate or extend
  `reference_library/reference_peaks_combined.csv` and it grows well past
  its current size, check it against GitHub's 100 MB per-file hard limit;
  anything close to that, or a large intermediate file such as a full
  ASICS pure-spectra table, belongs in [Git LFS](https://git-lfs.com/) or
  an external repository (Zenodo, OSF), linked from this README instead.
- **Generated output.** Everything under `output/` (CSVs, figures) is
  reproducible from the code and your data; it does not need to be
  committed, and the `.gitignore` already excludes it. This is different
  from `example_data/` and `example_results/` (the simulated demo under
  "Try the pipeline first" above), which are small, static, and meant to
  be committed.
- **Absolute local file paths and real sample names from your own project.**
  Double-check `config_control_vs_sample.R` and any shell history you paste
  into an issue or discussion before pushing; this repository's own scripts
  no longer contain any, but a fresh config file you edit locally for your
  own data will.
- **Credentials of any kind.** Not used by this pipeline, but worth a
  one-line reminder if you add any integration later.

## License

This repository has two different licenses, covering two different
things, and it's worth keeping them straight before you publish.

**The code in `src/`, `nmrpick/`, and `nmrmatch/`** is covered by the MIT
`LICENSE` file at the repo root (each package directory also carries its
own copy, as its language's packaging convention expects). MIT lets
anyone use, copy, modify, and redistribute the code,
including commercially, as long as they keep the original copyright
notice; it disclaims warranty, so you aren't liable if something in the
code causes a problem for someone else's use of it. It does not require
anyone who builds on the code to open-source their own changes, unlike
GPL. The copyright line in `LICENSE` is already filled in under Shrinidhi
Bhat, DMST, Indian Institute of Technology Madras; check only whether
your target journal asks for a different license, some do. An unlicensed
public repository is, by default, "all rights reserved", and others
cannot legally reuse the code even if they can read it, so some license
file is worth having either way.

Writing MIT code that calls ASICS, a GPL (>= 2)-licensed Bioconductor
package, is standard, widespread practice in R: your own code stays under
whatever license you choose, since you aren't distributing ASICS's
source yourself, just calling a separately-installed package at runtime.
This is the generally accepted reading of how GPL applies to calling code
rather than embedded code, not a settled point of law, so if your
institution has a technology-transfer or legal office, this is worth a
quick check with them given the paper this accompanies.

**The two CSVs in `reference_library/`** are not covered by that MIT
license. Part of their content comes from HMDB under CC BY-NC 4.0, a
non-commercial-only term, so they carry separate, more restrictive terms
of their own; see `reference_library/DATA_LICENSE.md` for the specifics
and the actual citations each source asks for.

## Citation

If you use this pipeline, please cite the accompanying paper (add citation
details here once available) alongside
[ASICS](https://bioconductor.org/packages/release/bioc/html/ASICS.html)
itself.

If your run also touches the reference library in `reference_library/`,
two more citations apply, each required by the source database itself,
not just good practice:

- **HMDB.** Cite the current HMDB paper, listed on
  [hmdb.ca/about](https://www.hmdb.ca/about), whenever you use data drawn
  from it (the About page lists which paper is current). See
  `reference_library/DATA_LICENSE.md` for the non-commercial term this
  data carries.
- **BMRB.** Acknowledge BMRB in any publication that uses it, per
  [bmrb.io's aims and policies](https://bmrb.io/bmrb/aims_and_policies.shtml),
  which asks for this explicitly.

To cite this pipeline's code itself, see `CITATION.cff` at the repo root,
or cite directly: Shrinidhi Bhat, DMST, Indian Institute of Technology
Madras.
