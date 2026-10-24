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

```bash
pip install pandas numpy
```

**Python, only if you use the optional HSQC add-on** (see "Optional: 2D
[1H,13C]-HSQC support" below). Nothing else in this pipeline needs these:
- nmrglue (reads the raw 2D Bruker spectrum)
- scipy (used by the 2D peak picker)
- matplotlib (used only to save each sample's HSQC contour plot)

```bash
pip install nmrglue scipy matplotlib
```

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
`match_r_feature_table.py` rather than in R: `--tolerance` (the ppm window
for matching a peak to a reference shift) and `--min-confidence` (the
confidence floor for keeping a candidate compound). Raise `--tolerance` if
your instrument shows more shift drift than the default; see that script's
own `--help` for the rest.

## Running the pipeline

From `src/`:

```r
source("run_peak_picking_control_vs_sample.R")
```

This writes `peak_list_control_vs_sample.csv` under `out_base`, and prints the
exact next command to run, filled in with your actual sample names. Then,
from a shell:

```bash
python match_r_feature_table.py \
    --feature-table ../output/control_vs_sample/peak_list_control_vs_sample.csv \
    --reference <path to your reference peak library CSV> \
    --tolerance 0.008 \
    --top-n 500 \
    --min-confidence 0.88 \
    --out ../output/control_vs_sample/compound_matches_control_vs_sample.csv
```

Then, from R, for the full statistics and figures:

```r
source("run_stats_and_plots_control_vs_sample.R")
```

or, for just a plain fold-change CSV with no R session and no figures:

```bash
python compute_fold_induction.py \
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

The scripts that built these (also in `src/`) are included for provenance
and for extending the library, not because you need to re-run them to use
this pipeline:

- `fetch_nmr_reference.py` and `build_bmrb_library.py` pull real peak data
  from ECMDB and BMRB's public APIs.
- `merge_reference_libraries.py` combines a BMRB-built and an HMDB-built
  library into one shared column set.
- `filter_reference_to_ecoli.py` is the ECMDB cross-reference filter that
  produced the `_ecoli_only` file above from the full one; this is the
  "peak filter" step, run as `python filter_reference_to_ecoli.py
  --reference reference_library/reference_peaks_combined.csv --out
  reference_library/reference_peaks_combined_ecoli_only.csv`.

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
   coverage/evidence/shift-score approach `utils.py` uses for 1D matching.

```bash
pip install nmrglue scipy matplotlib   # HSQC-only dependencies, see Requirements
python hsqc_analysis.py --sample-id sample_treatment_01 \
    --path /path/to/your/staged/bruker/samples/sample_treatment_01-HSQC/pdata/1 \
    --out-dir ../output
python build_bmrb_library.py --hsqc --out ../output/reference_peaks_hsqc_bmrb_library.csv
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
| `config_control_vs_sample.R` | The one file to edit for a new dataset: data location, group assignment, detection knobs, output location. |
| `peak_picking_functions.R` | Reduces a raw TopSpin peak-picking export to a smaller set of real peaks, then aligns that list across samples. An older, simpler alternative to `baseline_peak_picking.R`'s detector. |
| `baseline_peak_picking.R` | The main peak detector: adaptive baseline fitting (airPLS), a valley-to-valley prominence-based peak finder, optional Lorentzian lineshape refinement for true peak area, and reference-guided merging of peak lists across samples. |
| `run_peak_picking_control_vs_sample.R` | Runs peak picking end to end: imports spectra, fits baselines once per sample, detects peaks with two complementary detectors, merges across samples, and writes the feature table the Python matcher reads. |
| `utils.py` | The compound-scoring core: for one candidate compound, matches its reference peaks against a sample's experimental peaks and combines coverage, the number of matched peaks, shift accuracy, and how prominent the matched peaks are in that sample's own spectrum into a single confidence score. |
| `peak_matcher.py` | Runs `utils.py`'s scoring across every sample and every candidate compound in a reference library, ranking candidates by confidence; a direct command-line entry point against a generic long-format experimental-peaks CSV. |
| `match_r_feature_table.py` | Reshapes the wide feature table `run_peak_picking_control_vs_sample.R` writes into the long format `peak_matcher.py` expects, then runs the matcher. This is the script to call after peak picking. |
| `compute_fold_induction.py` | Turns the matcher's per-sample output into a per-compound fold-change table between your two groups, with no R session and no figures. |
| `run_stats_and_plots_control_vs_sample.R` | Reads the matcher's output and produces the statistics and figures: per-compound Welch t-test, Wilcoxon rank-sum, BH-FDR correction, fold-change, a volcano plot, a group-averaged heatmap, individual-replicate boxplots, a confidence-distribution histogram, a detected-compound-count bar chart, sample correlation, PCA, and a stacked relative-composition chart. |
| `fetch_nmr_reference.py` | Pulls real NMR peak data from ECMDB for the reference library. See "The reference peak library" above. |
| `build_bmrb_library.py` | Pulls real NMR peak data from BMRB's metabolomics collection for the reference library. |
| `merge_reference_libraries.py` | Combines a BMRB-built and an HMDB-built reference library into one shared column set. |
| `filter_reference_to_ecoli.py` | The peak-filter step: keeps only the rows confirmed as real *E. coli* metabolites via an ECMDB cross-reference, dropping everything else. |
| `peak_picking_2d.py` | Optional. Reads a raw 2D [1H,13C]-HSQC Bruker spectrum and picks cross-peaks. Only used by `hsqc_analysis.py`; the 1D pipeline never imports it. |
| `hsqc_analysis.py` | Optional. Per-sample HSQC pipeline: reads the 2D spectrum via `peak_picking_2d.py`, picks cross-peaks, saves a peak-list CSV and a contour plot. See "Optional: 2D HSQC support" above. |
| `hsqc_matcher.py` | Optional. Matches `hsqc_analysis.py`'s picked cross-peaks against an HSQC reference library using `utils.py`'s `score_compound_2d`, ranking candidates by confidence the same way `peak_matcher.py` does for 1D. |

## Where to change things for a different condition or dataset

Almost everything is `config_control_vs_sample.R`: point `data_dir` at your
staged samples, list your samples and their two groups in
`sample_group_map`, adjust the detection knobs if your spectra are cleaner
or noisier than the dataset this pipeline was tuned on, and set `out_base`
to wherever you want results written. No runner script should need editing.

Outside that file, the only dataset-specific choices are command-line flags
on `match_r_feature_table.py` (`--tolerance` for shift drift,
`--min-confidence` for how permissive to be about weak matches) and the
reference library CSV path itself, which you supply.

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
  the second is far stronger evidence; see `utils.py`'s `score_compound()`
  docstring for the full rationale and a worked example.

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
  committed, and the `.gitignore` already excludes it.
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

**The code in `src/`** is covered by the MIT `LICENSE` file at the repo
root. MIT lets anyone use, copy, modify, and redistribute the code,
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
