# nmrmatch

Reference-library matching, confidence scoring, and fold-change computation
for 1D NMR metabolomics. Companion package to
[`nmrpick`](../nmrpick), which handles baseline correction, peak detection,
and group-comparison statistics.

## Install

From a local clone (editable install, for development):

```bash
pip install -e ./nmrmatch
```

Directly from GitHub:

```bash
pip install "git+https://github.com/shrinidhi1997/Ecoli-nmr-metabolomics-pipeline.git#subdirectory=nmrmatch"
```

## What's in it

- `nmrmatch.scoring` - tolerance-based peak matching and the compound
  confidence score (`score_compound`, plus the 2D NMR variants).
- `nmrmatch.matching` - matches one or many experimental peak lists against
  a reference library (`match_sample`, `match_all_samples`).
- `nmrmatch.fold_change` - group-wise fold-induction from a compound-matches
  table.
- `nmrmatch.reflib` - one-off utilities used to build/curate the reference
  peak library itself (BMRB bulk pull, E. coli filtering, library merging,
  an ECMDB scrape template). Not part of the per-sample matching pipeline.
- `nmrmatch.cli` - the end-to-end driver that matches an R-exported feature
  table against the reference library.

## Usage

```python
from nmrmatch import score_compound, match_all_samples
```

Each former standalone script is also installed as a command-line tool:

| Command | Originally |
|---|---|
| `nmrmatch-match` | `peak_matcher.py` |
| `nmrmatch-match-feature-table` | `match_r_feature_table.py` |
| `nmrmatch-fold-change` | `compute_fold_induction.py` |
| `nmrmatch-fetch-reference` | `fetch_nmr_reference.py` |
| `nmrmatch-filter-ecoli` | `filter_reference_to_ecoli.py` |
| `nmrmatch-merge-libraries` | `merge_reference_libraries.py` |
| `nmrmatch-build-bmrb` | `build_bmrb_library.py` |

Run any of them with `--help` for its options.
