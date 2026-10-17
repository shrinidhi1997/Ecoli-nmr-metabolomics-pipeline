"""
utils.py
Shared low-level helpers for tolerance-based NMR peak matching.

Nothing in this file talks to the network or to ECMDB. It only operates on
plain lists/arrays of chemical shifts, so it can be unit tested with made-up
numbers and reused regardless of where the reference peak table came from.
"""

from __future__ import annotations
import numpy as np


def gaussian_match_weight(delta_ppm: float, tolerance_ppm: float) -> float:
    """
    Convert a ppm distance between an experimental peak and a reference peak
    into a match weight between 0 and 1.

    A hit exactly on-shift gets weight 1.0. A hit right at the edge of the
    tolerance window gets a small but non-zero weight. Anything further than
    `tolerance_ppm` should not be passed in at all (callers filter that out
    with `candidate_matches` below).

    Using a smooth falloff (rather than a flat 1/0 inside the window) avoids
    ties between "just barely inside tolerance" and "dead on" matches, which
    matters once you rank multiple candidate compounds against each other.
    """
    sigma = tolerance_ppm / 2.0 if tolerance_ppm > 0 else 1e-9
    return float(np.exp(-0.5 * (delta_ppm / sigma) ** 2))


def candidate_matches(experimental_ppm: np.ndarray, reference_ppm: np.ndarray,
                       tolerance_ppm: float):
    """
    For every reference peak, find the closest experimental peak within
    `tolerance_ppm`. Returns a list of (ref_index, exp_index, weight) for
    every reference peak that found a hit. A reference peak with no
    experimental peak nearby is simply absent from the output (unmatched).
    """
    out = []
    if len(experimental_ppm) == 0 or len(reference_ppm) == 0:
        return out
    exp = np.asarray(experimental_ppm, dtype=float)
    for r_idx, r_ppm in enumerate(reference_ppm):
        diffs = np.abs(exp - r_ppm)
        e_idx = int(np.argmin(diffs))
        best_delta = float(diffs[e_idx])
        if best_delta <= tolerance_ppm:
            w = gaussian_match_weight(best_delta, tolerance_ppm)
            out.append((r_idx, e_idx, w))
    return out


def score_compound(experimental_ppm, experimental_intensity, reference_ppm,
                    reference_intensity=None, tolerance_ppm: float = 0.008):
    """
    Score one candidate compound against one experimental peak list.

    Returns a dict with:
      n_ref            number of reference peaks for this compound
      n_exp             number of experimental peaks supplied
      n_matched        how many reference peaks found an experimental hit
      n_distinct_matched  how many DIFFERENT experimental peaks those hits
                         actually landed on (<= n_matched, <= n_exp)
      coverage          n_matched / n_ref  (did we see everything this
                         compound should show?)
      purity            n_matched / n_exp  (how much of the experimental
                         spectrum does this compound explain? only meaningful
                         when n_exp is small / peaks are from a near-pure
                         fraction, not a whole-mixture spectrum)
      shift_score       mean Gaussian match weight over matched peaks
                         (rewards close, not just within-tolerance, shifts)
      evidence_score     1 - 1/n_distinct_matched, saturating reward for the
                         ABSOLUTE count of distinct experimental peaks this
                         compound's reference peaks landed on, not just the
                         coverage ratio, and NOT the raw reference-peak hit
                         count either (see below for why that distinction
                         matters). Why this exists: confirmed on a real run
                         against an amino-acid standard mixture, a
                         2-reference-peak compound that matched BOTH of its
                         peaks (coverage = 1.0) was consistently outranking
                         real compounds matched on far more peaks (e.g.
                         L-arginine at 27 of 28, agmatine at 20 of 20) purely
                         because coverage alone can't tell "2 of 2" from
                         "27 of 28" apart -- both read as a perfect ratio,
                         even though 27 independent matching peaks is
                         overwhelmingly stronger evidence against coincidence
                         than 2 is. evidence_score gives 2 distinct matched
                         peaks only 0.5, 4 distinct matched peaks 0.75, 10
                         distinct matched peaks 0.9, 20 distinct matched
                         peaks 0.95 -- real, increasing credit for more
                         corroborating peaks, with diminishing returns so a
                         500-peak library entry can't completely dominate on
                         count alone.

                         Built from n_distinct_matched rather than the raw
                         n_matched count, because candidate_matches() looks
                         up, independently for EACH reference peak, whichever
                         experimental peak is nearest -- so several reference
                         peaks can legitimately all land on the SAME one
                         experimental peak. This is common and not a data
                         error: confirmed on the real E. coli reference
                         library, sugars like maltose and lactose, and other
                         multiplet-heavy compounds, list dozens of
                         closely-spaced sub-peaks (fine multiplet structure,
                         or several literature measurements merged under one
                         compound_id) within a typical matching tolerance of
                         each other (43% of the real E. coli-filtered library
                         has at least one such cluster; Isomaltose alone has
                         45 reference-peak pairs within 0.008 ppm of each
                         other out of 93 total peaks). Against only 3 real
                         experimental peaks sitting inside one such cluster,
                         raw n_matched can come back as 10 or more, even
                         though only 1-3 actual, physically distinct peaks
                         were ever observed -- rewarding evidence_score on
                         that raw count would score reference-library
                         duplication, not independent experimental
                         corroboration, which is exactly the kind of
                         structural bias this term exists to avoid in the
                         first place. n_distinct_matched counts each
                         experimental peak once no matter how many reference
                         peaks happen to land on it, so it can never exceed
                         n_exp.
      peak_prominence    mean percentile rank of the MATCHED peaks' own
                         intensity within this sample's full peak list (0-1,
                         1.0 = matched on the biggest peaks in the whole
                         spectrum, near 0 = matched only on peaks barely
                         above this sample's own noise floor). None if no
                         experimental intensity was supplied. This is a
                         percentile WITHIN one sample's own peak list, not a
                         comparison to other samples or other compounds, so
                         it stays meaningful across spectra with very
                         different absolute intensity scales.
      intensity_corr    Pearson correlation between matched reference and
                         experimental intensities, if intensities were given
                         for both sides; otherwise None
      matched_intensity_mean   mean of THIS sample's own experimental
                         intensities at the matched peaks (not a correlation,
                         an actual value, in whatever units the experimental
                         intensity column was in -- raw peak height or, if
                         you ran refine_peaks_lorentzian() upstream and fed
                         its area column through as "intensity", a real
                         integrated peak area). This is what a qNMR-style
                         fold-change calculation between two groups needs:
                         intensity_corr tells you whether the right RELATIVE
                         pattern of peaks showed up at all (shape check);
                         this tells you how much signal was actually there
                         in this one sample, which is the number you take
                         the ratio of between groups. None if no experimental
                         intensity was supplied.
      confidence        single 0-1 combined score used for ranking, see below

    confidence = 0.30*coverage + 0.25*evidence_score + 0.25*shift_score +
    0.20*peak_prominence (renormalized over just the first three, dividing
    by 0.80, when peak_prominence is unavailable), then folds in
    intensity_corr as a bonus the same way it always did. This replaces the
    older 0.65*coverage + 0.35*shift_score formula, which is why a compound
    with few reference peaks will usually show a LOWER confidence number
    now than it used to, even when it's the correct identification -- that
    compound genuinely only has weak (few-peak) evidence behind it in this
    reference library, and a lower number is the honest reflection of that,
    not a bug. It still won't disappear from your results with --top-n set
    generously high the way this project always runs it; it just won't
    silently outrank a much-better-evidenced compound at the top of the
    list anymore.
    """
    reference_ppm = np.asarray(reference_ppm, dtype=float)
    experimental_ppm = np.asarray(experimental_ppm, dtype=float)
    n_ref = len(reference_ppm)
    n_exp = len(experimental_ppm)

    matches = candidate_matches(experimental_ppm, reference_ppm, tolerance_ppm)
    n_matched = len(matches)
    n_distinct_matched = len({e_idx for _, e_idx, _ in matches})

    coverage = n_matched / n_ref if n_ref else 0.0
    purity = n_matched / n_exp if n_exp else 0.0
    shift_score = float(np.mean([w for _, _, w in matches])) if matches else 0.0
    evidence_score = (1.0 - 1.0 / n_distinct_matched) if n_distinct_matched > 0 else 0.0

    intensity_corr = None
    matched_intensity_mean = None
    peak_prominence = None
    if matches and experimental_intensity is not None:
        experimental_intensity = np.asarray(experimental_intensity, dtype=float)
        if len(experimental_intensity) == n_exp:
            exp_vals_all = [experimental_intensity[e] for _, e, _ in matches]
            if n_exp > 1:
                # Percentile rank of each matched peak's OWN intensity among
                # every peak detected in this one sample (ties count as
                # "at or below", so the single biggest peak in the sample
                # still gets exactly 1.0, not something just under it).
                ranks = [float(np.mean(experimental_intensity <= v)) for v in exp_vals_all]
                peak_prominence = float(np.mean(ranks))
            # Mean, not sum: matched peaks can have different multiplicities
            # (e.g. a CH3 triplet peak vs. a CH singlet peak for the same
            # compound), so summing would mix physically different
            # quantities. Mean of the matched peaks' own intensities is a
            # reasonable automated default for screening many compounds at
            # once; for a compound you care about precisely, go back and
            # pick ONE clean, non-overlapping peak by hand instead of
            # trusting this average.
            matched_intensity_mean = float(np.mean(exp_vals_all))

        if reference_intensity is not None:
            reference_intensity = np.asarray(reference_intensity, dtype=float)
            if len(experimental_intensity) == n_exp and len(reference_intensity) == n_ref:
                ref_vals = [reference_intensity[r] for r, _, _ in matches]
                exp_vals = [experimental_intensity[e] for _, e, _ in matches]
                if len(matches) >= 3 and np.std(ref_vals) > 0 and np.std(exp_vals) > 0:
                    intensity_corr = float(np.corrcoef(ref_vals, exp_vals)[0, 1])

    # Weighted blend: coverage (did we see this compound's own signature),
    # evidence_score (how much independent corroborating peak-count evidence
    # backs that up), and shift quality, each weighted close to evenly, plus
    # peak_prominence (were the matched peaks real signal or near-noise)
    # when intensity data is available. Renormalize over the first three
    # (0.30+0.25+0.25=0.80) when peak_prominence can't be computed, rather
    # than silently treating a missing term as a zero. intensity_corr then
    # folds in as a bonus/penalty exactly as it did before this change.
    if peak_prominence is not None:
        confidence = (0.30 * coverage + 0.25 * evidence_score +
                      0.25 * shift_score + 0.20 * peak_prominence)
    else:
        confidence = (0.30 * coverage + 0.25 * evidence_score + 0.25 * shift_score) / 0.80
    if intensity_corr is not None:
        confidence = 0.85 * confidence + 0.15 * max(0.0, intensity_corr)
    confidence = float(np.clip(confidence, 0.0, 1.0))

    return {
        "n_ref": n_ref,
        "n_exp": n_exp,
        "n_matched": n_matched,
        "n_distinct_matched": n_distinct_matched,
        "coverage": round(coverage, 4),
        "purity": round(purity, 4),
        "shift_score": round(shift_score, 4),
        "evidence_score": round(evidence_score, 4),
        "peak_prominence": round(peak_prominence, 4) if peak_prominence is not None else None,
        "intensity_corr": round(intensity_corr, 4) if intensity_corr is not None else None,
        "matched_intensity_mean": round(matched_intensity_mean, 6) if matched_intensity_mean is not None else None,
        "confidence": round(confidence, 4),
    }


# ---------------------------------------------------------------------------
# OPTIONAL: 2D [1H,13C]-HSQC cross-peak matching. Nothing above this line is
# touched or imported differently for 1D use -- these are purely additive
# and only used by hsqc_matcher.py / hsqc_analysis.py.
# ---------------------------------------------------------------------------

def gaussian_match_weight_2d(delta_h_ppm: float, delta_c_ppm: float,
                              tolerance_h_ppm: float, tolerance_c_ppm: float) -> float:
    """
    2D generalisation of gaussian_match_weight: each dimension's offset is
    first normalised by its OWN tolerance (1H and 13C need very different
    ppm windows -- a typical HSQC has ~100x the spectral width in 13C that
    it does in 1H), then combined as a radius in that normalised space,
    same falloff shape as the 1D case. A hit exactly on both shifts gets
    weight 1.0; a hit sitting right on the edge of the elliptical
    tolerance region gets a small, non-zero weight.
    """
    sigma_h = tolerance_h_ppm / 2.0 if tolerance_h_ppm > 0 else 1e-9
    sigma_c = tolerance_c_ppm / 2.0 if tolerance_c_ppm > 0 else 1e-9
    r2 = (delta_h_ppm / sigma_h) ** 2 + (delta_c_ppm / sigma_c) ** 2
    return float(np.exp(-0.5 * r2))


def candidate_matches_2d(experimental_h: np.ndarray, experimental_c: np.ndarray,
                          reference_h: np.ndarray, reference_c: np.ndarray,
                          tolerance_h_ppm: float, tolerance_c_ppm: float):
    """
    For every reference cross-peak (h_ppm, c_ppm), find the closest
    experimental cross-peak inside the elliptical tolerance region defined
    by tolerance_h_ppm/tolerance_c_ppm (i.e. normalised radius <= 1).
    Returns a list of (ref_index, exp_index, weight); a reference cross-peak
    with nothing nearby in BOTH dimensions at once is absent (unmatched).
    """
    out = []
    if len(experimental_h) == 0 or len(reference_h) == 0:
        return out
    exp_h = np.asarray(experimental_h, dtype=float)
    exp_c = np.asarray(experimental_c, dtype=float)
    for r_idx in range(len(reference_h)):
        r_h, r_c = reference_h[r_idx], reference_c[r_idx]
        dh = exp_h - r_h
        dc = exp_c - r_c
        radius2 = (dh / (tolerance_h_ppm or 1e-9)) ** 2 + (dc / (tolerance_c_ppm or 1e-9)) ** 2
        e_idx = int(np.argmin(radius2))
        if radius2[e_idx] <= 1.0:
            w = gaussian_match_weight_2d(float(dh[e_idx]), float(dc[e_idx]),
                                          tolerance_h_ppm, tolerance_c_ppm)
            out.append((r_idx, e_idx, w))
    return out


def score_compound_2d(experimental_h, experimental_c, reference_h, reference_c,
                       tolerance_h_ppm: float = 0.03, tolerance_c_ppm: float = 0.4):
    """
    2D counterpart of score_compound, for [1H,13C]-HSQC cross-peaks, matched
    on (h_ppm, c_ppm) pairs inside an elliptical tolerance region instead of
    a single ppm axis. Not currently called from anywhere in this pipeline
    (match_r_feature_table.py/peak_matcher.py are 1D-only); kept here ready
    for a future HSQC extension, so its scoring is kept in step with
    score_compound's fix below rather than left on the old formula.

    coverage/shift_score/confidence use the same evidence_score idea as
    score_compound (see that docstring for the motivating real-data case): a
    compound matched on 2 of 2 cross-peaks should not automatically outrank
    one matched on considerably more, just because a ratio alone can't tell
    "2 of 2" apart from "20 of 20". There is no peak_prominence term here,
    since this function takes no intensity values at all (cross-peak
    intensities are less standardised across instruments/processing than 1D
    peak heights, so this was never used as a scoring input for HSQC) --
    confidence here is coverage/evidence_score/shift_score renormalized the
    same way score_compound falls back when its own peak_prominence is
    unavailable. evidence_score is built from n_distinct_matched (distinct
    experimental cross-peaks matched), not the raw reference-cross-peak hit
    count, for the same reason as score_compound's own n_distinct_matched;
    see that function's docstring.
    """
    reference_h = np.asarray(reference_h, dtype=float)
    reference_c = np.asarray(reference_c, dtype=float)
    experimental_h = np.asarray(experimental_h, dtype=float)
    experimental_c = np.asarray(experimental_c, dtype=float)
    n_ref = len(reference_h)
    n_exp = len(experimental_h)

    matches = candidate_matches_2d(experimental_h, experimental_c, reference_h, reference_c,
                                    tolerance_h_ppm, tolerance_c_ppm)
    n_matched = len(matches)
    n_distinct_matched = len({e_idx for _, e_idx, _ in matches})

    coverage = n_matched / n_ref if n_ref else 0.0
    purity = n_matched / n_exp if n_exp else 0.0
    shift_score = float(np.mean([w for _, _, w in matches])) if matches else 0.0
    evidence_score = (1.0 - 1.0 / n_distinct_matched) if n_distinct_matched > 0 else 0.0

    confidence = float(np.clip(
        (0.30 * coverage + 0.25 * evidence_score + 0.25 * shift_score) / 0.80,
        0.0, 1.0
    ))

    return {
        "n_ref": n_ref,
        "n_exp": n_exp,
        "n_matched": n_matched,
        "n_distinct_matched": n_distinct_matched,
        "coverage": round(coverage, 4),
        "purity": round(purity, 4),
        "shift_score": round(shift_score, 4),
        "evidence_score": round(evidence_score, 4),
        "confidence": round(confidence, 4),
    }
