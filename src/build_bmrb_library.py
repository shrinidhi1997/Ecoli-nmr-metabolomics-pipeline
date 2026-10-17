"""
build_bmrb_library.py

Downloads real experimental 1H NMR peak data in bulk from BMRB's
metabolomics collection (bmrb.io/metabolomics) via their documented REST
API (api.bmrb.io), rather than one compound page at a time.

WHAT THIS DOES vs. WHAT IT DOESN'T
  * BMRB's metabolomics collection covers several hundred well-
    characterized small molecules (central metabolism, amino acids,
    nucleotides, sugars, common osmolytes, etc.) with REAL experimental
    spectra, many collected in the same D2O/phosphate buffer/pH 7.4/DSS
    protocol used for real user samples validated against this pipeline.
  * It is NOT all ~3,755 ECMDB compounds. Nobody publishes bulk NMR peak
    data for the long tail of obscure E. coli pathway intermediates.
    What you get here is the realistic ceiling: the compounds actually
    likely to show up as identifiable peaks in a real extract.
  * BMRB's metabolomics collection is NOT organism-specific -- it mixes
    in human, plant, and food-science compounds (e.g. eucalyptol, a plant
    terpene, has shown up in a real bulk pull). Use build_kegg_ecoli_list.py
    to flag which compounds are actually linked to E. coli metabolism.

REAL API SHAPE (confirmed via live --dump-one-entry / --dump-peak-data runs,
not guessed):
  Entries come back as {"bmse00000X": {"entry_id": ..., "saveframes": [...]}}.
  Each saveframe has "category" (lowercase_with_underscores, e.g.
  "entry_information", "assigned_chemical_shifts", "spectral_peak_list"),
  "tags" (a list of [key, value] pairs), and "loops" (each with "category"
  as "_Capitalized_With_Underscore", "tags" as column names, "data" as rows).

  A confirmed real quirk this handles: most compounds have SEVERAL 1D-1H
  spectral_peak_list saveframes (a dilute run, a concentrated run, a plain
  one), and the concentration-series ones (name contains "mM") only carry
  a _Spectral_transition_char loop, not _Peak_char -- extract_1h_peaks()
  prefers the plain experiment's _Peak_char loop, and falls back to
  _Spectral_transition_char (from whichever experiment is available) when
  no plain one exists, rather than silently returning nothing for those.

  A compound with NO spectral_peak_list at all (seen for bmse000001, which
  only has assigned_chemical_shifts) falls back further to that saveframe's
  _Atom_chem_shift loop, filtered to isotope-1 hydrogen only.

  A confirmed real "hidden name" quirk: entry_information's _Entry.Title
  is usually a clean name (e.g. "L_leucine" for bmse000042), but for some
  sub-collections it is only an internal catalog code -- a real example a
  user pasted directly from BMRB's native text format, bmse010280, has
  Title "lignin_cw_compound_1026" (a USDA lignin/cell-wall reference
  compound with no chemical meaning in that string at all). The real name
  lives in the chem_comp saveframe's _Chem_comp_common_name loop instead.
  extract_common_names() reads that loop; resolve_compound_name() decides
  which of the two to use for display (see its docstring for why this is
  not simply "always prefer the synonym loop").
"""

from __future__ import annotations
import argparse
import json
import re
import time
from pathlib import Path

import pandas as pd
import requests

API_ROOT = "https://api.bmrb.io/v2"
HEADERS = {"Application": "ecoli-nmr-pipeline (contact: set-your-email-here)"}


def get_json(url: str, params: dict | None = None, timeout: int = 30):
    resp = requests.get(url, params=params, headers=HEADERS, timeout=timeout)
    resp.raise_for_status()
    return resp.json()


def list_metabolomics_entries() -> list[str]:
    data = get_json(f"{API_ROOT}/list_entries", params={"database": "metabolomics"})
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        for key in ("metabolomics", "data", "entries", "ids"):
            if key in data and isinstance(data[key], list):
                return data[key]
    raise RuntimeError(
        f"Unrecognised list_entries response shape: {str(data)[:300]}. "
        "Open a browser to https://api.bmrb.io/v2/list_entries?database=metabolomics "
        "to see the real shape and adjust list_metabolomics_entries()."
    )


def _tags_to_dict(tags: list) -> dict:
    """A saveframe's 'tags' is a list of [key, value] pairs -> dict."""
    return {k: v for k, v in tags}


def _loop_to_records(loop: dict) -> list[dict]:
    """A loop's 'tags' (column names) + 'data' (rows) -> list of row-dicts."""
    cols = loop.get("tags", [])
    return [dict(zip(cols, row)) for row in loop.get("data", [])]


def _find_loops(saveframe: dict, category_candidates: list[str]) -> list[dict]:
    return [loop for loop in saveframe.get("loops", [])
            if loop.get("category", "") in category_candidates]


_NAME_NOISE_RE = re.compile(r"[^a-z0-9]")


def _normalize_name_for_comparison(name: str) -> str:
    """Lowercase and strip everything but letters/digits, so 'L_leucine',
    'L-Leucine', and 'L Leucine' all compare equal regardless of
    punctuation/spacing differences."""
    return _NAME_NOISE_RE.sub("", name.lower())


def extract_common_names(entry_json: dict, entry_id: str) -> list[str]:
    """
    Real rows from the chem_comp saveframe's _Chem_comp_common_name loop
    (columns confirmed via a live fetch of bmse000042/L-Leucine: "Name",
    "Type", "Entry_ID", "Comp_ID" -- Type is "synonym" for every row seen).
    Placeholder values ("?", "to be defined", etc., confirmed real on the
    neighbouring _Chem_comp_systematic_name field for some entries) are
    dropped.

    Confirmed real shape: bmse000042/L-Leucine has 22 synonym rows (common
    compounds are heavily cross-referenced). bmse010280 -- a USDA
    lignin/cell-wall model compound the user pasted directly from BMRB's
    native NMR-STAR text -- has exactly ONE row, and that row holds the
    ONLY real descriptive name available anywhere in the entry: this
    compound's _Entry.Title and _Chem_comp.Name are both just the internal
    catalog code "lignin_cw_compound_1026".
    """
    entry = entry_json.get(entry_id, entry_json)
    names: list[str] = []
    for sf in entry.get("saveframes", []):
        if sf.get("category") != "chem_comp":
            continue
        for loop in _find_loops(sf, ["_Chem_comp_common_name"]):
            for rec in _loop_to_records(loop):
                syn = str(rec.get("Name", "")).strip()
                if syn and syn.lower() not in ("?", "n/a", "none", "to be defined"):
                    names.append(syn)
    return names


def resolve_compound_name(title: str, common_names: list[str], entry_id: str) -> str:
    """
    Picks the best display name for a compound, given the raw
    entry_information _Entry.Title (what extract_1h_peaks()/
    extract_hsqc_peaks() used alone, before this fix) and the chem_comp
    saveframe's common-name synonym list (see extract_common_names(),
    never read before this fix).

    CONFIRMED REAL BUG this fixes: the user pasted a real raw BMRB entry
    (bmse010280) whose Title is just an internal catalog code,
    "lignin_cw_compound_1026" -- not a chemical name at all -- while the
    real descriptive name ("3-(4-acetoxy-3-methoxyphenyl) acrylic acid
    4-[1,3-diacetoxy-2-(2-methoxyphenoxy) propyl] phenyl ester") sat
    unread in the common-name loop the whole time.

    The fix is NOT "always prefer the common-name loop" -- for the
    overwhelmingly common case (a well-known compound like L-Leucine)
    Title is already clean ("L_leucine"), and blindly taking the synonym
    loop's FIRST row would make things WORSE: bmse000042's real first row
    is "L-Norvaline, 4-methyl-", a genuine but obscure synonym, not the
    clean name Title already gave. Instead: trust Title whenever it
    matches (case/punctuation stripped) something already in the synonym
    list -- i.e. whenever it IS a real recognized name. Only when Title
    matches none of the synonyms -- the actual signal that it's an
    internal placeholder, not a chemical name -- fall back to the synonym
    list (first entry, for lack of a way to rank them; for bmse010280
    there is only one entry, which is exactly the real name needed).
    """
    title = (title or "").strip()
    if title and common_names:
        norm_title = _normalize_name_for_comparison(title)
        if any(_normalize_name_for_comparison(s) == norm_title for s in common_names):
            return title
        return common_names[0]
    if title:
        return title
    if common_names:
        return common_names[0]
    return entry_id


def extract_1h_peaks(entry_json: dict, entry_id: str) -> tuple[str, list[float], str]:
    """
    Returns (compound_name, list_of_1H_shifts, inchikey_or_empty).
    Prefers real observed peak positions (a spectral_peak_list saveframe's
    _Peak_char loop) over assigned-atom shifts, since peak positions are
    what a picked experimental spectrum actually shows; falls back to
    assigned_chemical_shifts's _Atom_chem_shift loop when no peak list is
    present for this entry (common for very simple/symmetric compounds).
    """
    entry = entry_json.get(entry_id, entry_json)  # API nests the entry under its own ID
    saveframes = entry.get("saveframes", [])

    raw_title = ""
    peaks: list[float] = []
    inchikey = ""
    candidate_peak_saveframes: list[dict] = []

    for sf in saveframes:
        category = sf.get("category", "")
        tags = _tags_to_dict(sf.get("tags", []))

        if category == "entry_information":
            raw_title = tags.get("Title") or raw_title

        if category in ("chem_comp", "entity"):
            inchikey = (tags.get("Standard_InChI_key") or tags.get("InChI_key")
                        or tags.get("InChI_code") or inchikey)

        if category == "spectral_peak_list":
            dim_loops = _find_loops(sf, ["_Spectral_dim"])
            dims = _loop_to_records(dim_loops[0]) if dim_loops else []
            is_1d_1h = (len(dims) == 1
                        and str(dims[0].get("Atom_type", "")).upper() == "H")
            if is_1d_1h:
                candidate_peak_saveframes.append(sf)

    # Real entries commonly have SEVERAL 1D-1H peak lists (a dilute run, a
    # concentrated run, a plain one), confirmed from a live --dump-peak-data
    # run. Concentration-series variants (name contains "mM") only carry a
    # _Spectral_transition_char loop, no _Peak_char; the plain experiment
    # (no concentration in its name) is the one with the curated _Peak_char
    # loop. Pick ONE experiment to use rather than merging all of them --
    # merging would pile up near-duplicate, slightly inconsistent shift
    # values under one compound instead of one clean set. Prefer the plain
    # experiment; fall back to whichever candidate exists if there isn't one.
    chosen_sf = None
    for sf in candidate_peak_saveframes:
        exp_name = str(_tags_to_dict(sf.get("tags", [])).get("Experiment_name", ""))
        if "mM" not in exp_name and "mg" not in exp_name.lower():
            chosen_sf = sf
            break
    if chosen_sf is None and candidate_peak_saveframes:
        chosen_sf = candidate_peak_saveframes[0]

    if chosen_sf is not None:
        for loop in _find_loops(chosen_sf, ["_Peak_char"]):
            for rec in _loop_to_records(loop):
                val = rec.get("Chem_shift_val")
                if val is not None:
                    try:
                        peaks.append(float(val))
                    except (TypeError, ValueError):
                        pass
        if not peaks:  # this experiment has no _Peak_char -- use raw transitions
            for loop in _find_loops(chosen_sf, ["_Spectral_transition_char"]):
                for rec in _loop_to_records(loop):
                    val = rec.get("Chem_shift_val")
                    if val is not None:
                        try:
                            peaks.append(float(val))
                        except (TypeError, ValueError):
                            pass

    if not peaks:
        for sf in saveframes:
            if sf.get("category") != "assigned_chemical_shifts":
                continue
            shift_loops = _find_loops(sf, ["_Atom_chem_shift"])
            for loop in shift_loops:
                for rec in _loop_to_records(loop):
                    atom_type = str(rec.get("Atom_type", "")).upper()
                    isotope = str(rec.get("Atom_isotope_number", ""))
                    if atom_type == "H" and isotope in ("1", ""):
                        val = rec.get("Val")
                        if val is not None:
                            try:
                                peaks.append(float(val))
                            except (TypeError, ValueError):
                                pass

    peaks = sorted(set(round(p, 4) for p in peaks))
    name = resolve_compound_name(raw_title, extract_common_names(entry_json, entry_id), entry_id)
    return name, peaks, inchikey


def extract_hsqc_peaks(entry_json: dict, entry_id: str) -> tuple[str, list[tuple[float, float]], str]:
    """
    OPTIONAL, HSQC-specific. Returns (compound_name, list_of_(h_ppm, c_ppm)
    cross-peaks, inchikey_or_empty) from a 2D [1H,13C] spectral_peak_list
    saveframe. Building the 1H-only reference library (extract_1h_peaks,
    above) does not call this and is completely unaffected by it.

    Confirmed real quirk (via a live --dump-peak-data run on a 2D entry):
    a 2D _Peak_char loop stores each cross-peak SPLIT across TWO rows, one
    per Spectral_dim_ID, not one row per peak:
        _Spectral_dim:  ID=1 Atom_type=H ...   ID=2 Atom_type=C ...
        _Peak_char:     Peak_ID=1 Spectral_dim_ID=1 Chem_shift_val=3.104
                        Peak_ID=1 Spectral_dim_ID=2 Chem_shift_val=39.337
    This groups _Peak_char rows by Peak_ID, uses _Spectral_dim to learn
    which Spectral_dim_ID is the proton and which is the carbon for THIS
    saveframe (dimension order isn't assumed to be fixed across entries),
    and pairs them back into one (h_ppm, c_ppm) cross-peak per Peak_ID. A
    Peak_ID missing either half (incomplete assignment) is dropped rather
    than guessed at.
    """
    entry = entry_json.get(entry_id, entry_json)
    saveframes = entry.get("saveframes", [])

    raw_title = ""
    inchikey = ""
    hsqc_candidates: list[tuple[dict, dict]] = []

    for sf in saveframes:
        category = sf.get("category", "")
        tags = _tags_to_dict(sf.get("tags", []))

        if category == "entry_information":
            raw_title = tags.get("Title") or raw_title
        if category in ("chem_comp", "entity"):
            inchikey = (tags.get("Standard_InChI_key") or tags.get("InChI_key")
                        or tags.get("InChI_code") or inchikey)

        if category == "spectral_peak_list":
            n_dims = str(tags.get("Number_of_spectral_dimensions", ""))
            dim_loops = _find_loops(sf, ["_Spectral_dim"])
            if not dim_loops:
                continue
            dims = _loop_to_records(dim_loops[0])
            atom_types = sorted(str(d.get("Atom_type", "")).upper() for d in dims)
            # Any 2D H+C peak list counts (HSQC, HMQC, HSQC-TOCSY, ...) --
            # what matters for matching is that it's a real 2D H-C
            # correlation, not the exact pulse-sequence name.
            if n_dims == "2" and atom_types == ["C", "H"]:
                hsqc_candidates.append((sf, dim_loops[0]))

    name = resolve_compound_name(raw_title, extract_common_names(entry_json, entry_id), entry_id)

    if not hsqc_candidates:
        return name, [], inchikey

    chosen_sf, chosen_dim_loop = hsqc_candidates[0]
    for sf, dim_loop in hsqc_candidates:
        exp_name = str(_tags_to_dict(sf.get("tags", [])).get("Experiment_name", ""))
        if "mM" not in exp_name and "mg" not in exp_name.lower():
            chosen_sf, chosen_dim_loop = sf, dim_loop
            break

    dim_records = _loop_to_records(chosen_dim_loop)
    dim_to_atom = {str(d.get("ID")): str(d.get("Atom_type", "")).upper() for d in dim_records}

    by_peak: dict[str, dict[str, float]] = {}
    for loop in _find_loops(chosen_sf, ["_Peak_char"]):
        for rec in _loop_to_records(loop):
            peak_id = str(rec.get("Peak_ID"))
            atom = dim_to_atom.get(str(rec.get("Spectral_dim_ID")))
            val = rec.get("Chem_shift_val")
            if atom not in ("H", "C") or val is None:
                continue
            try:
                by_peak.setdefault(peak_id, {})[atom] = float(val)
            except (TypeError, ValueError):
                continue

    cross_peaks = sorted(
        (v["H"], v["C"]) for v in by_peak.values() if "H" in v and "C" in v
    )
    return name, cross_peaks, inchikey


def extract_db_ids(entry_json: dict, entry_id: str) -> dict[str, str]:
    """
    Returns {'kegg_id':..., 'chebi_id':..., 'pubchem_id':...}, each '' if
    not cross-referenced by this BMRB entry's chem_comp saveframe.

    CONFIRMED REAL SCHEMA (live fetch of bmse000042 / L-Leucine's
    '_Chem_comp_db_link' loop, columns: Author_supplied, Database_code,
    Accession_code, Accession_code_type, Entry_mol_name,
    Entry_relation_type, Entry_ID, Comp_ID):

      KEGG:     Database_code == "KEGG" (exact). Accession_code is the
                compound ID directly, e.g. "C00123".
      ChEBI:    Database_code == "CHEBI" -- ALL CAPS in the real data, not
                "ChEBI" as written everywhere else including ChEBI's own
                branding. Accession_code is the bare numeric ID, e.g.
                "15603" (matches ECMDB's chebi_id field format directly,
                e.g. "16763" -- no "CHEBI:" prefix on either side).
      PubChem:  Database_code == "PubChem". BMRB lists MULTIPLE PubChem
                rows per compound mixing two different ID namespaces via
                Accession_code_type: "sid" (substance ID) and "cid"
                (compound ID). bmse000042 has three "sid" rows (85164890,
                148996, 3423) and exactly one "cid" row (6106). ECMDB's own
                pubchem_id field (e.g. "58") is a compound ID too, so only
                Accession_code_type == "cid" rows are usable here --
                comparing a "sid" row against ECMDB's CID would either
                never match (a silent coverage loss, since a real match
                exists but under the cid row this code would have ignored)
                or, worse, coincidentally collide with an unrelated
                compound's ID in the other namespace. Checking
                Accession_code_type, not just Database_code, is what
                prevents both failure modes.

    Replaces the narrower extract_kegg_id() (kept below as a thin wrapper
    for compatibility) now that ChEBI/PubChem give a second and third exact
    cross-reference path into ECMDB alongside KEGG -- see
    load_ecmdb_crossref()/is_ecoli_by_ecmdb() for where these get used.
    """
    entry = entry_json.get(entry_id, entry_json)
    ids = {"kegg_id": "", "chebi_id": "", "pubchem_id": ""}
    for sf in entry.get("saveframes", []):
        if sf.get("category") != "chem_comp":
            continue
        for loop in _find_loops(sf, ["_Chem_comp_db_link"]):
            for rec in _loop_to_records(loop):
                code = str(rec.get("Database_code", "")).strip().upper()
                acc = str(rec.get("Accession_code", "")).strip()
                acc_type = str(rec.get("Accession_code_type", "")).strip().lower()
                if not acc:
                    continue
                if code == "KEGG" and not ids["kegg_id"]:
                    ids["kegg_id"] = acc
                elif code == "CHEBI" and not ids["chebi_id"]:
                    ids["chebi_id"] = acc
                elif code == "PUBCHEM" and acc_type == "cid" and not ids["pubchem_id"]:
                    ids["pubchem_id"] = acc
    return ids


def extract_kegg_id(entry_json: dict, entry_id: str) -> str:
    """Thin wrapper over extract_db_ids(), kept for compatibility with
    existing callers/tests that only want the KEGG ID. See extract_db_ids()
    for the confirmed real schema this reads."""
    return extract_db_ids(entry_json, entry_id)["kegg_id"]


def load_ecmdb_crossref(ecmdb_library_csv: str | None) -> dict[str, set[str]]:
    """
    Loads compound_library.csv (build_library.py, from a real downloaded
    ecmdb.json) into a dict of ID sets: kegg_id, chebi_id, pubchem_id,
    inchikey. Because ECMDB is BY DEFINITION an E. coli-specific compound
    list (curated for E. coli K-12 MG1655), any BMRB compound whose own id
    in any of these fields appears in the matching ECMDB set is
    authoritatively confirmed to belong in an E. coli library -- this is
    what is_ecoli_by_ecmdb() checks.

    kegg_id, chebi_id and pubchem_id are all exact-ID matches, each
    confirmed as a real, directly-comparable field on both sides via a
    live BMRB fetch (see extract_db_ids()'s docstring for the exact
    Database_code/Accession_code_type quirks that make this safe,
    especially PubChem's sid-vs-cid distinction). Three independent exact
    IDs instead of one catches BMRB compounds that have (say) a ChEBI link
    but no KEGG one, without any name-matching or fuzzy comparison.

    inchikey is checked too but is weaker in practice: a live check of a
    real BMRB entry found NO genuine InChIKey field at all (no
    Standard_InChI_key/InChI_key tag present) -- only a full InChI STRING
    under InChI_code, which extract_1h_peaks()/extract_hsqc_peaks() have
    historically stored under the "inchikey" name for lack of a better
    field. That string will not equal ECMDB's real InChIKey hash, so this
    path only helps on the (apparently rare) BMRB entry that does carry a
    genuine Standard_InChI_key/InChI_key tag -- kept as a free, harmless
    extra check, not a primary one.
    """
    empty = {"kegg_id": set(), "chebi_id": set(), "pubchem_id": set(), "inchikey": set()}
    if not ecmdb_library_csv:
        return empty
    df = pd.read_csv(ecmdb_library_csv)
    out = dict(empty)
    for col in out:
        if col in df.columns:
            out[col] = set(df[col].dropna().astype(str))
        else:
            print(f"WARNING: --ecmdb-library file has no '{col}' column, "
                  f"{col}-based is_ecoli_metabolite matching will be skipped.")
    return out


def is_ecoli_by_ecmdb(ids: dict[str, str],
                       ecmdb_crossref: dict[str, set[str]]) -> bool | None:
    """
    True if any of this compound's own kegg_id/chebi_id/pubchem_id/inchikey
    (ids, as returned by extract_db_ids() plus an 'inchikey' key) is
    confirmed present in the matching ECMDB set (see load_ecmdb_crossref);
    False if a crossref was loaded and none matched; None if no
    --ecmdb-library was given at all (the caller should treat that as "not
    checked", not "confirmed absent").
    """
    if not any(ecmdb_crossref.values()):
        return None
    for field in ("kegg_id", "chebi_id", "pubchem_id", "inchikey"):
        value = ids.get(field, "")
        if value and value in ecmdb_crossref.get(field, ()):
            return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="../output/reference_peaks_bmrb_library.csv")
    ap.add_argument("--cache-dir", default="../output/bmrb_cache",
                     help="per-entry JSON cached here so an interrupted run can resume")
    ap.add_argument("--limit", type=int, default=None,
                     help="only fetch the first N entries -- use this for a first test run")
    ap.add_argument("--delay", type=float, default=0.3,
                     help="seconds to wait between requests (politeness / avoids rate limiting)")
    ap.add_argument("--dump-one-entry", action="store_true",
                     help="print the first entry's raw JSON and exit, for verifying key names")
    ap.add_argument("--dump-peak-data", action="store_true",
                     help="scan entries until one has real peak/shift data, print just that "
                          "saveframe's loops (tags + a few rows), and exit -- more useful than "
                          "--dump-one-entry for checking extract_1h_peaks()'s field names")
    ap.add_argument("--ecmdb-library", default=None,
                     help="path to compound_library.csv from build_library.py (built "
                          "from a real ecmdb.json), for is_ecoli_metabolite flagging "
                          "against ECMDB's authoritative E. coli compound list -- "
                          "matched by exact KEGG Compound ID, ChEBI ID, or PubChem "
                          "Compound ID (any one matching is enough), with InChIKey as "
                          "a weaker secondary check (optional; see load_ecmdb_crossref())")
    ap.add_argument("--hsqc", action="store_true",
                     help="OPTIONAL. Build a 2D [1H,13C]-HSQC reference library "
                          "(schema: compound_id, compound_name, h_ppm, c_ppm, "
                          "inchikey, is_ecoli_metabolite) instead of the default "
                          "1H-only one. Leave this off for ordinary 1D use -- it "
                          "changes nothing about the default behaviour. If --out "
                          "is left at its default, it is switched to "
                          "reference_peaks_hsqc_bmrb_library.csv automatically.")
    args = ap.parse_args()

    if args.hsqc and args.out == ap.get_default("out"):
        args.out = "../output/reference_peaks_hsqc_bmrb_library.csv"

    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    print("Fetching entry list...")
    entry_ids = list_metabolomics_entries()
    print(f"Found {len(entry_ids)} metabolomics entries.")
    if args.limit:
        entry_ids = entry_ids[: args.limit]

    if args.dump_one_entry:
        first = get_json(f"{API_ROOT}/entry/{entry_ids[0]}")
        print(json.dumps(first, indent=2)[:4000])
        return

    if args.dump_peak_data:
        for entry_id in entry_ids[: max(args.limit or 20, 20)]:
            entry_json = get_json(f"{API_ROOT}/entry/{entry_id}")
            entry = entry_json.get(entry_id, entry_json)
            found_any = False
            for sf in entry.get("saveframes", []):
                if sf.get("category") in ("spectral_peak_list", "assigned_chemical_shifts"):
                    found_any = True
                    print(f"=== {entry_id}, saveframe category={sf.get('category')} ===")
                    print("tags:", _tags_to_dict(sf.get("tags", [])))
                    for loop in sf.get("loops", []):
                        print(f"  loop category={loop.get('category')}")
                        print(f"  columns: {loop.get('tags')}")
                        print(f"  first 3 rows: {loop.get('data', [])[:3]}")
            if found_any:
                return
        print("No entry in the first batch had a spectral_peak_list or "
              "assigned_chemical_shifts saveframe -- try a larger --limit.")
        return

    ecmdb_crossref = load_ecmdb_crossref(args.ecmdb_library)

    rows = []
    for i, entry_id in enumerate(entry_ids):
        cache_file = cache_dir / f"{entry_id}.json"
        try:
            if cache_file.exists():
                entry_json = json.loads(cache_file.read_text())
            else:
                entry_json = get_json(f"{API_ROOT}/entry/{entry_id}")
                cache_file.write_text(json.dumps(entry_json))
                time.sleep(args.delay)
        except requests.RequestException as e:
            print(f"  [{i+1}/{len(entry_ids)}] {entry_id}: request failed ({e}), skipping")
            continue

        db_ids = extract_db_ids(entry_json, entry_id)
        kegg_id = db_ids["kegg_id"]

        if args.hsqc:
            try:
                name, cross_peaks, inchikey = extract_hsqc_peaks(entry_json, entry_id)
            except Exception as e:
                print(f"  [{i+1}/{len(entry_ids)}] {entry_id}: parse failed ({e}), skipping")
                continue
            if not cross_peaks:
                print(f"  [{i+1}/{len(entry_ids)}] {entry_id} ({name}): "
                      "no HSQC cross-peaks found, skipping")
                continue
            is_ecoli = is_ecoli_by_ecmdb({**db_ids, "inchikey": inchikey}, ecmdb_crossref)
            for h_ppm, c_ppm in cross_peaks:
                rows.append({
                    "compound_id": entry_id, "compound_name": name,
                    "h_ppm": h_ppm, "c_ppm": c_ppm,
                    "inchikey": inchikey, "kegg_id": kegg_id,
                    "chebi_id": db_ids["chebi_id"], "pubchem_id": db_ids["pubchem_id"],
                    "is_ecoli_metabolite": is_ecoli,
                })
            print(f"  [{i+1}/{len(entry_ids)}] {entry_id} ({name}): "
                  f"{len(cross_peaks)} HSQC cross-peaks")
            continue

        try:
            name, peaks, inchikey = extract_1h_peaks(entry_json, entry_id)
        except Exception as e:
            print(f"  [{i+1}/{len(entry_ids)}] {entry_id}: parse failed ({e}), skipping")
            continue

        if not peaks:
            print(f"  [{i+1}/{len(entry_ids)}] {entry_id} ({name}): no 1H peaks found, skipping")
            continue

        is_ecoli = is_ecoli_by_ecmdb({**db_ids, "inchikey": inchikey}, ecmdb_crossref)
        for shift in peaks:
            rows.append({
                "compound_id": entry_id, "compound_name": name, "nucleus": "1H",
                "shift_ppm": shift, "inchikey": inchikey, "kegg_id": kegg_id,
                "chebi_id": db_ids["chebi_id"], "pubchem_id": db_ids["pubchem_id"],
                "is_ecoli_metabolite": is_ecoli,
            })
        print(f"  [{i+1}/{len(entry_ids)}] {entry_id} ({name}): {len(peaks)} peaks"
              + (f", KEGG={kegg_id}" if kegg_id else ""))

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    n_compounds = df["compound_id"].nunique() if not df.empty else 0
    kind = "HSQC cross-peaks" if args.hsqc else "peaks"
    print(f"\nWrote {len(df)} {kind} across {n_compounds} compounds to {args.out}")
    if any(ecmdb_crossref.values()) and not df.empty:
        confirmed = df[df["is_ecoli_metabolite"] == True]["compound_id"].nunique()
        print(f"{confirmed} of those compounds cross-referenced to a confirmed ECMDB entry "
              f"(by KEGG ID, ChEBI ID, PubChem ID, and/or InChIKey).")


if __name__ == "__main__":
    main()
