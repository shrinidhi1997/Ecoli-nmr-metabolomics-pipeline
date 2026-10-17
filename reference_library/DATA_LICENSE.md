# License for the files in this folder

This reference library was compiled by Shrinidhi Bhat, DMST, Indian
Institute of Technology Madras.

The `LICENSE` file at the root of this repository (MIT) covers the code in
`src/`. It does not cover `reference_peaks_combined.csv` and
`reference_peaks_combined_ecoli_only.csv` in this folder. Those two files
are compiled from two outside databases with their own, different terms,
checked on 2026-10-17:

**HMDB** (the Human Metabolome Database, hmdb.ca). Per
[hmdb.ca/about](https://www.hmdb.ca/about), HMDB's data is released under
Creative Commons Attribution-NonCommercial (CC BY-NC) 4.0. That license
permits non-commercial use and redistribution with attribution, but
commercial use or redistribution requires permission from HMDB directly
(contact given on that page). `reference_peaks_combined.csv` contains
24,659 of its 60,913 rows sourced from HMDB, and the HMDB share is larger
in the *E. coli*-filtered file: 5,772 of 8,607 rows. Cite the current HMDB
paper (listed on the About page) if you use data from these files.

**BMRB** (the Biological Magnetic Resonance Data Bank, bmrb.io). BMRB's
[Aims and Policies](https://bmrb.io/bmrb/aims_and_policies.shtml) and
[metabolomics data policy](https://bmrb.io/metabolomics/data_policy.shtml)
pages say the data is "made publicly available free of charge" and ask
that users "acknowledge BMRB in publications." Neither page states an
explicit license or explicit redistribution terms one way or the other.
Treat BMRB's share of this data as available for attributed academic
reuse, the way it is commonly used in the field, rather than as formally
public-domain or formally CC-licensed, since BMRB itself hasn't said
either.

**What this means in practice.** Because HMDB's CC BY-NC 4.0 term applies
to part of this data, these two CSVs as a whole are available for
non-commercial, attributed use and redistribution, not the unrestricted
reuse MIT grants for the code. If you, a reader, or a reviewer wants to
use this data commercially, the HMDB-sourced rows need permission from
HMDB directly; the rest of this repository's code has no such
restriction.

This is a summary for practical use, written by an AI assistant from the
two databases' own public pages, not a legal opinion. If anything here
matters for your institution's policy or your journal's data-availability
requirements, check with your institution's legal or technology-transfer
office before publishing, particularly given the mixed licensing terms.
