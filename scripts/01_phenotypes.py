"""Build the master participants table (one row per subject, all sites) + data dictionary.

Outputs (in --out, default <rawdata>/):
  participants.tsv    master table (BIDS-style; the single source of truth for phenotypes)
  participants.json   data dictionary
  code/missingness_by_site.tsv
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from tinception import config
from tinception.phenotypes import harmonize

DICTIONARY = {
    "participant_id": "Unique BIDS label sub-<site><original id>.",
    "site": "Recruiting site key (config/sites.yaml).",
    "original_id": "Subject ID as used by the site / legacy files.",
    "legacy_id": "Tinception v1 ID (sub-SSNNN) when available.",
    "group": "CO = no tinnitus, TI = chronic tinnitus.",
    "group4": "CTRL (normal hearing, no tinnitus), HLOS (hearing loss, no tinnitus), TINN (tinnitus, normal hearing), TNHL (tinnitus + hearing loss).",
    "group4_source": "'recruitment' when defined by the site's design, 'pta4' when derived from PTA4 > 25 dB HL.",
    "hearing_aid": "Hearing-aid user (cogtail only).",
    "sex": "F / M. Legacy 0/1 coding mapped 0=F, 1=M (verified by eTIV at all sites).",
    "age": {"Description": "Age at scan.", "Units": "years"},
    "thr_<ear>_<freq>": {"Description": "Pure-tone threshold; ear R/L or B (binaural value reported by the site or mean of R and L).", "Units": "dB HL"},
    "PTA_reported": {"Description": "PTA as provided by the site (definition varies by site; not used when an audiogram exists).", "Units": "dB HL"},
    "PTA4": {"Description": "Binaural mean threshold at 0.5/1/2/4 kHz; falls back to PTA_reported where no audiogram exists (see pta_source).", "Units": "dB HL"},
    "PTA_HF": {"Description": "Binaural mean threshold at 4 and 8 kHz (primary high-frequency measure).", "Units": "dB HL"},
    "PTA_HF468": {"Description": "Binaural mean threshold at 4/6/8 kHz (subset of sites).", "Units": "dB HL"},
    "PTA_EHF": {"Description": "Binaural mean over available 10–16 kHz thresholds (>=2 required).", "Units": "dB HL"},
    "PTA4_better_ear": {"Description": "min(PTA4 right, PTA4 left).", "Units": "dB HL"},
    "PTA4_asymmetry": {"Description": "|PTA4 right - PTA4 left|.", "Units": "dB"},
    "HL": "PTA4 > 25 dB HL.",
    "pta_source": "computed_pta4 | reported.",
    "THI": "Tinnitus Handicap Inventory (0–100).",
    "TFI": "Tinnitus Functional Index (0–100).",
    "BDI": "Beck Depression Inventory.",
    "BAI": "Beck Anxiety Inventory.",
    "TPFQ": "Tinnitus Primary Function Questionnaire.",
    "duration_months": {"Description": "Tinnitus duration.", "Units": "months"},
    "distress_instrument": "Instrument used for distress_raw (THI preferred over TFI).",
    "distress_raw": "THI or TFI total score.",
    "distress_z": "distress_raw z-scored within instrument across patients.",
    "affect_z": "Mean within-site z-score of available depression/anxiety scales.",
    "tinnitus": "1 = TI, 0 = CO.",
    "phenotype_status": "ok | no_image (phenotype row without T1) | no_phenotype (image without phenotype row).",
    "legacy_drop": "Excluded in Tinception v1 (DROPS folder, reason unrecorded); re-evaluated by 03_qc.",
    "fs_dir": "Existing FreeSurfer 8.0.0 subject directory (relative to the FreeSurfer root), if any.",
    "has_T1w": "T1w present in rawdata.", "has_T2w": "Any T2w present.",
    "has_dwi": "DWI present.", "has_rest": "Resting-state fMRI present.",
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=None, help="output directory (default: rawdata)")
    ap.add_argument("--sites", nargs="*", default=None, help="subset of site keys")
    args = ap.parse_args()

    out = args.out or config.paths().rawdata
    (out / "code").mkdir(parents=True, exist_ok=True)

    idmap_file = out / "code" / "id_map.tsv"
    id_map = pd.read_csv(idmap_file, sep="\t", na_values="n/a") if idmap_file.exists() else None
    if id_map is None:
        print("note: code/id_map.tsv not found (run 02_bids first) -> no image availability columns")
    df = harmonize.build_master(args.sites, id_map)
    df.to_csv(out / "participants.tsv", sep="\t", index=False, na_rep="n/a")
    with open(out / "participants.json", "w") as f:
        json.dump(DICTIONARY, f, indent=2)
    miss = harmonize.missingness_report(df)
    miss.to_csv(out / "code" / "missingness_by_site.tsv", sep="\t")

    print(f"wrote {len(df)} participants from {df.site.nunique()} sites -> {out}")
    print(df.groupby(["site", "group"]).size().unstack(fill_value=0).to_string())
    print(miss.to_string())


if __name__ == "__main__":
    main()
