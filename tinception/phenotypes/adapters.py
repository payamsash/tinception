"""Per-source phenotype adapters.

Each adapter returns one row per subject in the common schema below. Values are harmonised
(sex as F/M, thresholds as thr_<ear>_<freq>). Derived scores are computed later, in harmonize.py.
New sites only need a new adapter plus an entry in config/sites.yaml.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from .. import config, ids

SCHEMA = [
    "participant_id", "site", "original_id", "group", "group4", "hearing_aid", "sex", "age",
    "THI", "TFI", "BDI", "BAI", "TPFQ", "duration_months", "PTA_reported", "legacy_id",
    "HADS_A", "HADS_D", "HADS", "STAI", "ISI", "education_years", "hyperacusis_khalfa",
]

_EAR_COL = re.compile(r"^(RH|LH)-(\d+)$")
# cogtail stores a single binaural threshold under 'RH-*' (see legacy docs); treat as binaural.
BINAURAL_ONLY_SITES = {"cogtail"}
# cogtail recruits only hearing-impaired people: NT = no tinnitus, TI = tinnitus; HL/HA = hearing
# loss without / with hearing aid. There are no normal-hearing controls at this site.
_COGTAIL_SUFFIX = re.compile(r"_+(NT|TI)-(HL|HA)$")


def _sex(v) -> object:
    coding = config.analysis()["phenotypes"]["sex_coding"]
    if isinstance(v, str):
        v = v.strip().upper()
        return v if v in ("F", "M") else np.nan
    if pd.isna(v):
        return np.nan
    return coding.get(int(v), np.nan)


def _finalise(df: pd.DataFrame) -> pd.DataFrame:
    for c in SCHEMA:
        if c not in df:
            df[c] = np.nan
    thr = sorted(c for c in df.columns if c.startswith("thr_"))
    return df[SCHEMA + thr]


def legacy_master(site: str) -> pd.DataFrame:
    """Per-site master_files/<site>_master.xlsx (the curated legacy phenotype sheets)."""
    p = config.paths()
    raw = pd.read_excel(p.legacy_masters / f"{site}_master.xlsx")
    df = pd.DataFrame(index=raw.index)
    df["original_id"] = [ids.canonical_original_id(site, s) for s in raw["subject_ID"]]
    df["site"] = site
    df["participant_id"] = [ids.bids_id(site, s) for s in df["original_id"]]
    df["group"] = raw["group"].str.strip().str.upper()
    df["sex"] = raw["sex"].map(_sex)
    df["age"] = pd.to_numeric(raw["age"], errors="coerce")
    df["THI"] = pd.to_numeric(raw.get("THI"), errors="coerce")
    df["PTA_reported"] = pd.to_numeric(raw.get("PTA"), errors="coerce")

    for col in raw.columns:
        if m := _EAR_COL.match(col):
            ear = "B" if site in BINAURAL_ONLY_SITES else {"RH": "R", "LH": "L"}[m.group(1)]
            df[f"thr_{ear}_{m.group(2)}"] = pd.to_numeric(raw[col], errors="coerce")

    if site == "cogtail":
        parts = raw["subject_ID"].astype(str).str.extract(_COGTAIL_SUFFIX)
        df["group4"] = np.where(parts[0] == "TI", "TNHL", np.where(parts[0] == "NT", "HLOS", None))
        df["hearing_aid"] = (parts[1] == "HA").where(parts[1].notna())

    legacy = _legacy_id_map()
    df["legacy_id"] = [legacy.get((site, str(s))) for s in raw["subject_ID"]]
    return _finalise(df)


def husain_bids(site: str) -> pd.DataFrame:
    """UIUC / WHASC BIDS datasets (Husain lab): participants.tsv + per-subject sessions.tsv (ses-A)."""
    root = config.paths().sources[site]
    raw = pd.read_csv(root / "participants.tsv", sep="\t", na_values=["n/a"])
    df = pd.DataFrame(index=raw.index)
    df["original_id"] = raw["participant_id"].str.removeprefix("sub-")
    df["site"] = site
    df["participant_id"] = [ids.bids_id(site, s) for s in df["original_id"]]
    df["group4"] = raw["group"].str.strip().str.upper()
    df["group"] = df["group4"].map({"CTRL": "CO", "HLOS": "CO", "TINN": "TI", "TNHL": "TI"})
    df["sex"] = raw["sex"].map(_sex)
    df["age"] = pd.to_numeric(raw["age"], errors="coerce")
    tinnitus = df["group"] == "TI"
    # onset_duration is coded 0 for non-tinnitus groups; keep it only for patients.
    df["duration_months"] = pd.to_numeric(raw["onset_duration"], errors="coerce").where(tinnitus)
    for col in raw.columns:
        if m := re.match(r"^PTA_(\d+)$", col):
            df[f"thr_B_{m.group(1)}"] = pd.to_numeric(raw[col], errors="coerce")

    sess = []
    for orig in raw["participant_id"]:
        f = root / orig / f"{orig}_sessions.tsv"
        s = pd.read_csv(f, sep="\t", na_values=["n/a"]) if f.exists() else pd.DataFrame()
        row = s[s["session_id"] == "ses-A"].iloc[0] if len(s) else pd.Series(dtype=float)
        sess.append({k: pd.to_numeric(row.get(k), errors="coerce") for k in ("BDI", "BAI", "TFI", "TPFQ")})
    sess = pd.DataFrame(sess, index=raw.index)
    # TFI/TPFQ are recorded as 0 for non-tinnitus groups: not a real score.
    sess[["TFI", "TPFQ"]] = sess[["TFI", "TPFQ"]].where(tinnitus)
    df = pd.concat([df, sess], axis=1)
    return _finalise(df)


def antinomics_bids_ids() -> list[str]:
    """Subject labels present in the ANTINOMICS BIDS dataset (phenotypes come from the legacy sheet)."""
    root = config.paths().sources.get("antinomics_bids")
    if root is None or not Path(root).exists():
        return []
    return sorted(d.name.removeprefix("sub-") for d in Path(root).glob("sub-*") if d.is_dir())


def _legacy_id_map() -> dict[tuple[str, str], str]:
    f = config.paths().legacy_masters / "master.csv"
    if not f.exists():
        return {}
    m = pd.read_csv(f)
    return {(r.site, str(r.subject_ID)): r.tinception_id for r in m.itertuples()}


def audicog(site: str) -> pd.DataFrame:
    """AUDICOG (Hobeika, Paris): Data_audicog_behavior.xlsx, one row per participant.
    Per-ear audiogram 125 Hz–16 kHz; OAE/DPOAE columns are not imported (not used for structural MRI)."""
    raw = pd.read_excel(config.paths().sources["audicog_pheno"])
    df = pd.DataFrame(index=raw.index)
    df["original_id"] = raw["Num_participant_IRM"].astype(int).astype(str)
    df["site"] = site
    df["participant_id"] = [ids.bids_id(site, s) for s in df["original_id"]]
    df["group"] = raw["Group"].str.strip().map({"Control": "CO", "Tinnitus": "TI"})
    df["sex"] = raw["Sex"].str.strip().map({"Female": "F", "Male": "M"})
    df["age"] = pd.to_numeric(raw["Age"], errors="coerce")
    df["education_years"] = pd.to_numeric(raw["Number_years_education"], errors="coerce")
    df["THI"] = pd.to_numeric(raw["THI"], errors="coerce")
    df["duration_months"] = pd.to_numeric(raw["Length_tinnitus_month"], errors="coerce")
    df["HADS_A"] = pd.to_numeric(raw["HADS_Anxiety"], errors="coerce")
    df["HADS_D"] = pd.to_numeric(raw["HADS_Depression"], errors="coerce")
    df["HADS"] = df["HADS_A"] + df["HADS_D"]
    df["STAI"] = pd.to_numeric(raw["STAI"], errors="coerce")
    df["ISI"] = pd.to_numeric(raw["Insomnia_score_questionnaire_ISI"], errors="coerce")
    df["hyperacusis_khalfa"] = pd.to_numeric(raw["Khalfa"], errors="coerce")
    for col in raw.columns:
        if m := re.match(r"^PTA_(Left|Right)_(\d+)$", col):
            df[f"thr_{m.group(1)[0]}_{m.group(2)}"] = pd.to_numeric(raw[col], errors="coerce")
    return _finalise(df)


ADAPTERS = {"legacy_master": legacy_master, "husain_bids": husain_bids, "audicog": audicog}


def load_site(site: str) -> pd.DataFrame:
    source = config.sites()[site]["source"]
    return ADAPTERS[source](site)


def load_all(sites: list[str] | None = None) -> pd.DataFrame:
    sites = sites or list(config.sites())
    frames = [load_site(s) for s in sites]
    return pd.concat(frames, ignore_index=True, sort=False)
