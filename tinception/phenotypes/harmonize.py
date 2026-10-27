"""Build the master participants table from all adapters and derive harmonised scores."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .. import config
from ..io import assert_unique_ids
from . import adapters, archive, audiometry

DERIVED = ["PTA4", "PTA_HF", "PTA_HF468", "PTA_EHF", "PTA4_better_ear", "PTA4_asymmetry", "HL", "pta_source"]


def within_zscore(values: pd.Series, groups: pd.Series, reference: pd.Series | None = None) -> pd.Series:
    """z-score ``values`` within each level of ``groups`` (e.g. per instrument).

    If ``reference`` (bool mask) is given, mean/SD come from the reference rows only.
    """
    ref = pd.Series(True, index=values.index) if reference is None else reference.astype(bool)
    out = pd.Series(np.nan, index=values.index, dtype=float)
    for g in groups.dropna().unique():
        m = groups == g
        mu, sd = values[m & ref].mean(), values[m & ref].std(ddof=1)
        if sd and np.isfinite(sd):
            out[m] = (values[m] - mu) / sd
    return out


DISTRESS_ORDER = ["THI", "TFI", "TQ_goebel"]  # preference order when several exist


def distress_score(df: pd.DataFrame) -> pd.DataFrame:
    """Instrument-independent tinnitus distress (patients only).

    THI, TFI and the Goebel-Hiller TQ are not interchangeable, so each is z-scored within instrument
    across all patients who completed it; the first available in DISTRESS_ORDER is used.
    """
    ti = df["group"] == "TI"
    instrument = pd.Series(None, index=df.index, dtype=object)
    raw = pd.Series(np.nan, index=df.index)
    for inst in DISTRESS_ORDER:
        if inst not in df:
            continue
        take = instrument.isna() & df[inst].notna() & ti
        instrument[take] = inst
        raw[take] = df.loc[take, inst]
    return pd.DataFrame({"distress_instrument": instrument,
                         "distress_raw": raw,
                         "distress_z": within_zscore(raw, instrument)})


# archive sources allowed to add subjects absent from the site master sheet, and which ones
ADD_SITES = {"antinomics"}
ADD_SITES_WHITELIST = {"dmxi"}  # antinomics subject scanned (BIDS) but missing from antinomics_master.xlsx
LAST_CONFLICTS: pd.DataFrame | None = None  # archive-vs-master disagreements of the last build

FILL_COLS = ["THI", "TFI", "BDI", "BAI", "HADS_A", "HADS_D", "HADS", "duration_months", "education_years",
             "age", "sex", "group"]


def merge_archive(df: pd.DataFrame, arc: pd.DataFrame | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Add archive clinical variables. Archive values only fill missing master values; where both exist
    and differ, the master value is kept and the pair is returned in the conflicts table."""
    arc = archive.load_archive() if arc is None else arc
    arc = arc.drop_duplicates(["site", "original_id"])
    # sites whose archive source is a complete participant export may ADD subjects missing from the
    # site master sheet (e.g. antinomics dmxi); only subjects with an image are relevant (see ADD_SITES)
    from .. import ids
    have = set(zip(df["site"], df["original_id"]))
    new = arc[arc["site"].isin(ADD_SITES) & arc["group"].notna()
              & np.array([(s, o) not in have for s, o in zip(arc["site"], arc["original_id"], strict=True)])]
    new = new[new["original_id"].isin(ADD_SITES_WHITELIST)]
    if len(new):
        add = pd.DataFrame({"site": new["site"].values, "original_id": new["original_id"].values})
        add["participant_id"] = [ids.bids_id(s, o) for s, o in zip(add.site, add.original_id)]
        df = pd.concat([df, add], ignore_index=True)
    m = df.merge(arc, on=["site", "original_id"], how="left", suffixes=("", "__arc"), validate="many_to_one")
    # rows whose audiogram is (partly) filled from the archive: archive-only columns keep their name in m
    from_arc = pd.Series(False, index=m.index)
    for c in (c for c in arc if c.startswith("thr_")):
        from_arc |= m[c].isna() & m[f"{c}__arc"].notna() if c in df else m[c].notna()
    conflicts = []
    for c in archive.ARCHIVE_COLS + sorted(c for c in arc if c.startswith("thr_")):
        a = f"{c}__arc"
        if a not in m:
            continue
        if c in FILL_COLS or c.startswith("thr_"):
            both = m[c].notna() & m[a].notna()
            diff = both & (pd.to_numeric(m[c], errors="coerce").sub(pd.to_numeric(m[a], errors="coerce")).abs() > 0.5
                           if m[c].dtype != object else m[c].astype(str) != m[a].astype(str))
            for r in m[diff].itertuples():
                conflicts.append({"participant_id": r.participant_id, "site": r.site, "variable": c,
                                  "master": getattr(r, c), "archive": getattr(r, a)})
            if not c.startswith("thr_"):
                m[f"{c}_source"] = np.where(m[c].notna(), "master", np.where(m[a].notna(), "archive", None))
            m[c] = m[c].fillna(m[a])
        else:
            m[c] = m[a] if c not in df else m[c].fillna(m[a])
        m = m.drop(columns=a)
    # AUDICOG's Khalfa score lives in its own column -> into the harmonised hyperacusis fields
    if "hyperacusis_khalfa" in m:
        k = m["hyperacusis_khalfa"].notna() & m["hyperacusis_raw"].isna()
        m.loc[k, "hyperacusis_raw"] = m.loc[k, "hyperacusis_khalfa"]
        m.loc[k, "hyperacusis_instrument"] = "KHALFA"
    m["hyperacusis_z"] = within_zscore(pd.to_numeric(m["hyperacusis_raw"], errors="coerce"),
                                       m["hyperacusis_instrument"])
    m["audiogram_source"] = np.where(from_arc, "archive", None)  # e.g. erlangen MetadatenfürZürich audiogram
    return m, pd.DataFrame(conflicts, columns=["participant_id", "site", "variable", "master", "archive"])


# legacy antinomics master columns that do not match the raw .mat source (0% agreement, 2026-10-15) and a
# mislabelled frequency; applied to an ear that is kept from the legacy master because it has no .mat file
ANTINOMICS_LEGACY_BAD_HZ = (250, 2000)
ANTINOMICS_LEGACY_RENAME_HZ = {12500: 12000}


def apply_audiogram_overrides(df: pd.DataFrame, aud: pd.DataFrame | None = None) -> pd.DataFrame:
    """Replace the antinomics audiogram with its raw-source thresholds (.mat files), per ear. Every antinomics
    ear without a .mat file keeps its legacy values minus ANTINOMICS_LEGACY_BAD_HZ (12500 relabelled 12000)."""
    df = df.copy()
    src = df["audiogram_source"] if "audiogram_source" in df else pd.Series(None, index=df.index, dtype=object)
    df["audiogram_source"] = src.where(src.notna(), np.where(df.filter(like="thr_").notna().any(axis=1),
                                                             "site_master", None))
    key = df["site"].eq("antinomics")
    if not key.any():
        return df
    if aud is None:
        try:
            aud = archive.antinomics_audiometry()
        except (KeyError, FileNotFoundError, OSError) as e:
            print(f"WARNING: antinomics .mat audiometry unavailable ({e}); legacy audiograms used without "
                  f"the bad {ANTINOMICS_LEGACY_BAD_HZ} Hz columns")
            aud = pd.DataFrame(columns=["original_id"])
    aud = aud.set_index("original_id")
    oid = df.loc[key, "original_id"].str.lower()
    from_mat = pd.Series(False, index=oid.index)
    from_legacy = pd.Series(False, index=oid.index)
    for ear in ("R", "L"):
        new_cols = [c for c in aud if c.startswith(f"thr_{ear}_")]
        has = (aud.reindex(oid.values)[new_cols].notna().any(axis=1).to_numpy() if new_cols
               else np.zeros(len(oid), bool))
        rows, legacy = oid.index[has], oid.index[~has]
        df.loc[rows, [c for c in df if c.startswith(f"thr_{ear}_")]] = np.nan
        for c in new_cols:
            if c not in df:
                df[c] = np.nan
            df.loc[rows, c] = aud.loc[oid[has].values, c].to_numpy()
        for hz in ANTINOMICS_LEGACY_BAD_HZ:
            if f"thr_{ear}_{hz}" in df:
                df.loc[legacy, f"thr_{ear}_{hz}"] = np.nan
        for a, b in ANTINOMICS_LEGACY_RENAME_HZ.items():
            ca, cb = f"thr_{ear}_{a}", f"thr_{ear}_{b}"
            if ca in df:
                if cb not in df:
                    df[cb] = np.nan
                df.loc[legacy, cb] = df.loc[legacy, cb].fillna(df.loc[legacy, ca])
                df.loc[legacy, ca] = np.nan
        from_mat.loc[rows] = True
        from_legacy.loc[legacy] |= df.loc[legacy, [c for c in df if c.startswith(f"thr_{ear}_")]].notna().any(axis=1)
    label = np.select([from_mat & from_legacy, from_mat, from_legacy],
                      ["antinomics_mat+legacy", "antinomics_mat", "site_master"], default=None)
    df.loc[oid.index, "audiogram_source"] = label
    return df


def apply_corrections(df: pd.DataFrame) -> pd.DataFrame:
    """Apply reviewed corrections from config/phenotype_corrections.yaml; log them per subject."""
    import yaml

    f = config.CONFIG_DIR / "phenotype_corrections.yaml"
    df = df.copy()
    df["corrections"] = ""
    if not f.exists():
        return df
    for c in yaml.safe_load(f.read_text()).get("corrections", []):
        hit = df["participant_id"] == c["participant_id"]
        if not hit.any():
            raise KeyError(f"correction for unknown subject {c['participant_id']}")
        old = df.loc[hit, c["variable"]].iloc[0]
        df.loc[hit, c["variable"]] = c["value"]
        df.loc[hit, "corrections"] += f"{c['variable']}: {old} -> {c['value']}; "
    return df


def affect_score(df: pd.DataFrame) -> pd.Series:
    """Mean of available within-instrument z-scores of depression/anxiety scales (all subjects)."""
    cols = [c for c in config.analysis()["phenotypes"]["affect_instruments"] if c in df and df[c].notna().any()]
    if not cols:
        return pd.Series(np.nan, index=df.index)
    zs = pd.concat([within_zscore(df[c], df["site"].where(df[c].notna())) for c in cols], axis=1)
    return zs.mean(axis=1)


IMAGE_COLS = ["legacy_drop", "fs_dir", "has_T1w", "has_T2w", "has_dwi", "has_rest"]


def add_image_availability(df: pd.DataFrame, id_map: pd.DataFrame) -> pd.DataFrame:
    """Outer-join image availability (from 02_bids code/id_map.tsv). Subjects with images but no
    phenotype row are kept (group unknown) so participants.tsv covers every BIDS subject."""
    img = id_map[["participant_id", "site", "original_id"] + IMAGE_COLS].rename(
        columns={"site": "_site", "original_id": "_oid"})
    out = df.merge(img, on="participant_id", how="outer", validate="one_to_one", indicator=True)
    out["site"] = out["site"].fillna(out.pop("_site"))
    out["original_id"] = out["original_id"].fillna(out.pop("_oid"))
    out["phenotype_status"] = out.pop("_merge").map(
        {"both": "ok", "left_only": "no_image", "right_only": "no_phenotype"}).astype(str)
    for c in ["has_T1w", "has_T2w", "has_dwi", "has_rest", "legacy_drop"]:
        out[c] = out[c].fillna(False).astype(bool)
    return out


def build_master(sites: list[str] | None = None, id_map: pd.DataFrame | None = None,
                 with_archive: bool = True) -> pd.DataFrame:
    df = adapters.load_all(sites)
    assert_unique_ids(df)
    global LAST_CONFLICTS
    LAST_CONFLICTS = None
    if with_archive:
        df, LAST_CONFLICTS = merge_archive(df)
    df = apply_audiogram_overrides(df)
    df = apply_corrections(df)
    if id_map is not None:
        df = add_image_availability(df, id_map)

    aud = audiometry.summarise(df, config.analysis()["audiometry"])
    # sites without audiograms (iccac, inditms, triple) only have a reported PTA
    aud["pta_source"] = np.where(aud["PTA4"].notna(), "computed_pta4",
                                 np.where(df["PTA_reported"].notna(), "reported", None))
    aud["PTA4"] = aud["PTA4"].fillna(df["PTA_reported"])
    aud["HL"] = (aud["PTA4"] > config.analysis()["audiometry"]["hl_threshold_db"]).astype("boolean").where(aud["PTA4"].notna())
    df = pd.concat([df, aud[DERIVED]], axis=1)

    # 4-group label from PTA4 at every site (user decision 2026-10-01); the site's own design labels
    # (uiuc/whasc/cogtail) are kept in group4_recruited. No PTA4 -> fall back to the recruited label.
    hl = df["HL"].astype("float")  # 1/0/NaN
    key = df["group"].str.cat(hl.map({0.0: "nh", 1.0: "hl"}), sep="_")
    derived4 = key.map({"CO_nh": "CTRL", "CO_hl": "HLOS", "TI_nh": "TINN", "TI_hl": "TNHL"})
    df["group4_recruited"] = df["group4"]
    df["group4_source"] = np.where(derived4.notna(), "pta4",
                                   np.where(df["group4"].notna(), "recruitment", None))
    df["group4"] = derived4.fillna(df["group4"])
    # sensitivity grouping: hearing loss defined by high-frequency PTA (4+8 kHz)
    hf_thr = config.analysis()["audiometry"]["hl_hf_threshold_db"]
    df["HL_HF"] = (df["PTA_HF"] > hf_thr).astype("boolean").where(df["PTA_HF"].notna())
    key_hf = df["group"].str.cat(df["HL_HF"].astype("float").map({0.0: "nh", 1.0: "hl"}), sep="_")
    df["group4_hf"] = key_hf.map({"CO_nh": "CTRL", "CO_hl": "HLOS", "TI_nh": "TINN", "TI_hl": "TNHL"})

    df = pd.concat([df, distress_score(df)], axis=1)
    df["affect_z"] = affect_score(df)
    df["tinnitus"] = (df["group"] == "TI").astype(int).where(df["group"].isin(["CO", "TI"]))
    return df


def missingness_report(df: pd.DataFrame, columns: list[str] | None = None) -> pd.DataFrame:
    columns = columns or ["age", "sex", "group", "PTA4", "PTA_HF", "PTA_EHF", "THI", "TFI",
                          "distress_z", "BDI", "BAI", "affect_z", "duration_months", "group4"]
    avail = df.groupby("site")[columns].apply(lambda g: g.notna().mean().round(2))
    avail.insert(0, "n", df.groupby("site").size())
    avail.insert(1, "n_CO", df[df.group == "CO"].groupby("site").size())
    avail.insert(2, "n_TI", df[df.group == "TI"].groupby("site").size())
    return avail.fillna({"n_CO": 0, "n_TI": 0})
