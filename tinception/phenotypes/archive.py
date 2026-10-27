"""Clinical variables from the legacy archive (Tinception_repo/metadata), one extractor per site.

Each extractor returns one row per archive subject with ``site``, ``original_id`` (already mapped to
the master ID) and any of the standard columns in ``ARCHIVE_COLS``. Mapping rules were verified
subject-by-subject (see plan file, 2026-10-01). Values only FILL missing master values; disagreements
are reported by ``conflicts`` and never overwrite the master.

Instrument notes
- hyperacusis: GÜF (0–45; iccac, neuropren, tinspect), Khalfa HQ (14 items 0–42; tinmeg, audicog),
  HQ-20 (20 items 0–2, 0–40; london_1), HQ of unstated scale (london_2, new-castle).
  Harmonised as hyperacusis_z within instrument (harmonize.py).
- loudness: loudness_vas100 (0–100 scale) vs loudness_nrs10 (0–10 numeric rating).
- laterality: left | right | bilateral | head | other; laterality_asym = L>R / R>L for asymmetric bilateral.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from .. import config

ARCHIVE_COLS = [
    "THI", "TFI", "TQ_goebel", "BDI", "BAI", "MDI", "HADS_A", "HADS_D", "HADS", "PHQ15", "PSQ_index",
    "hyperacusis_raw", "hyperacusis_instrument", "IHS", "duration_months", "duration_note",
    "laterality", "laterality_asym", "laterality_raw", "loudness_vas100", "loudness_nrs10",
    "handedness", "education_years", "education_level", "age", "sex", "group",
]


def _d():
    return config.paths().phenotypes_raw


def _num(s) -> pd.Series:
    return pd.to_numeric(pd.Series(s), errors="coerce")


def _frame(site: str, ids, **cols) -> pd.DataFrame:
    df = pd.DataFrame({"original_id": list(ids), **{k: list(v) for k, v in cols.items()}})
    df.insert(0, "site", site)
    return df


def _months_from_text(x) -> float:
    """'24 Monaten' -> 24, '4 Jahre' / '3 years' / '>15 yrs' -> years*12, 'Ca 34' -> 34."""
    if pd.isna(x):
        return np.nan
    s = str(x).lower().replace(",", ".")
    m = re.search(r"(\d+(?:\.\d+)?)", s)
    if not m:
        return np.nan
    v = float(m.group(1))
    if re.search(r"jahr|year|yr", s):
        v *= 12
    return v


def _lat_text(x) -> tuple[str | float, str | float]:
    """Free-text laterality (English/German) -> (laterality, asymmetry)."""
    if pd.isna(x):
        return np.nan, np.nan
    s = str(x).strip().lower()
    if s in ("l>r", "r>l", "r~>l"):
        return "bilateral", "R>L" if s.startswith("r") else "L>R"
    if s in ("l", "left"):
        return "left", np.nan
    if s in ("r", "right"):
        return "right", np.nan
    if "kopf" in s or "head" in s or "central" in s or "centre" in s:
        return "head", np.nan
    if "beiden" in s or "both" in s or "equal" in s or "bilat" in s:
        asym = "R>L" if "recht" in s or "more right" in s else "L>R" if "link" in s or "more left" in s else np.nan
        return "bilateral", asym
    if "more left" in s or "mainly left" in s:
        return "bilateral", "L>R"
    if "more right" in s or "mainly right" in s:
        return "bilateral", "R>L"
    if "link" in s or "left" in s:
        return "left", np.nan
    if "recht" in s or "right" in s:
        return "right", np.nan
    return "other", np.nan


# ---------------------------------------------------------------------------------------- sites
def iccac() -> pd.DataFrame:
    import pyreadstat
    df, _ = pyreadstat.read_sav(str(_d() / "Unipark Combined Data T1.sav"))
    master = pd.read_excel(config.paths().legacy_masters / "iccac_master.xlsx")["subject_ID"].astype(str)
    lut = {m.removeprefix("ICCAC_").lower(): m for m in master}
    ids = [lut.get(str(p).strip().lower()) for p in df["PatID"]]

    def items(prefix, n, minus=1):
        cols = [f"{prefix}{i:02d}" for i in range(1, n + 1)]
        return df[cols].apply(pd.to_numeric, errors="coerce") - minus

    thi = df[[f"THI_{i:02d}" for i in range(1, 26)]].apply(pd.to_numeric, errors="coerce").replace({1: 4, 2: 2, 3: 0})
    tfi = items("TFI_", 25).mean(axis=1, skipna=True) * 10
    bdi_cols = [c for c in df if re.fullmatch(r"BDI_[A-RTU]", c)]
    bdi_s = (_num(df["BDI_S_Vorläufig"]) - 1).where(_num(df["BDI_S_Beurteilung"]) != 1, 0)
    bdi = (df[bdi_cols].apply(pd.to_numeric, errors="coerce") - 1).sum(axis=1, min_count=len(bdi_cols)) + bdi_s
    bai = items("BAI_", 21).sum(axis=1, min_count=21)
    guef_items = df[[f"GÜF_{i:02d}" for i in range(1, 16)]].apply(pd.to_numeric, errors="coerce") - 1
    # prorate up to 2 missing items (mean imputation), else missing
    guef = (guef_items.mean(axis=1) * 15).where(guef_items.notna().sum(axis=1) >= 13).round(1)
    lat_codes = {1: ("right", np.nan), 2: ("bilateral", "R>L"), 3: ("bilateral", np.nan), 4: ("left", np.nan),
                 5: ("bilateral", "L>R"), 6: ("head", np.nan), 7: ("other", np.nan)}
    lat = [lat_codes.get(int(v), (np.nan, np.nan)) if pd.notna(v) else (np.nan, np.nan) for v in _num(df["TONF_04"])]
    hand = df["DEMOG_03_Händigkeit"].astype(str).str.strip().str.upper().map({"R": "R", "L": "L", "B": "B"})
    return _frame("iccac", ids,
                  THI=thi.sum(axis=1, min_count=25), TFI=tfi.round(1), BDI=bdi, BAI=bai,
                  hyperacusis_raw=guef, hyperacusis_instrument=["GUEF"] * len(df),
                  duration_months=[_months_from_text(x) for x in df["TONF_01"]],
                  laterality=[a for a, _ in lat], laterality_asym=[b for _, b in lat],
                  laterality_raw=_num(df["TONF_04"]), loudness_vas100=_num(df["TONF_07"]),
                  handedness=hand).dropna(subset=["original_id"])


def neuropren() -> pd.DataFrame:
    D = _d()
    ti = pd.read_excel(D / "neuropren_Datensatz_Tinnitus_Komplett.xlsx")
    co = pd.read_excel(D / "Neuropren_Datensatz_Kontrollen_Komplett.xlsx")
    gk = pd.read_excel(D / "Neuropren_Fragebögen_Tinnitus_gekürzt.xlsx")
    master = pd.read_excel(config.paths().legacy_masters / "neuropren_master.xlsx")["subject_ID"].astype(str)
    lut = {}
    for m in master:
        if k := re.search(r"^(?:S\d+_)?([A-Za-z]{4})_MRT", m):
            lut[k.group(1).lower()] = m
    fix = {"xkzv": "xkvz"}

    def key(c):
        c = str(c).strip().lower()
        return lut.get(fix.get(c, c))

    loud = dict(zip(gk["code"].map(key), _num(gk["Tinnitus Lautstärke Skala 1-100"])))
    hand = {"rechts": "R", "links": "L", "beidseitig": "B"}
    rows = []
    for d, is_ti in ((ti, True), (co, False)):
        ids = d["code"].map(key)
        lat = [_lat_text(x) for x in (d["Tinnitus Lokalisation"] if is_ti else [np.nan] * len(d))]
        rows.append(_frame(
            "neuropren", ids,
            THI=_num(d["THI Score"]) if is_ti else [np.nan] * len(d),
            TFI=_num(d["TFI_Score"]) if is_ti else [np.nan] * len(d),
            hyperacusis_raw=_num(d["Geräuschüberempfindlichkeitsfragebogen Score"]),
            hyperacusis_instrument=["GUEF"] * len(d),
            PSQ_index=_num(d["PSQ Score zwischen 0 und 100"]),
            duration_months=_num(d["Tinnitusdauer (Monate)"]) if is_ti else [np.nan] * len(d),
            laterality=[a for a, _ in lat], laterality_asym=[b for _, b in lat],
            laterality_raw=d["Tinnitus Lokalisation"] if is_ti else [np.nan] * len(d),
            loudness_vas100=[loud.get(i) for i in ids] if is_ti else [np.nan] * len(d),
            loudness_nrs10=_num(d["tin_loudness_normal"]) if is_ti else [np.nan] * len(d),
            handedness=d["Händigkeit"].astype(str).str.strip().str.lower().map(hand),
            education_level=d["Bildungsstand"]))
    return pd.concat(rows, ignore_index=True).dropna(subset=["original_id"])


def tinspect() -> pd.DataFrame:
    d = pd.read_excel(_d() / "TINSPECT_masterfile_final.xlsx", sheet_name="all_final")
    master = pd.read_excel(config.paths().legacy_masters / "tinspect_master.xlsx")["subject_ID"].astype(str)
    lut = {m.split("_")[1]: m for m in master}
    ids = d["id"].astype(str).str.strip().map(lut)
    dur = pd.to_numeric(d["duration"].replace({"n.a.": np.nan}), errors="coerce")
    return _frame("tinspect", ids, BDI=_num(d["bdi"]), BAI=_num(d["bai"]),
                  hyperacusis_raw=_num(d["guef"]), hyperacusis_instrument=["GUEF"] * len(d),
                  PSQ_index=_num(d["psq_total"]), duration_months=dur,
                  # codes 1/2/3 have no codebook (1≈left, 2≈right, 3≈bilateral, one contradiction) -> raw only
                  laterality_raw=_num(d["tinnitus_side"])).dropna(subset=["original_id"])


def erlangen() -> pd.DataFrame:
    d = pd.read_excel(_d() / "MetadatenfürZürich.xlsx")
    thr = {}
    for c in d.columns:  # per-ear audiogram L125..L16k / R125..R16k; '.' = missing
        if m := re.fullmatch(r"([LR])(\d+)(k?)", str(c)):
            hz = int(m.group(2)) * (1000 if m.group(3) else 1)
            thr[f"thr_{m.group(1)}_{hz}"] = _num(d[c].replace(".", np.nan))
    return _frame("erlangen", d["ID"].astype(str).str.strip(),
                  TQ_goebel=_num(d["TF gesamt"]).replace(99, np.nan), **thr)


def tinmeg() -> pd.DataFrame:
    D = _d()
    q21 = pd.read_csv(D / "Tinmeg_Questionnaire Data_2021.csv", sep=None, engine="python", dtype={"StudieID": str})
    q23 = pd.read_csv(D / "Tinmeg_Questionnaire Data_2023.csv", sep=None, engine="python", dtype={"StudieID": str})
    idmap = pd.read_excel(D / "TinMeg3_Payam.xlsx", dtype=str)
    lut23 = dict(zip(idmap["quest_id"].str.strip(), idmap["MEG_ID"].str.strip()))
    q21["original_id"] = q21["StudieID"].str.strip()
    q23["original_id"] = q23["StudieID"].str.strip().map(lut23)
    rows = []
    for q in (q21, q23):
        q = q.replace(999, np.nan)

        def items(cols, q=q):
            return q[cols].apply(pd.to_numeric, errors="coerce")

        hq = (items([f"HQ_{i}" for i in range(1, 15)]) - 1).sum(axis=1, min_count=14)
        ha = items([c for c in q if re.fullmatch(r"HADS_\d+_A", c)]).sum(axis=1, min_count=7)
        hd = items([c for c in q if re.fullmatch(r"HADS_\d+_D", c)]).sum(axis=1, min_count=7)
        psq = items([f"PSQ_{i}" for i in range(1, 31)])
        for i in (1, 7, 10, 13, 17, 21, 25, 29):
            psq[f"PSQ_{i}"] = 5 - psq[f"PSQ_{i}"]
        psq_index = ((psq.sum(axis=1, min_count=30) - 30) / 90 * 100).round(1)
        thi_cols = [f"VAR05_{i}" for i in range(1, 26)]
        thi = items(thi_cols).sum(axis=1, min_count=25) if all(c in q for c in thi_cols) else np.nan
        rows.append(_frame("tinmeg", q["original_id"], THI=thi if np.ndim(thi) else [np.nan] * len(q),
                           hyperacusis_raw=hq, hyperacusis_instrument=["KHALFA"] * len(q),
                           HADS_A=ha, HADS_D=hd, HADS=ha + hd, PSQ_index=psq_index,
                           handedness=_num(q["TSCHQ_3"]).map({1: "R", 2: "L", 3: "B"})))
    return pd.concat(rows, ignore_index=True).dropna(subset=["original_id"]).drop_duplicates("original_id")


def cogtail() -> pd.DataFrame:
    k = pd.read_excel(_d() / "Masterfile2.0.xlsx", sheet_name="Kognition", header=2)
    master = pd.read_excel(config.paths().legacy_masters / "cogtail_master.xlsx")["subject_ID"].astype(str)
    lut = {m.split("_")[0].lower(): m for m in master}
    keys = k["VP-code"].astype(str).str.strip().str.split().str[0].str.lower()
    out = _frame("cogtail", keys.map(lut), education_years=_num(k["Ausbildung (Jahre)"]))
    return out.dropna(subset=["original_id"]).drop_duplicates("original_id")


def london_1() -> pd.DataFrame:
    d = pd.read_excel(_d() / "Questionnaire Summary london.xls", header=None)
    T = d.iloc[:, 2:].T
    T.columns = d[0].ffill().astype(str) + "|" + d[1].astype(str)
    ids = []
    for x in d.iloc[0, 2:]:
        s = str(x).strip()
        if s == "AGu":
            ids.append("AG_e1")
        elif m := re.match(r"^(\w+)\s*\(cont\.?\)", s):
            ids.append(f"{m.group(1)}_C_e1")
        else:
            ids.append(f"{s}_e1")

    def row(i):
        return d.iloc[i, 2:].tolist()

    is_ctrl = [i.endswith("_C_e1") for i in ids]
    hq = [np.nan if (c and v == 0) else v for c, v in zip(is_ctrl, _num(row(44)))]  # 'HQ data missing' -> 0
    lat = [_lat_text(x) for x in row(19)]
    dur = [np.nan if c else _months_from_text(x) for c, x in zip(is_ctrl, row(17))]  # controls: matched initials
    return _frame("london_1", ids, THI=_num(row(43)), HADS_A=_num(row(45)), HADS_D=_num(row(46)),
                  HADS=_num(row(45)) + _num(row(46)), PHQ15=_num(row(47)),
                  hyperacusis_raw=hq, hyperacusis_instrument=["HQ20"] * len(ids),
                  duration_months=dur, laterality=[a for a, _ in lat], laterality_asym=[b for _, b in lat],
                  laterality_raw=row(19), loudness_nrs10=_num(row(21)),
                  handedness=pd.Series(row(14)).astype(str).str.strip().str.upper().map({"R": "R", "L": "L"}))


def new_castle() -> pd.DataFrame:
    d = pd.read_excel(_d() / "Questionnaire Summary NCL.xls", header=2).iloc[:, :22]
    d = d[d["Number"].notna()]
    ids = d["Number"].astype("int64").astype(str)
    lat = [_lat_text(x) for x in d["T Lat"]]
    is_ti = _num(d["T"]) == 1
    return _frame("new-castle", ids, THI=_num(d["THI"]),
                  hyperacusis_raw=_num(d["HQ"]), hyperacusis_instrument=["HQ_NCL"] * len(d),
                  # unit of 'T Dur' not stated (2–29): assumed YEARS -> months; flagged in duration_note
                  duration_months=(_num(d["T Dur"]) * 12).where(is_ti),
                  duration_note=np.where(is_ti & _num(d["T Dur"]).notna(), "unit assumed years", ""),
                  laterality=[a for a, _ in lat], laterality_asym=[b for _, b in lat], laterality_raw=d["T Lat"],
                  loudness_nrs10=_num(d["T Ld Ov"]))


def london_2() -> pd.DataFrame:
    d = pd.read_excel(_d() / "demographics_19Ppts_EachGroup.xlsx", header=None)
    out = []
    for data_rows, id_row, is_ti in ((range(1, 17), 33, True), (range(37, 43), 59, False)):
        # row label sits in column 0 or column 1 depending on the row
        block = {str(d.iloc[r, 0] if pd.notna(d.iloc[r, 0]) else d.iloc[r, 1]).strip(): d.iloc[r, 2:21].tolist()
                 for r in data_rows}
        mri = d.iloc[id_row, 2:21].tolist()
        ids = [f"{str(m).strip()}_FIL" if pd.notna(m) else None for m in mri]

        def get(label, block=block, ids=ids):
            for k, v in block.items():
                if label.lower() in k.lower():
                    return v
            return [np.nan] * len(ids)

        lat = [_lat_text(x) for x in get("Side of tinnitus")] if is_ti else [(np.nan, np.nan)] * len(ids)
        out.append(_frame(
            "london_2", ids,
            TFI=_num(get("Tinnitus functional index")) if is_ti else [np.nan] * len(ids),
            hyperacusis_raw=_num(get("Hyperacusis Questionnaire")),
            hyperacusis_instrument=["HQ_L2"] * len(ids), IHS=_num(get("I H S")),
            duration_months=(_num(get("tinnitus duration (years)")) * 12) if is_ti else [np.nan] * len(ids),
            laterality=[a for a, _ in lat], laterality_asym=[b for _, b in lat],
            laterality_raw=get("Side of tinnitus") if is_ti else [np.nan] * len(ids),
            loudness_nrs10=_num(get("Average tinnitus loudness")) if is_ti else [np.nan] * len(ids),
            handedness=pd.Series(get("handedness")).astype(str).str.strip().str.upper().str[:1].map(
                {"R": "R", "L": "L"})))
    return pd.concat(out, ignore_index=True).dropna(subset=["original_id"])


def triple() -> pd.DataFrame:
    p = pd.read_excel(_d() / "mir_triple_patients_descr.xlsx")
    nr = next(c for c in p if c.strip() == "Nr.")
    ids = [f"Triple_{int(n):02d}_patients" for n in p[nr]]
    return _frame("triple", ids, MDI=_num(p["MDI"]))


def inditms() -> pd.DataFrame:
    d = pd.read_excel(_d() / "indiTMS_descriptives_payam.xlsx")
    fix = {"2fjeu": "2fjen", "8kmc7": "8mkc7"}
    ids = ["indiTMS_" + fix.get(str(x).strip(), str(x).strip()) for x in d["test_session_ID"]]

    def tfi(v):  # two TFI cells were auto-converted to dates by Excel (23.2 -> 2026-02-23)
        if isinstance(v, pd.Timestamp | np.datetime64) or hasattr(v, "month"):
            v = pd.Timestamp(v)
            return float(f"{v.day}.{v.month}")
        return pd.to_numeric(v, errors="coerce")

    return _frame("inditms", ids, MDI=_num(d["MDI"]), TFI=[tfi(v) for v in d["TFI"]])


def _first_number(x) -> float:
    """'60', '50-60' (-> 55), 'ca. 16-17' (-> 16.5), '14 Jahre' -> number; else NaN."""
    if pd.isna(x):
        return np.nan
    nums = [float(n) for n in re.findall(r"\d+(?:\.\d+)?", str(x).replace(",", "."))]
    if not nums:
        return np.nan
    m = re.search(r"(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)", str(x))
    return (float(m.group(1)) + float(m.group(2))) / 2 if m else nums[0]


def antinomics() -> pd.DataFrame:
    """ANTINOMICS / TIDE Zurich Unipark export (project 995463), coding from the scanned codebook
    (Questionnaire_TIDE_Zurich_Codebook.pdf; -77/-66/-99 = not shown/empty).
    HADS items are in Unipark order (std 1,8,2,9,...) with item-specific scoring; HQ is a 3-point
    version (Ja/Manchmal/Nein -> 2/1/0, total 0-28) and is kept as its own hyperacusis instrument."""
    f = config.paths().sources["antinomics_questionnaire"]
    d = pd.read_csv(f, sep=None, engine="python")
    d = d.replace([-77, -66, -99, "-77", "-66", "-99"], np.nan)
    d["oid"] = d["Studiennummer"].astype(str).str.strip().str.strip("[]").str.lower()
    t = pd.to_datetime(d["datetime"], format="%Y-%m-%d %H:%M:%S")  # Unipark ISO timestamps; fail if that changes
    d = (d.assign(_t=t).sort_values("_t", kind="stable")
         .drop_duplicates("oid", keep="last").drop(columns="_t"))  # latest submission

    def items(prefix, n):
        return d[[f"{prefix}_{i}" for i in range(1, n + 1)]].apply(pd.to_numeric, errors="coerce")

    thi = items("THI", 25).replace({1: 4, 2: 2, 3: 0}).sum(axis=1, min_count=25)
    tfi = ((items("TFI", 25) - 1).mean(axis=1) * 10).where(items("TFI", 25).notna().sum(axis=1) >= 19).round(1)
    hq = (3 - items("HQ", 14)).sum(axis=1, min_count=14)
    hads_raw = items("HADS", 14)
    rev = {1: True, 2: True, 3: False, 4: False, 5: True, 6: True, 7: False, 8: True, 9: True, 10: False,
           11: True, 12: True, 13: False, 14: False}  # True: score = 4 - code, False: code - 1
    hads = pd.DataFrame({i: (4 - hads_raw[f"HADS_{i}"]) if r else (hads_raw[f"HADS_{i}"] - 1)
                         for i, r in rev.items()})
    ha = hads[[1, 4, 5, 8, 9, 12, 13]].sum(axis=1, min_count=7)
    hd = hads[[2, 3, 6, 7, 10, 11, 14]].sum(axis=1, min_count=7)
    psq = items("PSQ", 30).replace(0, np.nan)
    for i in (1, 4, 6, 8, 12, 14, 16, 19):  # positively worded in this version
        psq[f"PSQ_{i}"] = 5 - psq[f"PSQ_{i}"]
    psq_index = ((psq.sum(axis=1, min_count=30) - 30) / 90 * 100).round(1)
    lat_codes = {1: ("right", np.nan), 2: ("bilateral", "R>L"), 3: ("bilateral", np.nan), 4: ("left", np.nan),
                 5: ("bilateral", "L>R"), 6: ("head", np.nan), 7: ("other", np.nan)}
    lat = [lat_codes.get(int(v), (np.nan, np.nan)) if pd.notna(v) else (np.nan, np.nan)
           for v in _num(d["TSCHQ_location"])]
    dur = d["TSCHQ_duration"].map(_months_from_text)
    esit = _num(d["ESIT_SQ_Monate"].map(_first_number)).fillna(_num(d["ESIT_SQ_Jahre"].map(_first_number)) * 12)
    dur = dur.fillna(esit)
    annett = d[[f"Annett_{i:02d}" for i in range(1, 13)]].apply(pd.to_numeric, errors="coerce")
    n_right, n_valid = (annett == 2).sum(axis=1), annett.notna().sum(axis=1)
    hand = np.where(n_valid == 0, None, np.where(n_right >= 9, "R", np.where(n_right <= 3, "L", "B")))
    age = _num(d["Alter"])
    edu = d["AusbildungJahre"].map(_first_number)
    tin = _num(d["Tinnitus_present"]).map({1: "TI", 4: "TI", 2: "CO", 5: "CO"})
    out = _frame("antinomics", d["oid"],
                 THI=thi, TFI=tfi, HADS_A=ha, HADS_D=hd, HADS=ha + hd, PSQ_index=psq_index,
                 hyperacusis_raw=hq, hyperacusis_instrument=["HQ14_3PT"] * len(d),
                 duration_months=dur.where(tin == "TI"),
                 laterality=[a for a, _ in lat], laterality_asym=[b for _, b in lat],
                 laterality_raw=_num(d["TSCHQ_location"]),
                 loudness_vas100=d["TSCHQ_loudness"].map(_first_number).where(lambda v: v <= 100),
                 handedness=hand, education_years=edu.where(edu.between(5, 30)),
                 age=age.where(age.between(15, 99)),
                 sex=_num(d["Geschlecht"]).map({1: "F", 2: "M"}))
    out["group"] = tin.values
    return out


def antinomics_audiometry(root=None) -> pd.DataFrame:
    """Per-ear pure-tone thresholds from the ANTINOMICS audiometry software output
    (<root>/<id>/<ID> <L|R> ANTINOMICS_absThreshold <date>.mat; betweenRuns.thresholds over
    betweenRuns.variableList1 in Hz). Repeated runs: the latest experiment.date per ear is used.
    These REPLACE the legacy master values (whose 250 and 2000 Hz columns do not match the source)."""
    import scipy.io as sio

    root = Path(root or config.paths().sources["antinomics_audiometry"])
    runs = {}  # (subject, ear) -> (date, thresholds, file name); a subject may have several folders
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        sub = d.name.split("_")[0].lower()
        for f in sorted(d.glob("*.mat")):
            if f.name.startswith("._"):
                continue
            m = sio.loadmat(str(f), squeeze_me=True, struct_as_record=False)
            ear = str(m["experiment"].ear)[0].upper()
            when = pd.to_datetime(str(m["experiment"].date), format="%d-%b-%Y %H_%M_%S", errors="coerce")
            when = pd.Timestamp.min if pd.isna(when) else when  # undated runs never beat dated ones
            b = m["betweenRuns"]
            thr = dict(zip(np.atleast_1d(b.variableList1).astype(int), np.atleast_1d(b.thresholds).astype(float),
                           strict=True))
            if not any(np.isfinite(v) for v in thr.values()):
                continue
            k = (sub, ear)
            if k not in runs or (when, f.name) > runs[k][::2]:
                runs[k] = (when, thr, f.name)
    rows = {}  # subjects without a usable .mat are absent -> they keep their master audiogram
    for (sub, ear), (_, thr, _) in runs.items():
        rows.setdefault(sub, {}).update({f"thr_{ear}_{hz}": v for hz, v in thr.items()})
    out = pd.DataFrame.from_dict(rows, orient="index")
    out.index.name = "original_id"
    out = out.reset_index()
    out.insert(0, "site", "antinomics")
    return out


EXTRACTORS = {"antinomics": antinomics, "iccac": iccac, "neuropren": neuropren, "tinspect": tinspect, "erlangen": erlangen,
              "tinmeg": tinmeg, "cogtail": cogtail, "london_1": london_1, "new-castle": new_castle,
              "london_2": london_2, "triple": triple, "inditms": inditms}


def load_archive() -> pd.DataFrame:
    frames = []
    for site, fn in EXTRACTORS.items():
        try:
            frames.append(fn())
        except Exception as e:  # noqa: BLE001 - one broken source must not block the rest; reported loudly
            warnings.warn(f"archive extractor {site} failed: {type(e).__name__}: {e}", stacklevel=2)
    df = pd.concat(frames, ignore_index=True)
    for c in ARCHIVE_COLS:
        if c not in df:
            df[c] = np.nan
    thr = sorted(c for c in df if c.startswith("thr_"))
    return df[["site", "original_id", *ARCHIVE_COLS, *thr]]
