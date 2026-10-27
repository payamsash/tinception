"""Audiometric summary measures computed identically for every site.

Threshold columns follow ``thr_<ear>_<freqHz>`` with ear in {R, L, B}. ``B`` is a binaural value:
either reported directly by the site (UIUC/WHASC, cogtail) or the mean of R and L.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

THR = re.compile(r"^thr_(R|L|B)_(\d+)$")


def threshold_columns(df: pd.DataFrame, ear: str) -> dict[int, str]:
    out = {}
    for c in df.columns:
        m = THR.match(c)
        if m and m.group(1) == ear:
            out[int(m.group(2))] = c
    return dict(sorted(out.items()))


def fill_binaural(df: pd.DataFrame) -> pd.DataFrame:
    """Add thr_B_<f> = mean(R, L) where both ears are measured and no binaural value exists."""
    df = df.copy()
    right, left = threshold_columns(df, "R"), threshold_columns(df, "L")
    for f in sorted(set(right) & set(left)):
        mean = df[[right[f], left[f]]].mean(axis=1, skipna=False)
        col = f"thr_B_{f}"
        df[col] = df[col].fillna(mean) if col in df else mean
    return df


def pta(df: pd.DataFrame, freqs: list[int], ear: str = "B", min_freqs: int | None = None) -> pd.Series:
    """Mean threshold over ``freqs``. Returns NaN unless at least ``min_freqs`` are measured
    (default: all of them). Frequencies the site never measured count as missing."""
    cols = threshold_columns(df, ear)
    min_freqs = len(freqs) if min_freqs is None else min_freqs
    block = pd.DataFrame({f: df[cols[f]] if f in cols else np.nan for f in freqs}, index=df.index)
    n = block.notna().sum(axis=1)
    return block.mean(axis=1).where(n >= min_freqs)


def summarise(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Return PTA4 / PTA_HF / PTA_EHF (binaural), better-ear PTA4, interaural asymmetry and HL flag."""
    df = fill_binaural(df)
    p4, phf, pehf = cfg["pta4_freqs"], cfg["pta_hf_freqs"], cfg["pta_ehf_freqs"]
    out = pd.DataFrame(index=df.index)
    out["PTA4"] = pta(df, p4, "B")
    out["PTA_HF"] = pta(df, phf, "B")
    out["PTA_HF468"] = pta(df, cfg["pta_hf468_freqs"], "B")
    out["PTA_EHF"] = pta(df, pehf, "B", min_freqs=2)
    r4, l4 = pta(df, p4, "R"), pta(df, p4, "L")
    out["PTA4_better_ear"] = pd.concat([r4, l4], axis=1).min(axis=1, skipna=False)
    out["PTA4_asymmetry"] = (r4 - l4).abs()
    out["HL"] = (out["PTA4"] > cfg["hl_threshold_db"]).astype("boolean").where(out["PTA4"].notna())
    return out
