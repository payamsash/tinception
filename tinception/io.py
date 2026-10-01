"""Strict, ID-keyed table operations. Never align tables by row position."""

from __future__ import annotations

import pandas as pd

ID = "participant_id"


def assert_unique_ids(df: pd.DataFrame, key: str = ID) -> None:
    dup = df[key][df[key].duplicated()].unique()
    if len(dup):
        raise ValueError(f"duplicate {key}: {list(dup)[:10]} (n={len(dup)})")


def merge_on_id(
    left: pd.DataFrame,
    right: pd.DataFrame,
    key: str = ID,
    how: str = "left",
    require: str = "none",
) -> pd.DataFrame:
    """Merge two tables on the subject key with one-to-one validation.

    require: 'none' | 'left' (every left ID must be in right) | 'both' (identical ID sets).
    """
    assert_unique_ids(left, key)
    assert_unique_ids(right, key)
    lids, rids = set(left[key]), set(right[key])
    if require in ("left", "both") and (missing := lids - rids):
        raise KeyError(f"{len(missing)} IDs missing on the right, e.g. {sorted(missing)[:5]}")
    if require == "both" and (extra := rids - lids):
        raise KeyError(f"{len(extra)} IDs only on the right, e.g. {sorted(extra)[:5]}")
    overlap = (set(left.columns) & set(right.columns)) - {key}
    if overlap:
        raise ValueError(f"overlapping columns would be suffixed: {sorted(overlap)}")
    return left.merge(right, on=key, how=how, validate="one_to_one")
