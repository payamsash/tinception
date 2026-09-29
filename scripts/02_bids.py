"""Build the unified BIDS dataset (rawdata/) from all sources.

Modes (combine as needed):
  (default)  plan only: write rawdata/code/id_map.tsv and print what would be copied
  --copy     copy images + sidecars (never moves or deletes sources; skips existing files)
  --verify   md5-compare every copied image with its source

Run order: 02_bids (plan) -> 01_phenotypes (uses id_map for availability flags) -> 02_bids --copy --verify
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from tinception import bids, config


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def id_map(sources: list[bids.ImageSource]) -> tuple[pd.DataFrame, pd.DataFrame]:
    img = pd.DataFrame([
        {"participant_id": s.participant_id, "site": s.site, "original_id": s.original_id,
         "datatype": s.datatype, "suffix": s.suffix, "source": str(s.path),
         "legacy_drop": s.legacy_drop, "fs_dir": s.fs_dir}
        for s in sources
    ])
    per_sub = img.groupby("participant_id").agg(
        site=("site", "first"), original_id=("original_id", "first"),
        legacy_drop=("legacy_drop", "max"),
        fs_dir=("fs_dir", lambda s: next((v for v in s if v), None)),
        has_T1w=("suffix", lambda s: (s == "T1w").any()),
        has_T2w=("suffix", lambda s: s.str.endswith("T2w").any()),
        has_dwi=("datatype", lambda s: (s == "dwi").any()),
        has_rest=("suffix", lambda s: s.str.contains("task-rest").any()),
        t1_source=("source", lambda s: next((v for v in s if v.endswith("T1w.nii.gz") or "/anat/" not in v), None)),
    ).reset_index()
    return per_sub, img


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rawdata", type=Path, default=None)
    ap.add_argument("--datatypes", nargs="*", default=["anat"],
                    help="datatypes to copy (default: anat; the paper is structural-only)")
    ap.add_argument("--copy", action="store_true")
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args()

    raw = args.rawdata or config.paths().rawdata
    all_src = bids.all_sources()
    per_sub, img = id_map(all_src)  # availability flags reflect everything at the source
    sources = [s for s in all_src if s.datatype in set(args.datatypes)]
    (raw / "code").mkdir(parents=True, exist_ok=True)
    per_sub.to_csv(raw / "code" / "id_map.tsv", sep="\t", index=False, na_rep="n/a")
    img.to_csv(raw / "code" / "image_sources.tsv", sep="\t", index=False, na_rep="n/a")

    print(f"{len(per_sub)} subjects, {len(img)} images")
    print(per_sub.groupby("site").agg(n=("participant_id", "size"), T1w=("has_T1w", "sum"),
                                      fs=("fs_dir", lambda s: s.notna().sum()),
                                      dropped_v1=("legacy_drop", "sum"), T2w=("has_T2w", "sum"),
                                      dwi=("has_dwi", "sum"), rest=("has_rest", "sum")).to_string())
    need_fs = per_sub[per_sub.has_T1w & per_sub.fs_dir.isna()]
    need_fs[["participant_id", "site", "original_id"]].to_csv(raw / "code" / "needs_freesurfer.tsv", sep="\t", index=False)
    print(f"\nneeds FreeSurfer: {len(need_fs)} -> code/needs_freesurfer.tsv")
    # andes legacy surfaces used -T2pial; list for the T1-only sensitivity rerun
    per_sub[per_sub.site == "andes"][["participant_id", "site", "original_id"]].to_csv(
        raw / "code" / "andes_t1only.tsv", sep="\t", index=False)

    if args.copy:
        with open(raw / "dataset_description.json", "w") as f:
            json.dump(bids.dataset_description(), f, indent=2)
        for s in tqdm(sources, desc="copy"):
            bids.copy_source(s, raw)

    if args.verify:
        bad = [s for s in tqdm(sources, desc="verify")
               if not s.target(raw).exists() or md5(s.path) != md5(s.target(raw))]
        pd.DataFrame({"source": [str(s.path) for s in bad]}).to_csv(raw / "code" / "verify_failures.tsv", sep="\t", index=False)
        print(f"verify: {len(sources) - len(bad)}/{len(sources)} OK" + (f", {len(bad)} FAILED" if bad else ""))


if __name__ == "__main__":
    main()
