"""Read-only inventory of the legacy data root -> cleanup_manifest.tsv for review.

Columns: path, class, action, dest, n_files, bytes, size_h, newest, note, approved.
The `approved` column is empty; set it to `yes` for rows you accept before running 00b_reorganise.
Also writes nifti_duplicates.tsv (T1 copies grouped by content fingerprint).
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
from pathlib import Path

import pandas as pd

from tinception import config
from tinception.inventory import classify, du, walk_entries


def human(n: float) -> str:
    for unit in ("B", "K", "M", "G", "T"):
        if n < 1024:
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}P"


def fingerprint(path: Path, full: bool) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        if full:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        else:  # size + first and last MiB: fast and sufficient to flag duplicates
            h.update(str(path.stat().st_size).encode())
            h.update(f.read(1 << 20))
            f.seek(max(0, path.stat().st_size - (1 << 20)))
            h.update(f.read())
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=None, help="data root (default: paths.root)")
    ap.add_argument("--out", type=Path, default=Path("cleanup"), help="output directory")
    ap.add_argument("--full-md5", action="store_true", help="full md5 instead of fast fingerprint")
    args = ap.parse_args()

    root = args.root or config.paths().root
    args.out.mkdir(parents=True, exist_ok=True)

    rows = []
    for rel in walk_entries(root):
        r = classify(rel)
        n, size, newest = du(root / rel)
        rows.append(dict(path=rel, **{"class": r.klass}, action=r.action, dest=r.dest, n_files=n,
                         bytes=size, size_h=human(size),
                         newest=dt.datetime.fromtimestamp(newest, tz=dt.UTC).date() if newest else "",
                         note=r.note, approved=""))
        print(f"{r.klass:12s} {r.action:8s} {human(size):>8s}  {rel}")
    man = pd.DataFrame(rows)
    man.to_csv(args.out / "cleanup_manifest.tsv", sep="\t", index=False)

    # duplicate T1s across the known copy locations
    nii = [p for d in ("vbm_per_site", "QC", "VBM") if (root / d).exists()
           for p in (root / d).rglob("*.nii.gz") if not p.name.startswith("._")]
    dup = pd.DataFrame({"path": [str(p.relative_to(root)) for p in nii],
                        "bytes": [p.stat().st_size for p in nii],
                        "fingerprint": [fingerprint(p, args.full_md5) for p in nii]})
    dup["n_copies"] = dup.groupby("fingerprint")["path"].transform("size")
    dup.sort_values(["fingerprint", "path"]).to_csv(args.out / "nifti_duplicates.tsv", sep="\t", index=False)

    summary = man.groupby(["class", "action"])["bytes"].sum().map(human)
    print("\n", summary.to_string())
    print(f"\nduplicate NIfTI groups: {(dup.n_copies > 1).sum()} files in groups with >1 copy")
    print(f"manifest -> {args.out / 'cleanup_manifest.tsv'}  (fill `approved`=yes to accept rows)")


if __name__ == "__main__":
    main()
