"""Apply an APPROVED cleanup manifest (from 00a_inventory) to the data root.

Only rows with approved == yes are touched. Dry-run by default; pass --execute to act.
Order: save MRIQC IQM jsons -> delete (children before parents) -> archive -> move.
RAW rows are never deleted here: raw originals are removed only after 02_bids --verify passes.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd

from tinception import config

ORDER = {"delete": 0, "archive": 1, "move": 2}


def _ignore_vanished(func, path, exc):
    # exFAT/macOS: AppleDouble '._*' companions disappear with their parent file
    if not isinstance(exc, FileNotFoundError):
        raise exc


def remove(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path, onexc=_ignore_vanished)
    else:
        path.unlink(missing_ok=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest", type=Path)
    ap.add_argument("--root", type=Path, default=None)
    ap.add_argument("--execute", action="store_true", help="actually modify the filesystem")
    ap.add_argument("--approve-all", action="store_true", help="treat every manifest row as approved")
    args = ap.parse_args()

    root = args.root or config.paths().root
    man = pd.read_csv(args.manifest, sep="\t", keep_default_na=False)
    if args.approve_all:
        man["approved"] = "yes"
    todo = man[(man.approved.str.lower() == "yes") & man.action.isin(ORDER)].copy()
    if (todo["class"] == "RAW").any() and (todo.action == "delete").any():
        raise SystemExit("refusing: RAW rows cannot be deleted by this script")
    todo["depth"] = todo.path.str.count("/")
    todo = todo.sort_values(by=["action", "depth"], key=lambda s: s.map(ORDER) if s.name == "action" else -s)

    mode = "EXECUTE" if args.execute else "DRY-RUN"
    qc = root / "qc_dir"
    if qc.exists() and ((todo.path == "qc_dir") & (todo.action == "delete")).any():
        dest = root / "derivatives" / "mriqc"
        jsons = [p for p in qc.rglob("*.json") if not p.name.startswith("._")]
        print(f"[{mode}] copy {len(jsons)} MRIQC jsons -> {dest}")
        if args.execute:
            for p in jsons:
                tgt = dest / p.relative_to(qc)
                tgt.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(p, tgt)

    for row in todo.itertuples():
        src = root / row.path
        if row.path.endswith("/*"):  # aggregated plain files of a directory
            files = [p for p in src.parent.iterdir() if p.is_file() and not p.name.startswith("._")]
            print(f"[{mode}] {row.action:7s} {len(files)} files in {src.parent.relative_to(root)}"
                  + (f" -> {row.dest}" if row.dest else ""))
            if args.execute:
                if row.action == "delete":
                    raise SystemExit("refusing to bulk-delete plain files; review manually")
                dst = root / row.dest
                dst.mkdir(parents=True, exist_ok=True)
                for p in files:
                    shutil.move(str(p), str(dst / p.name))
            continue
        if not src.exists():
            print(f"[{mode}] skip (missing) {row.path}")
            continue
        if row.action == "delete":
            print(f"[{mode}] delete  {row.path}  ({row.size_h})")
            if args.execute:
                remove(src)
        else:
            dst = root / row.dest
            print(f"[{mode}] {row.action:7s} {row.path} -> {row.dest}")
            if args.execute:
                if dst.exists():
                    raise SystemExit(f"destination exists: {dst}")
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(src), str(dst))


if __name__ == "__main__":
    main()
