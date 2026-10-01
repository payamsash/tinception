"""Inventory and classification of the legacy derivatives folder (Step A of the v2 plan).

Rules map a path (relative to the data root) to a class and a proposed action. Nothing here
modifies the filesystem; execution is done by scripts/00b_reorganise.py from an approved manifest.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path

CLASSES = ("RAW", "KEEP-DERIV", "ARCHIVE", "REGENERABLE", "DUPLICATE", "JUNK", "REVIEW")


@dataclass(frozen=True)
class Rule:
    pattern: str          # regex on the POSIX relative path
    klass: str
    action: str           # keep | move | archive | delete
    dest: str = ""        # destination relative to root (for move/archive)
    note: str = ""


RULES: list[Rule] = [
    Rule(r"(^|/)\._[^/]*$|(^|/)\.DS_Store$", "JUNK", "delete", note="macOS metadata"),
    Rule(r"^vbm_per_site/[^/]+/struc$", "REGENERABLE", "delete", note="FSL-VBM intermediates"),
    Rule(r"^vbm_per_site/[^/]+_t2/\*$", "RAW", "keep", note="original T2 -> copied into rawdata/ by 02_bids"),
    Rule(r"^vbm_per_site/[^/]+/\*$", "RAW", "keep", note="original T1 -> copied into rawdata/ by 02_bids; delete only after md5 verification"),
    Rule(r"^vbm_per_site/[^/]+/[^/]+$", "REGENERABLE", "delete", note="non-image leftovers (logs, slicesdir)"),
    Rule(r"^subjects_fs_dir$", "KEEP-DERIV", "move", "derivatives/freesurfer-8.0.0", "FreeSurfer 8.0.0 (all 701 subjects)"),
    Rule(r"^ssa_data$", "KEEP-DERIV", "move", "derivatives/mind_legacy", "MIND matrices, Schaefer-1000"),
    Rule(r"^qc_dir$", "REGENERABLE", "delete", note="MRIQC html/svg; IQM jsons are archived first (see qc_dir_jsons)"),
    Rule(r"^QC$", "DUPLICATE", "delete", note="byte-identical copy of vbm_per_site T1s"),
    Rule(r"^VBM$", "DUPLICATE", "delete", note="T1 copies + FSL-VBM outputs"),
    Rule(r"^(ssa_results|vbm_norm|sbm_norm|SBM|MBM|plots)$", "REGENERABLE", "delete"),
    Rule(r"^subcortical_roi/norm_models[^/]*$", "REGENERABLE", "delete"),
    Rule(r"^GWAS/raw$", "REGENERABLE", "delete", note="downloadable summary statistics"),
    Rule(r"^(subcortical_roi|GWAS)/[^/]+$", "ARCHIVE", "archive", "archive/legacy_results_2026/{parent}/{name}",
         "small final legacy results"),
    Rule(r"^(vbm_results|biotypes_gmm|biotypes|PRS|phenotypes|VBM_design)$", "ARCHIVE", "archive",
         "archive/legacy_results_2026/{name}", "small final legacy results"),
    Rule(r"^(rawdata|derivatives|archive|results)$", "KEEP-DERIV", "keep", note="v2 layout"),
]

# directories whose children are classified individually; plain files inside them are
# aggregated into one pseudo-entry "<dir>/*"
DESCEND = re.compile(r"^(|vbm_per_site|vbm_per_site/[^/]+|subcortical_roi|GWAS)$")


def classify(rel: str) -> Rule:
    for r in RULES:
        if re.search(r.pattern, rel):
            parent, _, name = rel.rpartition("/")
            dest = r.dest.format(name=name, parent=parent).removesuffix("/*")
            return Rule(r.pattern, r.klass, r.action, dest, r.note)
    return Rule("", "REVIEW", "keep", note="no rule matched: review manually")


def du(path: Path) -> tuple[int, int, float]:
    """(n_files, apparent bytes, newest mtime) for a file, a directory tree or a "<dir>/*" entry."""
    if path.name == "*":
        stats = [p.lstat() for p in path.parent.iterdir()
                 if p.is_file() and not p.name.startswith("._")]
        return len(stats), sum(s.st_size for s in stats), max((s.st_mtime for s in stats), default=0.0)
    if path.is_file() or path.is_symlink():
        st = path.lstat()
        return 1, st.st_size, st.st_mtime
    n = size = 0
    newest = 0.0
    for dirpath, _, files in os.walk(path):
        for f in files:
            try:
                st = os.lstat(os.path.join(dirpath, f))
            except OSError:
                continue
            n += 1
            size += st.st_size
            newest = max(newest, st.st_mtime)
    return n, size, newest


def walk_entries(root: Path, rel: str = ""):
    """Yield relative paths to classify: children of DESCEND dirs, everything else as a unit."""
    base = root / rel if rel else root
    children = sorted(base.iterdir(), key=lambda p: p.name)
    aggregate = rel != "" and any(c.is_file() and not c.name.startswith("._") for c in children)
    if aggregate:
        yield f"{rel}/*"
    for child in children:
        crel = f"{rel}/{child.name}" if rel else child.name
        if aggregate and child.is_file() and not child.name.startswith("._"):
            continue
        rule = classify(crel)
        if child.is_dir() and DESCEND.match(crel) and rule.klass != "JUNK":
            yield from walk_entries(root, crel)
        else:
            yield crel
