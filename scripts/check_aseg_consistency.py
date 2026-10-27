"""QC: compare thalamus/hippocampus volumes in aseg.presurf.mgz vs final aseg.mgz for every subject.
They should agree closely; large differences mean the FreeSurfer outputs are internally inconsistent."""
import sys
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

fs = Path(sys.argv[1])
out = Path(sys.argv[2])
labels = {"L_thal": 10, "R_thal": 49, "L_hippo": 17, "R_hippo": 53}
rows = []
for d in sorted(fs.glob("sub-*")):
    if "_failed" in d.name:
        continue
    r = {"participant_id": d.name}
    for tag in ("aseg.presurf", "aseg"):
        f = d / "mri" / f"{tag}.mgz"
        if not f.exists():
            continue
        img = nib.load(str(f))
        a = np.asarray(img.dataobj)
        v = float(np.prod(img.header.get_zooms()))
        for k, lab in labels.items():
            r[f"{tag}_{k}"] = (a == lab).sum() * v
        r[f"{tag}_mtime"] = pd.Timestamp(f.stat().st_mtime, unit="s").date()
    rows.append(r)
df = pd.DataFrame(rows)
for k in labels:
    df[f"ratio_{k}"] = df[f"aseg_{k}"] / df[f"aseg.presurf_{k}"]
df.to_csv(out, sep="\t", index=False)
print(df.filter(like="ratio_").describe().round(3).to_string())
