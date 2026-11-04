"""Frozen pyment embeddings (analysis plan 6.11 B). Runs INSIDE the pinned pyment-predict container (python 3.13,
TensorFlow, pyment 5.1.0; no tinception import): builds each published model exactly as pyment-predict does, cuts it
at the bottleneck (global pooling over the last 64-channel convolution, `<name>_top_pool`, i.e. the input of the
prediction heads) and writes one 64-d embedding per subject from the same `mri/crop.mgz` the brain-age predictions
used (created by pyment-predict; subjects without a crop are skipped and listed).

  python 19b_pyment_embed.py <fastsurfer root> <out dir>
  -> <out dir>/embeddings_{multi-2025,reg-2025}.tsv (participant_id, emb00..emb63), embeddings_skipped.txt
"""

from __future__ import annotations

import os
import sys

import nibabel as nib
import numpy as np
import pandas as pd

MODELS = {"multi-2025": "sfcn-multi", "reg-2025": "sfcn-reg"}


def main() -> None:
    root, out = sys.argv[1], sys.argv[2]
    import tensorflow as tf
    from pyment.models.sfcn import sfcn_factory

    folders = sorted(f for f in os.listdir(root) if f.startswith("sub-") and os.path.isdir(os.path.join(root, f)))
    have = [f for f in folders if os.path.isfile(os.path.join(root, f, "mri", "crop.mgz"))]
    skipped = sorted(set(folders) - set(have))
    with open(os.path.join(out, "embeddings_skipped.txt"), "w") as fh:
        fh.write("\n".join(skipped) + ("\n" if skipped else ""))
    enc = {}
    for weights, arch in MODELS.items():
        model = sfcn_factory(arch)(weights=weights)
        pool = [layer for layer in model.layers if layer.name.endswith("_top_pool")]
        if len(pool) != 1:
            raise SystemExit(f"{weights}: expected one *_top_pool layer, found {[p.name for p in pool]}")
        enc[weights] = tf.keras.Model(model.inputs, pool[0].output)
    rows = {w: [] for w in MODELS}
    for f in have:   # same input as pyment-predict: crop.mgz -> get_fdata() -> batch axis
        img = np.expand_dims(nib.load(os.path.join(root, f, "mri", "crop.mgz")).get_fdata(), 0)
        for w, e in enc.items():
            rows[w].append([f.removesuffix("_T1w"), *map(float, e.predict(img, verbose=0)[0])])
    for weights, r in rows.items():
        dim = len(r[0]) - 1 if r else 0
        df = pd.DataFrame(r, columns=["participant_id", *[f"emb{i:02d}" for i in range(dim)]])
        if df.participant_id.duplicated().any() or not np.isfinite(df.iloc[:, 1:].to_numpy()).all():
            raise SystemExit(f"{weights}: duplicated subjects or non-finite embeddings")
        tmp = os.path.join(out, f".embeddings_{weights}.tsv")
        df.to_csv(tmp, sep="\t", index=False)
        os.replace(tmp, os.path.join(out, f"embeddings_{weights}.tsv"))
        print(weights, df.shape, "skipped", len(skipped), flush=True)


if __name__ == "__main__":
    main()
