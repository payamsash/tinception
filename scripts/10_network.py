"""Network context (analysis plan 6.9 B; Lariviere et al. 2020): hub and epicentre models of group effect maps and
of individual deviation maps, using ENIGMA Toolbox normative HCP connectivity (DK 68 + 14 subcortical).

Maps: H1 main-model coefficients (adjusted d) for tinnitus, PTA4 and PTA_HF beyond PTA4 (DK thickness + aseg; tinnitus
also DK area), and the H3 confound-regressed Model B bootstrap-ratio maps (patients: core-set mode 1, hearing-only
modes 1 and 2). Individual maps: CO-reference deviations (DK thickness) of patients and held-out controls.
Outputs in <results>/10_network/: hub.tsv, epicentres.tsv, subject_epicentres.tsv, summary.json.

Example (cluster, repo root):  sbatch slurm/10_network.sbatch
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as st
from statsmodels.stats.multitest import multipletests

from tinception import config
from tinception.stats import group as g
from tinception.stats import network as net

ID = "participant_id"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None, help="default: <results>/10_network")
    args = ap.parse_args()
    cfg, p = config.analysis(), config.paths()
    ncfg, seed = cfg["network"], cfg["seed"]
    out = args.out or p.results / "10_network"
    out.mkdir(parents=True, exist_ok=True)
    conn = net.connectomes()
    coords, hemi, names = g.parcel_centroids(p.freesurfer / "fsaverage")
    if sorted(names) != sorted(conn["fc_ctx"].index):
        raise SystemExit("DK parcels of fsaverage and ENIGMA connectomes differ")
    spins_for = net.spin_cache(coords, hemi, names, seed)

    q_thr = cfg["group"]["fdr_q"]
    maps = net.effect_maps(p.results, cfg["h3"], ncfg)

    # group maps: hubs and epicentres
    hub, epi = [], []
    for name, (ctx, sub) in maps.items():
        x = ctx.reindex(names).dropna()
        if len(x) < len(names):
            print(f"{name}: {len(names) - len(x)} cortical regions absent (failed calibration), spins on the rest")
        sp = spins_for(tuple(x.index), ncfg["spin_n"])
        hub += [{"map": name, "level": "cortex", "n_regions": len(x), **h} for h in net.hub_test(x, conn, sp)]
        if sub is not None and len(sub):
            if len(sub) < len(net.SUBCORTEX):
                print(f"{name}: subcortical hub test on {len(sub)}/{len(net.SUBCORTEX)} regions")
            hub += [{"map": name, "level": "subcortex", **h}
                    for h in net.subcortical_hub_test(sub, conn, ncfg["spin_n"], seed)]
        for kind in ("fc", "sc"):
            e = net.epicentres(x, net.seed_profiles(list(x.index), conn, kind), sp)
            epi.append(e.assign(map=name, connectivity=kind))
        print(f"{name}: done", flush=True)
    hub, epi = pd.DataFrame(hub), pd.concat(epi, ignore_index=True)
    hub.to_csv(out / "hub.tsv", sep="\t", index=False)
    epi.to_csv(out / "epicentres.tsv", sep="\t", index=False)

    # individual deviation maps: share of patients vs held-out controls with an epicentre at each seed
    master = pd.read_csv(args.master or p.root / "tinception_master.tsv", sep="\t", na_values="n/a",
                         low_memory=False)[[ID, "site", "group"]]
    subj = []
    for model in ncfg["deviation_models"]:
        d = p.derivatives / "normative" / "CO" / model
        ok = pd.read_csv(d / "feature_calibration.tsv", sep="\t", index_col=0).ok
        kn = [n for n in names if ok.get(f"thick__{n}", False)]
        sp = spins_for(tuple(kn), ncfg["subject_spin_n"])
        z = pd.read_csv(d / "z_all.tsv", sep="\t", na_values="n/a").merge(master, on=ID)
        z = z.dropna(subset=[f"thick__{n}" for n in kn])
        groups = {"patients": z[z.group.eq("TI")], "controls": z[z.reference.astype(bool)]}
        if min(len(v) for v in groups.values()) == 0:
            print(f"SKIP subject epicentres {model}: {({k: len(v) for k, v in groups.items()})}")
            continue
        profs = {kind: net.seed_profiles(kn, conn, kind) for kind in ("fc", "sc")}
        flags = {k: net.subject_epicentres(v[[f"thick__{n}" for n in kn]].to_numpy(float), profs, sp,
                                           ncfg["alpha_subject"]) for k, v in groups.items()}
        for kind, prof in profs.items():
            for s_i, s_name in enumerate(prof.index):
                fp, fc_ = flags["patients"][kind][:, s_i], flags["controls"][kind][:, s_i]
                a, b = np.abs(fp) > 0, np.abs(fc_) > 0
                tab = [[a.sum(), (~a).sum()], [b.sum(), (~b).sum()]]
                subj.append({"model": model, "connectivity": kind, "seed": s_name,
                             "prop_patients": a.mean(), "prop_controls": b.mean(),
                             "prop_patients_pos": (fp > 0).mean(), "prop_patients_neg": (fp < 0).mean(),
                             "p_fisher": st.fisher_exact(tab).pvalue,
                             "p_binom_patients": st.binomtest(int(a.sum()), len(a), ncfg["alpha_subject"],
                                                              alternative="greater").pvalue,
                             "n_patients": len(a), "n_controls": len(b)})
        print(f"subject epicentres {model}: done", flush=True)
    if not subj:
        raise SystemExit("no subject epicentres computed (every deviation model skipped); hub/epicentres written")
    subj = pd.DataFrame(subj)
    for col in ("p_fisher", "p_binom_patients"):  # BH over seeds within model x connectivity
        subj[col.replace("p_", "q_")] = subj.groupby(["model", "connectivity"])[col].transform(
            lambda x: multipletests(x, method="fdr_bh")[1])
    subj.to_csv(out / "subject_epicentres.tsv", sep="\t", index=False)

    summary = {"hub": hub.to_dict("records"),
               "epicentres_q05": epi[epi.q < q_thr][["map", "connectivity", "seed", "r", "p_spin", "q"]]
               .to_dict("records"),
               "subject_patients_vs_controls_q05": subj[subj.q_fisher < q_thr][
                   ["model", "connectivity", "seed", "prop_patients", "prop_controls", "q_fisher"]].to_dict("records"),
               "subject_patients_above_alpha_q05": int((subj.q_binom_patients < q_thr).sum()),
               "subject_mean_prop": subj.groupby(["model", "connectivity"])[["prop_patients", "prop_controls"]]
               .mean().reset_index().to_dict("records")}
    (out / "summary.json").write_text(json.dumps(g.json_clean(summary), indent=2, allow_nan=False, default=str))
    print(json.dumps(g.json_clean(summary), indent=2, default=str))


if __name__ == "__main__":
    main()
