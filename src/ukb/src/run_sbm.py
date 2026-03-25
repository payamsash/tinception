import os
import shutil
import subprocess
import zipfile
from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm

from nipype import Node, Workflow
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.freesurfer import ReconAll, SurfaceTransform


def run_qcache_nipype(subject_ids, subjects_dir, work_dir, n_procs=5):
    infosource = Node(
        IdentityInterface(fields=["subject_id"]),
        name="infosource",
    )
    infosource.iterables = [("subject_id", subject_ids)]

    qcache = Node(ReconAll(), name="qcache")
    qcache.inputs.subjects_dir = str(subjects_dir)
    qcache.inputs.directive = "qcache"
    qcache.inputs.args = "-measure thickness -fwhm 10"

    wf = Workflow(name="fs_qcache", base_dir=str(work_dir))
    wf.connect(infosource, "subject_id", qcache, "subject_id")

    return wf.run(plugin="MultiProc", plugin_args={"n_procs": n_procs})


def run_surf_to_fsaverage5_nipype(subject_ids, subjects_dir, measure="thickness", fwhm=10):
    subjects_dir = Path(subjects_dir)

    for subject_id in subject_ids:
        for hemi in ["lh", "rh"]:
            in_file = subjects_dir / subject_id / "surf" / f"{hemi}.{measure}.fwhm{fwhm}.fsaverage.mgh"
            out_file = subjects_dir / subject_id / "surf" / f"{hemi}.{measure}.fwhm{fwhm}.fsaverage5.mgh"

            if not in_file.exists():
                print(f"missing input: {in_file}")
                continue

            sxfm = SurfaceTransform()
            sxfm.inputs.subjects_dir = str(subjects_dir)
            sxfm.inputs.hemi = hemi
            sxfm.inputs.source_subject = "fsaverage"
            sxfm.inputs.target_subject = "fsaverage5"
            sxfm.inputs.source_file = str(in_file)
            sxfm.inputs.out_file = str(out_file)

            print(f"running {subject_id} {hemi}")
            sxfm.run()


def ensure_dx_folder(remote_folder):
    subprocess.run(["dx", "mkdir", "-p", remote_folder], check=True)


def upload_hemi_outputs(subject_ids, subject_dirs, out_project_dir):
    for sid, subj_dir in zip(subject_ids, subject_dirs):
        remote_folder = f"{out_project_dir}/{sid}/surf"
        ensure_dx_folder(remote_folder)

        for hemi in ["lh", "rh"]:
            for mode in ["fsaverage", "fsaverage5"]:
                fname = subj_dir / "surf" / f"{hemi}.thickness.fwhm10.{mode}.mgh" 
                if fname.exists():
                    subprocess.run(
                        [
                            "dx",
                            "upload",
                            str(fname),
                            "--path",
                            f"{remote_folder}/",
                        ],
                        check=True,
                    )

def main(init_subject, end_subject):
    # ---- config ----
    license_fname = "license.txt"
    subprocess.run(["dx", "download", license_fname], check=True)
    os.environ["FS_LICENSE"] = "/opt/notebooks/license.txt"

    df_fname = "./results/ukb_vol.csv"
    fs_code = "20263"
    extension = "2_0.zip"
    batch_size = 8
    n_procs = 8
    out_project_dir = "/SBM"

    workdir = Path.cwd()
    subjects_dir = workdir / "subjects"
    subjects_dir.mkdir(exist_ok=True)

    nipype_work_dir = workdir / "nipype_work"
    nipype_work_dir.mkdir(exist_ok=True)

    os.environ["SUBJECTS_DIR"] = str(subjects_dir)

    # ---- input table ----
    subprocess.run(["dx", "download", df_fname, "-o", "ukb_vol.csv"], check=True)
    df = pd.read_csv("ukb_vol.csv")
    subjects = df["eid"].astype(str).tolist()[init_subject:end_subject]

    # ---- fsaverage ----
    for mode in ["fsaverage", "fsaverage5"]:
        fsaverage_fname = f"{mode}.zip"
        fsaverage_zip = subjects_dir / f"{mode}.zip"
        subprocess.run(["dx", "download", fsaverage_fname, "-o", str(fsaverage_zip)], check=True)
        with zipfile.ZipFile(fsaverage_zip, "r") as zf:
            zf.extractall(subjects_dir)
        fsaverage_zip.unlink(missing_ok=True)

    batch_subject_ids = []
    batch_subject_dirs = []
    batch_zip_files = []

    for subject in tqdm(subjects):
        folder_id = subject[:2]
        zip_name = f"{subject}_{fs_code}_{extension}"
        dx_path = f"/Bulk/Brain MRI/T1/{folder_id}/{zip_name}"

        local_zip = workdir / zip_name
        extract_dir = workdir / f"extract_{subject}"
        subject_id = str(subject).split("_")[0]
        final_subject_dir = subjects_dir / subject_id

        if local_zip.exists():
            local_zip.unlink()
        if extract_dir.exists():
            shutil.rmtree(extract_dir)
        if final_subject_dir.exists():
            shutil.rmtree(final_subject_dir)

        extract_dir.mkdir(exist_ok=True)

        try:
            subprocess.run(["dx", "download", dx_path, "-o", str(local_zip)], check=True)

            with zipfile.ZipFile(local_zip, "r") as zf:
                zf.extractall(extract_dir)

            extracted_fs_dir = extract_dir / "FreeSurfer"
            if not extracted_fs_dir.exists():
                print(f"missing extracted FreeSurfer dir for {subject}")
                local_zip.unlink(missing_ok=True)
                shutil.rmtree(extract_dir, ignore_errors=True)
                continue

            shutil.move(str(extracted_fs_dir), str(final_subject_dir))
            shutil.rmtree(extract_dir, ignore_errors=True)

            # remove the downloaded zip immediately after successful move
            local_zip.unlink(missing_ok=True)

            batch_subject_ids.append(subject_id)
            batch_subject_dirs.append(final_subject_dir)
            batch_zip_files.append(local_zip)

            if len(batch_subject_ids) == batch_size:
                print(f"running qcache with nipype on {batch_subject_ids}")

                try:
                    run_qcache_nipype(
                        subject_ids=batch_subject_ids,
                        subjects_dir=subjects_dir,
                        work_dir=nipype_work_dir,
                        n_procs=n_procs,
                    )
                except Exception as e:
                    print(f"qcache batch failed for {batch_subject_ids}: {e}")

                try:
                    run_surf_to_fsaverage5_nipype(
                        subject_ids=batch_subject_ids,
                        subjects_dir=subjects_dir
                    )
                except Exception as e:
                    print(f"surf_to_surf batch failed for {batch_subject_ids}: {e}")    

                upload_hemi_outputs(
                    subject_ids=batch_subject_ids,
                    subject_dirs=batch_subject_dirs,
                    out_project_dir=out_project_dir,
                )

                for subj_dir in batch_subject_dirs:
                    if subj_dir.exists():
                        shutil.rmtree(subj_dir, ignore_errors=True)

                for zf in batch_zip_files:
                    if zf.exists():
                        zf.unlink()

                batch_subject_ids = []
                batch_subject_dirs = []
                batch_zip_files = []

        except subprocess.CalledProcessError as e:
            print(f"error for subject {subject}: {e}")
            local_zip.unlink(missing_ok=True)
            shutil.rmtree(extract_dir, ignore_errors=True)
            shutil.rmtree(final_subject_dir, ignore_errors=True)

    if batch_subject_ids:
        print(f"running final qcache batch with nipype on {batch_subject_ids}")

        try:
            run_qcache_nipype(
                subject_ids=batch_subject_ids,
                subjects_dir=subjects_dir,
                work_dir=nipype_work_dir,
                n_procs=min(n_procs, len(batch_subject_ids)),
            )
        except Exception as e:
            print(f"final qcache batch failed for {batch_subject_ids}: {e}")

        try:
            run_surf_to_fsaverage5_nipype(
                subject_ids=batch_subject_ids,
                subjects_dir=subjects_dir
            )
        except Exception as e:
            print(f"surf_to_surf batch failed for {batch_subject_ids}: {e}") 

        upload_hemi_outputs(
            subject_ids=batch_subject_ids,
            subject_dirs=batch_subject_dirs,
            out_project_dir=out_project_dir,
        )

        for subj_dir in batch_subject_dirs:
            if subj_dir.exists():
                shutil.rmtree(subj_dir, ignore_errors=True)

        for zf in batch_zip_files:
            if zf.exists():
                zf.unlink()

    print("done")
    
if __name__ == "__main__":
    main(init_subject=23, end_subject=100)