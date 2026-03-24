import os
import shutil
import subprocess
import zipfile
from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm


def main():
    df = pd.read_csv("ukb_vol.csv")
    subjects = df["eid"].astype(str).head(20).tolist()

    fs_code = "20263"
    extension = "2_0.zip"
    batch_size = 10

    workdir = Path.cwd()
    subjects_dir = workdir / "subjects"
    subjects_dir.mkdir(exist_ok=True)

    out_project_dir = "/SBM"
    os.environ["SUBJECTS_DIR"] = str(subjects_dir)

    subprocess.run(["dx", "mkdir", "-p", out_project_dir], check=True)

    fsaverage_zip = Path("fsaverage.zip")
    with zipfile.ZipFile(fsaverage_zip, "r") as zf:
        zf.extractall(subjects_dir)

    batch_subject_ids = []
    batch_subject_dirs = []
    batch_zip_files = []

    for subject in tqdm(subjects):
        folder_id = subject[:2]
        zip_name = f"{subject}_{fs_code}_{extension}"
        dx_path = f"/Bulk/Brain MRI/T1/{folder_id}/{zip_name}"

        local_zip = workdir / zip_name
        extract_dir = workdir / f"extract_{subject}"
        subject_id = subject.split("_")[0]
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

            batch_subject_ids.append(subject_id)
            batch_subject_dirs.append(final_subject_dir)
            batch_zip_files.append(local_zip)

            if len(batch_subject_ids) == batch_size:
                print(f"running qcache in parallel on {batch_subject_ids}")

                cmd = (
                    f"parallel --jobs {batch_size} "
                    f"'recon-all -sd {subjects_dir} -s {{}} -qcache -measure thickness -fwhm 10' "
                    f"::: {' '.join(batch_subject_ids)}"
                )
                subprocess.run(cmd, shell=True, check=True)

                for sid, subj_dir in zip(batch_subject_ids, batch_subject_dirs):
                    lh = subj_dir / 'surf' / 'lh.thickness.fwhm10.fsaverage.mgh'
                    rh = subj_dir / 'surf' / 'rh.thickness.fwhm10.fsaverage.mgh'

                    for hemi in [lh, rh]:
                        if hemi.exists():
                            subprocess.run(
                                [
                                    "dx", "upload", str(hemi),
                                    "--path", f"{out_project_dir}/{sid}/surf/{hemi.name}"
                                ],
                                check=True
                            )
                        else:
                            print(f"missing output: {hemi}")

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
        print(f"running final parallel batch on {batch_subject_ids}")

        cmd = (
            f"parallel --jobs {len(batch_subject_ids)} "
            f"'recon-all -sd {subjects_dir} -s {{}} -qcache -measure thickness -fwhm 10' "
            f"::: {' '.join(batch_subject_ids)}"
        )
        subprocess.run(cmd, shell=True, check=True)

        for sid, subj_dir in zip(batch_subject_ids, batch_subject_dirs):
            lh = subj_dir / 'surf' / 'lh.thickness.fwhm10.fsaverage.mgh'
            rh = subj_dir / 'surf' / 'rh.thickness.fwhm10.fsaverage.mgh'

            for hemi in [lh, rh]:
                if hemi.exists():
                    subprocess.run(
                        [
                            "dx", "upload", str(hemi),
                            "--path", f"{out_project_dir}/{sid}/surf/{hemi.name}"
                        ],
                        check=True
                    )
                else:
                    print(f"missing output: {hemi}")

        for subj_dir in batch_subject_dirs:
            if subj_dir.exists():
                shutil.rmtree(subj_dir, ignore_errors=True)

        for zf in batch_zip_files:
            if zf.exists():
                zf.unlink()

    print("done")


if __name__ == "__main__":
    main()


dx run app-swiss-army-knife \
  -iin="run_qcache.py" \
  -iin="results/ukb_vol.csv" \
  -iin="fsaverage.zip" \
  -iimage_file="freesurfer_qcache.tar.gz" \
  -icmd='set -euo pipefail; python3 run_qcache.py' \
  --instance-type mem3_ssd1_v2_x16 \
  --watch -y


Loading the Docker image container-J6B38V8J8Xp9yqqp8XppK2z4:dxjupyterlab-image-processing_5.2.0.tar.gz