#!/bin/bash -l
# Install MATLAB Runtime R2019b (MCR v97) on the share; it is bind-mounted into the FreeSurfer
# image at $FREESURFER_HOME/MCRv97 so segmentHA_T1.sh / segmentThalamicNuclei.sh (the tools used
# for the 701 legacy subjects, outputs *-T1.v22 / v13) can run. Same URL as FS's fs_install_mcr.
#SBATCH --job-name=mcr
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/mcr_%j.out
set -euo pipefail
DEST=/shares/tinnitus.orl.med.uzh/containers/MCR_R2019b
TMP=/scratch/$USER/mcr_install; rm -rf "$TMP"; mkdir -p "$TMP" && cd "$TMP"
curl -sSL https://ssd.mathworks.com/supportfiles/downloads/R2019b/Release/9/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2019b_Update_9_glnxa64.zip -o installer.zip
unzip -q installer.zip
./install -mode silent -agreeToLicense yes -destinationFolder "$DEST"
ls "$DEST"/v97/bin/glnxa64/libmwlaunchermain.so && echo "MCR OK -> $DEST/v97"
rm -rf "$TMP"
