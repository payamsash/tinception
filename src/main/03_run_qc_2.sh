#!/bin/bash

export FREESURFER_HOME=/usr/local/freesurfer/8.0.0
export PATH=/usr/lib/mrtrix3/bin:$PATH
export PATH=/home/ubuntu/fsl/bin:$PATH
export PATH=/home/ubuntu/ants-2.6.2/ants-2.6.2/bin:$PATH

mriqc ./bids_data ./qc_dir participant --participant-label sub-01001
