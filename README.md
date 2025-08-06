# Tinception

### Normative Gradient Deviations Reveal Individualized and Shared Brain Signatures of Tinnitus

### Abstract

Tinnitus is a heterogeneous condition lacking consistent neural markers across individuals. Using high-resolution structural MRI, we examined both shared and individual-specific brain alterations in tinnitus. Voxel- and surface-based morphometry revealed reduced gray matter in the subcallosal area and increased volume in the bilateral putamen, alongside cortical thinning and surface area reductions in attention-related regions. Subcortical analyses showed increased volume in the accessory basal amygdala and thalamic nuclei, as well as the parabrachial complex in the brainstem—part of the attention-arousal network. Gradient-based mapping uncovered a shift in macroscale cortical organization, separating auditory/somatosensory and default mode networks, and revealed atypical deviation patterns in individuals with tinnitus. Embedding these deviations showed distinct but overlapping clusters for controls and tinnitus participants. Our findings highlight both convergent and heterogeneous neuroanatomical profiles in tinnitus and suggest the potential of normative deviation charts for individualized assessment.

## Project Structure
```plaintext
.
├── paper.       # Codes to create paper figures
├── roi          # Codes for region of interest analysis
├── sbm          # Codes for surface based morphometry
├── ssa.         # Codes for gradient analysis on MIND networks
├── tools        # gereral tools used in different parts
├── vbm          # Codes for voxel based morphometry
├── venv         # Virtual environment for the project

```

## Installation
- Ensure your Python 3.9+ installation path is defined in your system's PATH.

- Clone the repository:
    ```bash
    git clone git@github.com:payamsash/tinception.git
    ```
### with venv
- Set up a virtual environment (optional but recommended):
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows: .\env\Scripts\activate
    pip install -r requirements.txt
    ```
### with Conda
- Create and activate a Conda environment:
   ```bash
    conda env create -f environment.yml
    conda activate rspv
