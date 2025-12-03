#!/usr/bin/env bash

# Preprocess PET templates: resampling, alignment, and masking into analysis space.
# Produces region-level or voxel-level PET maps for downstream weighting.

### python code to register FS atlases to PET map (incl. first register fsaverage to MNI), incl. resampling GM mask. 
WD=/cluster/projects/p33/ukbio_users/dennisva/mbi
FS=/cluster/software/EL9/amd/zen/easybuild/software/FreeSurfer/5.3.0-centos6_x86_64
petmap=source-castrillon2023_desc-cmrglc_space-MNI152_res-3mm_feature.nii.gz

module load Python/3.13.1-GCCcore-14.2.0
#pip install --user typing_extensions nibabel nilearn pandas antspyx

python3 ${WD}/scripts/functions/prep_fsaverage_to_pet.py \
  --pet ${WD}/neuromaps/${petmap} \
  --mni ${WD}/neuromaps/MNI152_T1_1mm_brain.nii.gz \
  --fsavg-dir ${FS}/subjects/fsaverage \
  --atlas-out ${WD}/neuromaps/aparc+aseg_inPETspace_3mm.nii.gz \
  --gm ${WD}/neuromaps/avg152T1_gray.nii.gz \
  --gm-out ${WD}/neuromaps/GMprob_inPETspace_3mm.nii.gz \
  --lut ${FS}/FreeSurferColorLUT.txt \
  --lut-out ${WD}/neuromaps/lut_aparc_aseg.csv
# (optional) add: --with-syn