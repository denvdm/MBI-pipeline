
"""
prep_fsaverage_to_pet.py
Register fsaverage volumetric atlas to MNI space, then resample to PET grid for R extraction.

Pipeline:
1) Read fsaverage T1/brain and aparc+aseg (label image).
2) Register fsaverage brain -> MNI T1 (rigid+affine, optional SyN).
3) Apply transforms to aparc+aseg with nearest-neighbour (stay in MNI space).
4) Resample transformed atlas to PET grid with nearest-neighbour.
5) (Optional) Resample GM probability map to PET grid with continuous interpolation.
6) (Optional) Export FreeSurfer LUT (label->name) to CSV.

Dependencies: ants (ANTsPy), nibabel, nilearn, (pandas if exporting LUT).

Example:
  python prep_fsaverage_to_pet.py \
    --pet FDG_template_3mm.nii.gz \
    --mni MNI152_T1_1mm.nii.gz \
    --fsavg-dir /path/to/freesurfer/subjects/fsaverage \
    --atlas-out aparc+aseg_inPETspace_3mm.nii.gz \
    --gm MNI152_GMprob_1mm.nii.gz \
    --gm-out GMprob_inPETspace_3mm.nii.gz \
    --lut /path/to/FreeSurferColorLUT.txt \
    --lut-out lut_aparc_aseg.csv \
    --with-syn  # optional; adds a deformable step for better alignment

Notes:
- This assumes fsaverage is already in an MNI-like space; we still run a concise registration for robustness.
- Use --with-syn if you want a deformable step; otherwise rigid+affine is usually enough for fsaverage->MNI.
"""

import argparse, os, sys, re, tempfile, shutil
import nibabel as nib
import numpy as np

import ants
from nilearn.image import resample_to_img

def parse_args():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--pet', required=True, help='PET template NIfTI (target grid), e.g., FDG_template_3mm.nii.gz')
    ap.add_argument('--mni', required=True, help='MNI T1 template NIfTI (e.g., MNI152_T1_1mm.nii.gz)')
    ap.add_argument('--fsavg-dir', required=True, help='FreeSurfer fsaverage SUBJECTS_DIR/fsaverage directory')
    ap.add_argument('--fsavg-brain', default=None, help='Override: fsaverage brain image (defaults to mri/brain.mgz or T1.mgz)')
    ap. add_argument('--fsavg-atlas', default=None, help='Override: fsaverage aparc+aseg label (defaults to mri/aparc+aseg.mgz)')
    ap.add_argument('--atlas-out', default='aparc+aseg_inPETspace_3mm.nii.gz', help='Output atlas resampled to PET grid')
    ap.add_argument('--gm', default=None, help='GM probability map (MNI space) to resample to PET grid (optional)')
    ap.add_argument('--gm-out', default='GMprob_inPETspace_3mm.nii.gz', help='Output GM prob resampled to PET grid')
    ap.add_argument('--lut', default=None, help='Path to FreeSurferColorLUT.txt (optional)')
    ap.add_argument('--lut-out', default='lut_aparc_aseg.csv', help='Where to write LUT CSV if --lut is given')
    ap.add_argument('--with-syn', action='store_true', help='Add deformable SyN step after affine')
    return ap.parse_args()

def ensure_exists(p):
    if not os.path.exists(p):
        raise FileNotFoundError(f"Missing file: {p}")
    return p

def pick_fsavg_paths(fsavg_dir, fsavg_brain=None, fsavg_atlas=None):
    if fsavg_brain is None:
        cand = [os.path.join(fsavg_dir, 'mri', 'brain.mgz'),
                os.path.join(fsavg_dir, 'mri', 'T1.mgz')]
        fsavg_brain = next((c for c in cand if os.path.exists(c)), None)
    if fsavg_atlas is None:
        fsavg_atlas = os.path.join(fsavg_dir, 'mri', 'aparc+aseg.mgz')
    if fsavg_brain is None or not os.path.exists(fsavg_brain):
        raise FileNotFoundError("Cannot find fsaverage brain image (mri/brain.mgz or mri/T1.mgz).")
    if not os.path.exists(fsavg_atlas):
        raise FileNotFoundError("Cannot find fsaverage aparc+aseg.mgz in mri/.")
    return fsavg_brain, fsavg_atlas

def load_as_ants(img_path):
    # ANTs can read NIfTI; for MGZ, convert via nibabel to a temp NIfTI
    if img_path.lower().endswith(('.nii', '.nii.gz')):
        return ants.image_read(img_path)
    else:
        nii_tmp = tempfile.NamedTemporaryFile(suffix='.nii.gz', delete=False)
        nii_tmp.close()
        img = nib.load(img_path)
        nib.save(nib.Nifti1Image(img.get_fdata().astype(np.float32), img.affine), nii_tmp.name)
        out = ants.image_read(nii_tmp.name)
        return out

def save_nifti_like(nparr, ref_img, out_path, dtype=np.float32):
    # Save numpy array with ref affine/header
    nii = nib.Nifti1Image(nparr.astype(dtype), ref_img.affine, ref_img.header)
    nib.save(nii, out_path)

def export_fs_lut(lut_path, out_csv):
    import pandas as pd
    rows = []
    with open(lut_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = re.split(r'\s+', s)
            if len(parts) >= 2 and parts[0].isdigit():
                idx = int(parts[0]); name = parts[1]
                rows.append({'label': idx, 'name': name})
    if not rows:
        raise RuntimeError("No label rows parsed from LUT; check the file format.")
    import pandas as pd
    df = pd.DataFrame(rows).drop_duplicates('label')
    df.to_csv(out_csv, index=False)
    print(f"[OK] Wrote LUT CSV: {os.path.abspath(out_csv)}  (n={len(df)})")

def main():
    args = parse_args()
    pet_img = nib.load(ensure_exists(args.pet))
    mni_img = nib.load(ensure_exists(args.mni))
    fs_brain_path, fs_atlas_path = pick_fsavg_paths(ensure_exists(args.fsavg_dir), args.fsavg_brain, args.fsavg_atlas)

    print(f"[INFO] fsaverage brain: {fs_brain_path}")
    print(f"[INFO] fsaverage atlas: {fs_atlas_path}")
    print(f"[INFO] MNI T1: {args.mni}")
    print(f"[INFO] PET template: {args.pet}")

    # --- Registration fsaverage brain -> MNI T1
    fs_brain_ants = load_as_ants(fs_brain_path)
    mni_ants = ants.image_read(args.mni)  # already NIfTI
    print("[INFO] Running rigid+affine registration (ANTs)...")
    reg_aff = ants.registration(fixed=mni_ants, moving=fs_brain_ants, type_of_transform='Affine')
    transforms = reg_aff['fwdtransforms']  # list of transforms to map fs->mni
    if args.with_syn:
        print("[INFO] Adding SyN deformable step...")
        reg_syn = ants.registration(fixed=mni_ants, moving=fs_brain_ants, type_of_transform='SyN', initial_transform=transforms)
        transforms = reg_syn['fwdtransforms']

    # --- Apply transforms to atlas (nearest-neighbour)
    print("[INFO] Applying transforms to aparc+aseg (NN)...")
    fs_atlas_ants = load_as_ants(fs_atlas_path)
    atlas_in_mni = ants.apply_transforms(fixed=mni_ants, moving=fs_atlas_ants,
                                         transformlist=transforms,
                                         interpolator='nearestNeighbor')
    atlas_in_mni_nii = nib.Nifti1Image(atlas_in_mni.numpy().astype(np.int16), mni_img.affine, mni_img.header)

    # --- Resample atlas to PET grid (NN)
    print("[INFO] Resampling transformed atlas -> PET grid (NN)...")
    atlas_resampled = resample_to_img(atlas_in_mni_nii, pet_img, interpolation='nearest')
    nib.save(atlas_resampled, args.atlas_out)
    print(f"[OK] Wrote atlas in PET grid: {os.path.abspath(args.atlas_out)}")

    # --- Optional: resample GM probability map to PET grid
    if args.gm:
        print("[INFO] Resampling GM probability map -> PET grid (continuous)...")
        gm_img = nib.load(ensure_exists(args.gm))
        gm_res = resample_to_img(gm_img, pet_img, interpolation='continuous')
        nib.save(gm_res, args.gm_out)
        print(f"[OK] Wrote GM prob in PET grid: {os.path.abspath(args.gm_out)}")
    else:
        print("[WARN] No GM map provided; you can run the R extractor with gm_prob=NULL.")

    # --- Optional LUT export
    if args.lut:
        try:
            export_fs_lut(args.lut, args.lut_out)
        except Exception as e:
            print(f"[WARN] LUT export failed: {e}")

    print("[DONE] Prep complete. Use outputs as inputs to the R extractor.")
    print(f"  atlas -> {args.atlas_out}")
    if args.gm:
        print(f"  gm    -> {args.gm_out}")
    if args.lut:
        print(f"  lut   -> {args.lut_out}")

if __name__ == '__main__':
    main()
