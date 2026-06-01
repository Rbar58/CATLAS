# Provenance — cat head mesh (`cat_head.msh`)

Everything needed to reproduce the delivered mesh from scratch. The pipeline is
deterministic (fixed random seeds), so a rerun reproduces the same model.

## Inputs

| Input | Role | Shape | Voxel (mm) | md5 (first 12) |
|-------|------|-------|-----------|----------------|
| `mmc.nii` (your whole-head T1) | subject MRI (fixed/target) | 192×192×144 | 0.4167×0.4167×0.40 | `1b32381cc7ef` |
| `CatT1avg.nii` (CATLAS) | atlas T1, brain-only (moving) | 150×192×96 | 0.5×0.5×0.5 | `83be7b4b02c3` |
| `TPM.nii` (CATLAS) | atlas tissue priors; brain mask = (GM+WM+CSF)≥0.5 | 150×192×96×4 | 0.5×0.5×0.5 | `a2b061074dbc` |

`CatT1avg.nii` and `TPM.nii` are tracked in this repo. `mmc.nii` is your
unpublished data and is intentionally **not** committed; supply it via `--mri`.

## Software (exact versions used)

| Package | Version |
|---------|---------|
| Python | 3.11.15 |
| SimNIBS | 4.6.0 |
| SimpleITK | 2.5.5 |
| numpy | 2.4.6 |
| scipy | 1.17.1 |
| scikit-learn | 1.8.0 |
| nibabel | 5.4.2 |
| antspyx (Stage-1 registration) | 0.6.3 |

OS: headless Linux x86_64. SimNIBS meshing was installed in a non-standard
(pip/wheel) way — see `install_simnibs_meshing.sh`. On your own machine the
normal SimNIBS installer is fine; the version is what matters.

## Pipeline (script: `realmri_pipeline.py`)

Reproducibility: the default ANTs `SyNRA` stage is essentially deterministic
for a fixed input/version; the GaussianMixture in Stage 2 uses `random_state=0`.
(The SimpleITK fallback uses sampling seeds 42 and 7.)

**Stage 1 — register CATLAS → subject** (`stage1_register_ants`, default)
- **N4 bias-field correction** of the subject MRI (`subject_n4.nii.gz`), reused
  by Stage 2 for the tissue split.
- **Anatomical header fix** (`anatomical_reorient`): the subject header had its
  A-P and S-I axes effectively swapped vs. the atlas, rotating the warped brain
  ~90°. The affine is rewritten so axis0 = L-R (−X), axis1 = A-P (+Y anterior),
  axis2 = S-I (+Z superior); data is untouched (header-only), so the voxel grid
  still matches `mmc.nii` and the mask is saved back on the original affine.
  Verified against anatomy in 3 orthogonal views.
- Brain mask in atlas space = `(TPM GM+WM+CSF) ≥ 0.5`, used as the *moving mask*
  (`mask_all_stages=True`) so every stage is driven by brain voxels only.
- **ANTs `SyNCC`** (cross-correlation metric, `reg_iterations=(120,120,80,40)`):
  rigid → affine → SyN deformable, on `CatT1avg.nii` (moving) → reoriented
  subject (fixed). Warp the atlas brain mask with `genericLabel`; keep largest
  component; close + fill holes. CC pulls the boundary onto the true tissue; the
  earlier `SyNRA` left the brain slightly small and ventrally shifted.
- Output: `brain_in_head.nii.gz` (atlas brain warped to subject; **≈28.7 mL**,
  physiologically correct for a cat).
- *Fallback* `stage1_register` (SimpleITK Similarity3D multi-start {0.6,0.8,1.0}
  + B-spline 6×6×6) is retained but less accurate (≈33 mL, over-includes
  ventrally); select with `--reg simpleitk`.

**Stage 2 — segment** (`stage2_segment`)
- **Head/scalp**: Gaussian-smooth (σ=0.8); air level estimated from the 8 FOV
  corners (14³ voxels each), threshold = max(mean+6σ, cornermax, 20); binary
  close; keep largest component; fill holes (3D + slice-wise per axis).
- **Brain**: warped atlas brain ∩ head, filled, largest component, closed.
- **Brain interior**: 2-class GaussianMixture on the subject's **N4-corrected**
  intensities (`subject_n4.nii.gz`) inside the brain → low mean = CSF (3), high
  mean = brain (2). GM/WM are **not** split (this T1 has ~no GM/WM contrast).
- **CSF rim**: 1-voxel dilation of brain not already brain → CSF (3).
- **Skull (APPROXIMATE)**: dilate (brain+CSF) by r = round(1.5 mm / mean voxel)
  voxels, intersect head, minus (brain+CSF) → bone (4). *Geometric shell — T1
  cannot image bone; use a CT for a real skull.*
- **Scalp**: remaining head voxels → 5.
- Output: `cat_tissues.nii.gz` (labels 2=brain 3=CSF 4=bone 5=scalp).

**Stage 3 — mesh** (`stage3_mesh`)
- `simnibs.mesh_tools.meshing.create_mesh(labels, affine, optimize=False,
  num_threads=4)` → `write_msh` → `cat_head.msh`.

## Output (delivered)

- `cat_head.msh` — tetrahedral mesh, **52,609 nodes / 288,718 tets**,
  region tags 2/3/4/5 (volumes) and 1002/1003/1004/1005 (surfaces).
- `cat_tissues.nii.gz` — the label volume it was meshed from.
- QC: `seg_qc.png`, `seg_qc_contour.png`, `mesh_qc_3d.png`.

## Reproduce

```bash
bash install_simnibs_meshing.sh    # headless meshing env (one-time)
export LD_LIBRARY_PATH="$(python3 -c 'import sys,os;print(os.path.join(sys.prefix,"lib"))'):$LD_LIBRARY_PATH"
python realmri_pipeline.py --mri /path/to/mmc.nii --atlas-dir .. --out ./out
# -> out/cat_tissues.nii.gz  and  out/cat_head.msh
```

Code state: branch `claude/modest-thompson-hGJI5`, commit `2e55b3f` (plus the
commit adding this file).
