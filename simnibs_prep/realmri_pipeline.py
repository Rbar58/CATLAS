#!/usr/bin/env python3
"""
Cat head E-field model from a REAL whole-head MRI + the CATLAS atlas.

This is the pipeline that actually built the delivered mesh. It replaces the
older CATLAS-only synthesis (build_head_segmentation.py) by anchoring the model
on a subject's real MRI:

  Stage 1  register   CATLAS (CatT1avg + TPM brain) -> the real head MRI
                      preferred path (ANTs): N4 bias correction -> ANATOMICAL
                      HEADER FIX (see below) -> SyNCC deformable, metric focused
                      on the brain via a moving mask
  Stage 2  segment    head/scalp from the real MRI; brain extent from the warped
                      atlas; CSF-vs-brain split from THIS subject's N4 intensities;
                      skull as a geometric shell  [APPROXIMATE: T1 cannot image
                      bone -- use a CT for an accurate skull]
  Stage 3  mesh       SimNIBS create_mesh -> tetrahedral .msh

ORIENTATION FIX (important): the delivered subject scan's NIfTI header had its
anterior-posterior and superior-inferior axes effectively swapped relative to
the atlas, which rotated the warped brain ~90 deg. Inspecting the anatomy fixed
the true voxel-axis meaning to (axis0 = L-R, axis1 = A-P anterior-high,
axis2 = S-I superior-high) and the affine is rewritten to match BEFORE
registration -- with both images anatomically consistent, SyN has no gross
rotation to fight. This mapping was determined for THIS subject; a different
scan may need a different (or no) correction -- always QC the orientation.

Labels (SimNIBS standard): 2=brain 3=CSF 4=bone 5=scalp
(GM/WM are intentionally merged into one brain compartment: in-vivo T1 here has
 negligible GM/WM contrast, so a split would be fabricated. Add WM later from a
 better scan, or by warping the CATLAS WM prior.)

Usage:
    python realmri_pipeline.py --mri /path/to/whole_head.nii \
        --atlas-dir /path/to/CATLAS --out ./out [--stage all]

Requires SimNIBS meshing available for Stage 3 (see install_simnibs_sandbox.sh).
"""
import argparse, os, time
import numpy as np, nibabel as nib
from scipy import ndimage


# ----------------------------------------------------------------------------
def stage1_register(mri, atlas_dir, out):
    import SimpleITK as sitk
    fixed  = sitk.ReadImage(mri, sitk.sitkFloat32)
    moving = sitk.ReadImage(os.path.join(atlas_dir, "CatT1avg.nii"), sitk.sitkFloat32)
    tpm = nib.load(os.path.join(atlas_dir, "TPM.nii")).get_fdata()
    brain = ((tpm[..., 0] + tpm[..., 1] + tpm[..., 2]) >= 0.5).astype(np.uint8)
    mmask = sitk.GetImageFromArray(np.transpose(brain, (2, 1, 0))); mmask.CopyInformation(moving)
    fixed_n, moving_n = sitk.Normalize(fixed), sitk.Normalize(moving)

    # --- similarity (rigid + single scale) so the brain can't shear/balloon ---
    def sim_run(scale0):
        init = sitk.CenteredTransformInitializer(
            fixed_n, moving_n, sitk.Similarity3DTransform(),
            sitk.CenteredTransformInitializerFilter.GEOMETRY)
        s = sitk.Similarity3DTransform(init)
        p = list(s.GetParameters()); p[6] = scale0; s.SetParameters(p)
        R = sitk.ImageRegistrationMethod()
        R.SetMetricAsMattesMutualInformation(50)
        R.SetMetricSamplingStrategy(R.RANDOM); R.SetMetricSamplingPercentage(0.2, 42)
        R.SetMetricMovingMask(mmask); R.SetInterpolator(sitk.sitkLinear)
        R.SetOptimizerAsRegularStepGradientDescent(2.0, 1e-4, 500, relaxationFactor=0.7)
        R.SetOptimizerScalesFromPhysicalShift()
        R.SetShrinkFactorsPerLevel([4, 2, 1]); R.SetSmoothingSigmasPerLevel([2, 1, 0])
        R.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
        R.SetInitialTransform(s, inPlace=False)
        return R.Execute(fixed_n, moving_n), R.GetMetricValue()

    best = None
    for sc in (0.6, 0.8, 1.0):
        t, m = sim_run(sc)
        if best is None or m < best[1]:
            best = (t, m, sc)
    sim = best[0]; print(f"  similarity: scale={sim.GetParameters()[6]:.3f} metric={best[1]:.4f}")

    # --- B-spline deformable refine on top of the similarity ---
    bspline = sitk.BSplineTransformInitializer(fixed_n, [6, 6, 6])
    R = sitk.ImageRegistrationMethod()
    R.SetMetricAsMattesMutualInformation(48)
    R.SetMetricSamplingStrategy(R.RANDOM); R.SetMetricSamplingPercentage(0.25, 7)
    R.SetMetricMovingMask(mmask); R.SetInterpolator(sitk.sitkLinear)
    R.SetMovingInitialTransform(sim); R.SetInitialTransform(bspline, inPlace=True)
    R.SetOptimizerAsLBFGSB(1e-5, 60, maximumNumberOfCorrections=5)
    R.SetShrinkFactorsPerLevel([3, 1]); R.SetSmoothingSigmasPerLevel([1, 0])
    R.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    bsp = R.Execute(fixed_n, moving_n)
    full = sitk.CompositeTransform([sim, bsp])
    wb = sitk.Resample(mmask, fixed, full, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
    out_brain = os.path.join(out, "brain_in_head.nii.gz")
    sitk.WriteImage(wb, out_brain)
    n = int(sitk.GetArrayFromImage(wb).sum())
    sp = np.prod(fixed.GetSpacing())
    print(f"  brain mask: {n} vox  ~{n*sp/1000:.1f} mL  -> {out_brain}")
    return out_brain


# ----------------------------------------------------------------------------
def anatomical_reorient(mri, out):
    """Rewrite the subject affine so the voxel axes carry their true anatomical
    meaning: axis0 = L-R (left = -X), axis1 = A-P (anterior = +Y), axis2 = S-I
    (superior = +Z). Returns the path to the reoriented (data-identical, header-
    only) image. NOTE: the axis assignment below was established by inspecting
    THIS subject's anatomy -- verify/adjust for other scans."""
    img = nib.load(mri); vol = img.get_fdata().astype(np.float32)
    sx, sy, sz = [float(z) for z in img.header.get_zooms()[:3]]
    M = np.array([[-sx, 0, 0], [0, sy, 0], [0, 0, sz]])
    A = np.eye(4); A[:3, :3] = M; A[:3, 3] = -M @ (np.array(vol.shape) / 2.0)
    path = os.path.join(out, "subject_reorient.nii.gz")
    nib.save(nib.Nifti1Image(vol, A), path)
    print(f"  reoriented header axcodes: {nib.aff2axcodes(A)}  -> {path}")
    return path


def stage1_register_ants(mri, atlas_dir, out, n4=True, reorient=True):
    """Preferred Stage 1: ANTs SyN focused on the brain via the atlas brain mask.
    Pre-steps that proved necessary on the real scan:
      * N4 bias-field correction (also reused for the Stage-2 tissue split), and
      * anatomical header fix (see module docstring) to remove a ~90 deg rotation.
    Then SyNCC (cross-correlation metric) -- more accurate at pulling the brain
    boundary onto the true tissue than SyNRA, which left the brain slightly small
    and ventrally shifted. Needs `pip install antspyx`."""
    import ants
    src = mri
    if n4:
        print("  N4 bias-field correction ...")
        fixed_n4 = ants.n4_bias_field_correction(ants.image_read(mri))
        src = os.path.join(out, "subject_n4.nii.gz")
        nib.save(nib.Nifti1Image(fixed_n4.numpy(), nib.load(mri).affine), src)
    if reorient:
        src = anatomical_reorient(src, out)

    fixed  = ants.image_read(src)
    moving = ants.image_read(os.path.join(atlas_dir, "CatT1avg.nii"))
    tpm = nib.load(os.path.join(atlas_dir, "TPM.nii")).get_fdata()
    bm = ((tpm[..., 0] + tpm[..., 1] + tpm[..., 2]) >= 0.5).astype("float32")
    bm_path = os.path.join(out, "atlas_brain.nii.gz")
    nib.save(nib.Nifti1Image(bm, nib.load(os.path.join(atlas_dir, "CatT1avg.nii")).affine), bm_path)
    moving_mask = ants.image_read(bm_path)
    reg = ants.registration(fixed=fixed, moving=moving, type_of_transform="SyNCC",
                            moving_mask=moving_mask, mask_all_stages=True,
                            reg_iterations=(120, 120, 80, 40))
    warped = ants.apply_transforms(fixed=fixed, moving=moving_mask,
                                   transformlist=reg["fwdtransforms"],
                                   interpolator="genericLabel")
    arr = warped.numpy() > 0.5
    lab, _ = ndimage.label(arr)
    if lab.max() > 0:
        arr = lab == (np.argmax(np.bincount(lab.flat)[1:]) + 1)
    arr = ndimage.binary_fill_holes(ndimage.binary_closing(arr, iterations=1))
    # the reorientation is header-only, so the voxel grid still matches the
    # original MRI -- save the mask on the original affine for downstream stages
    out_brain = os.path.join(out, "brain_in_head.nii.gz")
    nib.save(nib.Nifti1Image(arr.astype("uint8"), nib.load(mri).affine), out_brain)
    sp = np.prod(nib.load(mri).header.get_zooms()[:3])
    print(f"  ANTs brain mask: {int(arr.sum())} vox  ~{int(arr.sum())*sp/1000:.1f} mL  -> {out_brain}")
    return out_brain


# ----------------------------------------------------------------------------
def stage2_segment(mri, brain_path, out):
    from sklearn.mixture import GaussianMixture
    img = nib.load(mri); vol = img.get_fdata().astype(np.float32); aff = img.affine
    vox = np.sqrt((aff[:3, :3] ** 2).sum(0))
    brain_warp = nib.load(brain_path).get_fdata() > 0.5

    sm = ndimage.gaussian_filter(vol, 0.8)
    # air estimated from the 8 FOV corners (head fills the frame, so the full
    # border isn't air); head = anything well above that, largest blob, filled
    k = 14
    corners = np.concatenate([sm[i:j, a:b, c:d].ravel()
              for i, j in [(0, k), (-k, None)] for a, b in [(0, k), (-k, None)]
              for c, d in [(0, k), (-k, None)]])
    thr = max(corners.mean() + 6*corners.std(), corners.max(), 20.0)
    fg = ndimage.binary_closing(sm > thr, iterations=2)
    lab, _ = ndimage.label(fg)
    head = ndimage.binary_fill_holes(lab == (np.argmax(np.bincount(lab.flat)[1:]) + 1))
    for ax in range(3):   # slice-wise hole fill catches FOV-spanning concavities
        head = np.swapaxes(ndimage.binary_fill_holes(np.swapaxes(head, 0, ax)), 0, ax)

    brain = ndimage.binary_fill_holes(brain_warp & head)
    lab, _ = ndimage.label(brain)
    if lab.max() > 0:
        brain = lab == (np.argmax(np.bincount(lab.flat)[1:]) + 1)
    brain = ndimage.binary_fill_holes(ndimage.binary_closing(brain, iterations=2))

    # CSF/brain split on the N4 bias-corrected intensities if Stage 1 produced
    # them (cleaner than the raw T1, which carries a shading gradient)
    n4_path = os.path.join(out, "subject_n4.nii.gz")
    gmm_vol = nib.load(n4_path).get_fdata().astype(np.float32) if os.path.exists(n4_path) else vol
    vals = gmm_vol[brain].reshape(-1, 1)
    gmm = GaussianMixture(2, n_init=3, random_state=0).fit(vals)
    order = np.argsort(gmm.means_.ravel())          # low->CSF, high->brain
    to_tissue = {order[0]: 3, order[1]: 2}
    seg = np.zeros(vol.shape, np.uint8)
    seg[brain] = np.vectorize(to_tissue.get)(gmm.predict(vals)).astype(np.uint8)

    rim = ndimage.binary_dilation(brain, iterations=1) & ~brain
    seg[rim & (seg == 0)] = 3                         # CSF rim
    brain_csf = seg > 0
    r = max(1, int(round(1.5 / float(np.mean(vox)))))
    bone = ndimage.binary_dilation(brain_csf, iterations=r) & head & ~brain_csf
    seg[bone & (seg == 0)] = 4                        # skull (APPROXIMATE)
    seg[head & (seg == 0)] = 5                        # scalp

    out_seg = os.path.join(out, "cat_tissues.nii.gz")
    nib.save(nib.Nifti1Image(seg, aff), out_seg)
    print("  tissue voxels:", {int(k): int((seg == k).sum()) for k in (2, 3, 4, 5)})
    print("  ->", out_seg)
    return out_seg


# ----------------------------------------------------------------------------
def stage3_mesh(seg_path, out):
    from simnibs.mesh_tools.meshing import create_mesh
    from simnibs.mesh_tools.mesh_io import write_msh
    img = nib.load(seg_path)
    labels = img.get_fdata().astype("uint16")
    t0 = time.time()
    msh = create_mesh(labels, img.affine, optimize=False, num_threads=4)
    out_msh = os.path.join(out, "cat_head.msh")
    write_msh(msh, out_msh)
    print(f"  meshed in {time.time()-t0:.0f}s  nodes={msh.nodes.nr} "
          f"tets={(msh.elm.elm_type==4).sum()}  -> {out_msh}")
    return out_msh


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mri", required=True, help="whole-head MRI (NIfTI)")
    ap.add_argument("--atlas-dir", default=".", help="dir with CatT1avg.nii + TPM.nii")
    ap.add_argument("--out", default="./out")
    ap.add_argument("--stage", default="all", choices=["all", "1", "2", "3"])
    ap.add_argument("--reg", default="ants", choices=["ants", "simpleitk"],
                    help="Stage-1 registration backend (ants is far more accurate)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    brain = os.path.join(a.out, "brain_in_head.nii.gz")
    seg = os.path.join(a.out, "cat_tissues.nii.gz")
    if a.stage in ("all", "1"):
        print(f"Stage 1: register CATLAS -> head ({a.reg})")
        brain = (stage1_register_ants if a.reg == "ants" else stage1_register)(a.mri, a.atlas_dir, a.out)
    if a.stage in ("all", "2"):
        print("Stage 2: segment");                 seg = stage2_segment(a.mri, brain, a.out)
    if a.stage in ("all", "3"):
        print("Stage 3: mesh (SimNIBS)");           stage3_mesh(seg, a.out)
    print("done.")
