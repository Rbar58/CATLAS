#!/usr/bin/env python3
"""
Cat head E-field model from a REAL whole-head MRI + the CATLAS atlas.

This is the pipeline that actually built the delivered mesh. It replaces the
older CATLAS-only synthesis (build_head_segmentation.py) by anchoring the model
on a subject's real MRI:

  Stage 1  register   CATLAS (CatT1avg + TPM brain) -> the real head MRI
                      (similarity, then B-spline deformable; metric focused on
                       the brain via a moving mask)
  Stage 2  segment    head/scalp from the real MRI; brain extent from the warped
                      atlas; CSF-vs-brain split from THIS subject's intensities;
                      skull as a geometric shell  [APPROXIMATE: T1 cannot image
                      bone -- use a CT for an accurate skull]
  Stage 3  mesh       SimNIBS create_mesh -> tetrahedral .msh

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
def stage2_segment(mri, brain_path, out):
    from sklearn.mixture import GaussianMixture
    img = nib.load(mri); vol = img.get_fdata().astype(np.float32); aff = img.affine
    vox = np.sqrt((aff[:3, :3] ** 2).sum(0))
    brain_warp = nib.load(brain_path).get_fdata() > 0.5

    sm = ndimage.gaussian_filter(vol, 0.8)
    fg = ndimage.binary_closing(sm > np.percentile(sm[sm > 0], 12), iterations=2)
    lab, _ = ndimage.label(fg)
    head = ndimage.binary_fill_holes(lab == (np.argmax(np.bincount(lab.flat)[1:]) + 1))

    brain = ndimage.binary_fill_holes(brain_warp & head)
    lab, _ = ndimage.label(brain)
    if lab.max() > 0:
        brain = lab == (np.argmax(np.bincount(lab.flat)[1:]) + 1)
    brain = ndimage.binary_fill_holes(ndimage.binary_closing(brain, iterations=2))

    vals = vol[brain].reshape(-1, 1)
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
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    brain = os.path.join(a.out, "brain_in_head.nii.gz")
    seg = os.path.join(a.out, "cat_tissues.nii.gz")
    if a.stage in ("all", "1"):
        print("Stage 1: register CATLAS -> head"); brain = stage1_register(a.mri, a.atlas_dir, a.out)
    if a.stage in ("all", "2"):
        print("Stage 2: segment");                 seg = stage2_segment(a.mri, brain, a.out)
    if a.stage in ("all", "3"):
        print("Stage 3: mesh (SimNIBS)");           stage3_mesh(seg, a.out)
    print("done.")
