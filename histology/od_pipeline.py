#!/usr/bin/env python3
"""
NPHS1 (nephrin) IHC — DAB optical-density pipeline
==================================================
Replicates the published NPHS1 image-analysis pipeline:
  RGB -> HED color deconvolution -> DAB channel -> optical density
  normalized 0-1 (1st-99th percentile) -> downsample x4 ->
  Gaussian smoothing (sigma=2) + Otsu threshold -> morphological
  filtering + border-object removal -> convex-hull ROI per glomerulus
  -> Mean_OD (mean DAB OD inside ROI).

Directory layout expected (one sub-folder per mouse, images inside):
    ROOT/
      AJcl1/  img1.tif img2.tif ...
      AJcl2/  ...
      BJclADR7/ ...

Outputs (written to --out, default ./od_out):
  od_per_glomerulus.csv   one row per detected glomerulus
  qc_overlays/<mouse>.png  detection overlay montage per mouse (visual QC)

Usage:
    python od_pipeline.py ROOT_DIR [--out OUT_DIR] [--folders AJcl1 AJcl2 ...]

Requirements: numpy scipy scikit-image tifffile pillow matplotlib
"""
import os
import sys
import glob
import csv
import argparse

import numpy as np
import tifffile
from PIL import Image
from skimage.color import rgb2hed
from skimage.filters import gaussian, threshold_otsu
from skimage.morphology import (remove_small_objects, binary_closing, disk,
                                convex_hull_image)
from skimage.measure import label, regionprops
from skimage.segmentation import clear_border
from skimage.transform import rescale
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Tunable parameters (match the published pipeline) ----------------------
DOWNSAMPLE = 4        # analysis at reduced resolution
GAUSS_SIGMA = 2       # smoothing before Otsu
MIN_AREA = 800        # min ROI area (px, in downsampled image)
MAX_ECC = 0.85        # eccentricity filter (drop elongated/tangential objects)
OD_PCTL = (1, 99)     # percentile range for OD normalization
IMG_EXT = (".tif", ".tiff", ".png", ".jpg", ".jpeg")


def load_rgb(path):
    try:
        a = tifffile.imread(path)
    except Exception:
        a = np.array(Image.open(path).convert("RGB"))
    if a.ndim == 2:
        a = np.stack([a] * 3, -1)
    if a.shape[-1] == 4:
        a = a[..., :3]
    return a


def dab_od(rgb):
    """DAB channel of HED deconvolution, normalized to 0-1 by percentile."""
    dab = rgb2hed(rgb)[:, :, 2]
    lo, hi = np.percentile(dab, OD_PCTL)
    return np.clip((dab - lo) / (hi - lo), 0, 1)


def analyze_image(path):
    """Return (downsampled_rgb, od_map, [glomerulus dicts], [bbox tuples])."""
    rgb = load_rgb(path)
    small = rescale(rgb, 1 / DOWNSAMPLE, channel_axis=2,
                    anti_aliasing=True, preserve_range=True).astype(np.uint8)
    od = dab_od(small)
    smooth = gaussian(od, sigma=GAUSS_SIGMA)
    mask = smooth > threshold_otsu(smooth)
    mask = binary_closing(mask, disk(3))
    mask = remove_small_objects(mask, min_size=MIN_AREA)
    mask = clear_border(mask)                    # drop border-touching glomeruli
    lbl = label(mask)

    results, boxes = [], []
    for rp in regionprops(lbl):
        if rp.area < MIN_AREA or rp.eccentricity > MAX_ECC:
            continue
        minr, minc, maxr, maxc = rp.bbox
        hull = convex_hull_image(lbl[minr:maxr, minc:maxc] == rp.label)
        mean_od = float(od[minr:maxr, minc:maxc][hull].mean())
        results.append({"area": int(rp.area), "mean_od": mean_od})
        boxes.append((minr, minc, maxr, maxc))
    return small, od, results, boxes


def save_overlay(mouse, montage, out_dir):
    """montage: list of (image_name, overlay_rgb, n_detected)."""
    n = len(montage)
    cols = min(4, n)
    rows = (n + cols - 1) // cols
    fig, ax = plt.subplots(rows, cols, figsize=(4 * cols, 3 * rows))
    ax = np.array(ax).reshape(-1)
    for k, (nm, ov, cnt) in enumerate(montage):
        ax[k].imshow(ov)
        ax[k].set_title(f"{nm} (n={cnt})", fontsize=8)
        ax[k].axis("off")
    for k in range(n, len(ax)):
        ax[k].axis("off")
    fig.suptitle(f"{mouse} — detected glomeruli", fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{mouse}.png"), dpi=90, bbox_inches="tight")
    plt.close()


def main():
    ap = argparse.ArgumentParser(description="NPHS1 DAB optical-density pipeline")
    ap.add_argument("root", help="Directory containing one sub-folder per mouse")
    ap.add_argument("--out", default="od_out", help="Output directory")
    ap.add_argument("--folders", nargs="*", default=None,
                    help="Specific mouse folders to process (default: all)")
    args = ap.parse_args()

    ov_dir = os.path.join(args.out, "qc_overlays")
    os.makedirs(ov_dir, exist_ok=True)

    folders = args.folders or sorted(
        d for d in os.listdir(args.root)
        if os.path.isdir(os.path.join(args.root, d)))

    rows = []
    for mouse in folders:
        d = os.path.join(args.root, mouse)
        imgs = sorted(f for f in glob.glob(os.path.join(d, "*"))
                      if f.lower().endswith(IMG_EXT))
        montage, gi = [], 0
        for img in imgs:
            small, od, res, boxes = analyze_image(img)
            ov = small.copy()
            for (minr, minc, maxr, maxc) in boxes:      # draw red bbox
                ov[minr:minr + 3, minc:maxc] = [255, 0, 0]
                ov[maxr - 3:maxr, minc:maxc] = [255, 0, 0]
                ov[minr:maxr, minc:minc + 3] = [255, 0, 0]
                ov[minr:maxr, maxc - 3:maxc] = [255, 0, 0]
            montage.append((os.path.basename(img), ov, len(res)))
            for r in res:
                gi += 1
                rows.append({"mouse": mouse, "image": os.path.basename(img),
                             "glom": gi, "area": r["area"],
                             "mean_od": round(r["mean_od"], 4)})
        if montage:
            save_overlay(mouse, montage, ov_dir)
        print(f"{mouse}: {gi} glomeruli", flush=True)

    csv_path = os.path.join(args.out, "od_per_glomerulus.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["mouse", "image", "glom", "area", "mean_od"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {len(rows)} glomeruli -> {csv_path}")
    print(f"QC overlays -> {ov_dir}")


if __name__ == "__main__":
    main()
