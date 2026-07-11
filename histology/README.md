# NPHS1 IHC optical-density pipeline

`od_pipeline.py` — automated glomerular ROI detection and DAB mean optical
density quantification from NPHS1 (nephrin) immunohistochemistry images.

Pipeline (from the script's own docstring): RGB -> HED color deconvolution ->
DAB channel -> optical density normalized 0-1 (1st-99th percentile) ->
downsample x4 -> Gaussian smoothing (sigma=2) + Otsu threshold ->
morphological filtering + border-object removal -> convex-hull ROI per
glomerulus -> Mean_OD (mean DAB OD inside ROI).

## Relation to the Supplementary package

This matches, step for step, the automated glomerular ROI detection pipeline
described in the manuscript Methods (Otsu thresholding of a Gaussian-smoothed,
downsampled image; sigma=2; morphological filtering; removal of
border-touching objects; convex-hull ROI; DAB mean optical density) and was
previously flagged in `submission_supplementary/SUPPLEMENTARY_MANIFEST.txt`
as **Figure S2 (NPHS1 ROI overlay)** and **Figure S5 (WT1/NPHS1 two-way
ANOVA)** context: at that time, no image-analysis code for this pipeline was
present anywhere in the project, only the manuscript's textual description of
its logic.

**This resolves the "no code" half of that gap, not the "no data" half.**
Running this script still requires the original IHC image files (one
sub-folder per mouse, per `--help`), which are not present in this
repository or the wider project. Its output (`od_per_glomerulus.csv`,
per-mouse QC overlay PNGs) is what Figure S2 would be rendered from, and what
a WT1/NPHS1 two-way ANOVA (Figure S5) would need per-glomerulus/per-mouse
numeric input from — but neither the source images nor a prior run's output
CSV exist anywhere in this project as of this commit. If those images become
available, this script can be run directly:

```
python od_pipeline.py ROOT_DIR --out od_out
```
