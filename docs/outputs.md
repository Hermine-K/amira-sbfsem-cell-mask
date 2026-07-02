# Outputs

In the output folder, CellFocus produces:

- `masked/`, `unmasked/`, `masks/`: three TIFF series (masked images, cropped
  unmasked images, binary mask).
- `qc/`: a quality-control panel (`qc_panel.png`) and five analysis graphs
  (colorblind-safe palette, readable axes at any slice count).
- `stats/`: three CSV metric files and a text `summary.txt`.

## Analysis graphs

1. Bounding box size per slice.
2. Size reduction per slice.
3. Cell centroid drift across slices (justifies the re-centering).
4. Intensity stability inside the cell (signal preservation).
5. Outlier detection (IQR) and potential canvas reduction.
