# amira-sbfsem-cell-mask

Preprocessing module (**Module 1**) for organelle segmentation on **SBF-SEM**
images, integrated into **Amira-Avizo 3D 2025.1**.

This module automatically isolates the cell on each slice of an SBF-SEM volume,
crops the images around the cell and discards the resin background. It greatly
reduces the data size (around 90%) and increases the proportion of informative
pixels, a prerequisite for a usable segmentation training.

Developed during an M2 Bioinformatics internship at the Centre de Biologie
Structurale (CBS), BSME team, University of Montpellier.

---

## Why this module

Raw SBF-SEM volumes are large (> 20 GB) and made up of 60–95% empty background.
During segmentation training, the random sampling of patches then falls mostly
on this background, and the model learns almost nothing useful. This module
solves the problem upstream: it segments the cell (Otsu thresholding +
morphology), re-centers it on a common canvas and removes the useless border,
while preserving the biological structures.

## What the pipeline does

For each slice of the volume:

1. **Otsu thresholding** on a downsampled version to separate the cell from the
   background (threshold determined automatically, no manual tuning).
2. **Morphological closing** and hole filling to obtain a full, continuous cell
   mask.
3. **Largest connected component** kept (the cell of interest).
4. **Bounding box** computed per slice, then a **common canvas** (the largest
   box plus a margin) to preserve the alignment of the 3D stack.
5. **Re-centering** of each cell on this canvas.
6. **Optional downsampling** of the output.

## Outputs

- Three TIFF series: masked images, cropped unmasked images, binary mask
- A quality-control (QC) panel
- Five analysis graphs
- Three CSV metric files

---

## Installation in Amira

Copy the three files to the following locations in your Amira installation
(adapt the root path to your setup):

| File | Destination |
|---|---|
| `src/sbfsem/cell_mask_core.py` | `<Amira>/share/python_modules/sbfsem/cell_mask_core.py` |
| `amira/script_objects/cell_mask.pyscro` | `<Amira>/share/python_script_objects/cell_mask.pyscro` |
| `amira/resources/cell_mask.rc` | `<Amira>/share/resources/cell_mask.rc` |

Restart Amira: a **"SBF-SEM Cell Mask v2"** button appears in the left panel.
Clicking it loads the module and shows its properties.

> Note: the path to `cell_mask.pyscro` is hard-coded in `cell_mask.rc`. If your
> Amira installation differs, update this path.

## Command-line usage

The pipeline core also runs standalone, without Amira:

```bash
python src/sbfsem/cell_mask_core.py --input /path/data --output /path/out --output-scale 2
```

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | Input folder of TIFFs (e.g. denoised stack) |
| `--output` | required | Output folder (must exist) |
| `--n-samples` | `0` | `0` = all slices; `>0` = test on N evenly-spaced slices |
| `--downsample` | `4` | Downsampling factor for mask computation |
| `--closing-radius` | `5` | Radius for morphological closing |
| `--padding` | `100` | Margin (pixels) around the maximum bounding box |
| `--output-scale` | `1` | `1` = full resolution, `2` = half, `4` = quarter |

---

## Dependencies

- Python 3
- numpy
- tifffile
- scikit-image
- scipy
- matplotlib
- pandas

Install: `pip install -r requirements.txt`

---

## Indicative performance

For 125 slices at 13664 x 13184 pixels, expect about 45 to 60 minutes of
processing. In Amira, the UI is frozen during processing (this is normal);
progress is shown in the status field and the console.

---

## Repository structure

```
amira-sbfsem-cell-mask/
├── README.md
├── LICENSE
├── requirements.txt
├── .gitignore
├── src/
│   └── sbfsem/
│       └── cell_mask_core.py      # pipeline logic (CLI + importable)
└── amira/
    ├── script_objects/
    │   └── cell_mask.pyscro        # Amira UI module
    └── resources/
        └── cell_mask.rc            # Amira left-panel button
```

---

## Context

M2 Bioinformatics internship, University of Montpellier.
Centre de Biologie Structurale (CBS), BSME team.
Supervisor: Patrick Bron.
