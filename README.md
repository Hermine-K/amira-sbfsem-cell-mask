# CellFocus

**Cell-focused preprocessing for SBF-SEM organelle segmentation.**

CellFocus isolates the cell on each slice of an SBF-SEM volume, re-centers it on
a common canvas, removes the resin background and reduces the data size by around
90% while preserving the biological structures. It runs both as a command-line
tool and as an integrated module inside **Amira-Avizo 3D 2025.1**.

Developed during an M2 Bioinformatics internship at the Centre de Biologie
Structurale (CBS), BSME team, University of Montpellier.

Repository: `amira-sbfsem-cell-mask` · Package: `cellfocus`

---

## Why this tool

Raw SBF-SEM volumes are large (> 20 GB) and made up of 60 to 95% empty
background. During segmentation training, random patch sampling then falls mostly
on this background and the model learns almost nothing useful. CellFocus solves
the problem upstream: it segments the cell (Otsu thresholding plus morphology),
re-centers it on a common canvas and removes the useless border, which greatly
increases the proportion of informative pixels.

## What the pipeline does

For each slice of the volume:

1. Otsu thresholding on a downsampled version to separate the cell from the
   background (threshold determined automatically). Both contrast polarities are
   supported: bright cell on dark background (default) or dark cell on bright
   background (`--cell-dark`).
2. Morphological closing and hole filling to obtain a full, continuous mask.
3. Largest connected component kept (the cell of interest).
4. Per-slice bounding box, then a common canvas (largest box plus a margin) to
   preserve the alignment of the 3D stack.
5. Re-centering of each cell on this canvas.
6. Optional output downsampling.

## Outputs

- Three TIFF series: masked images, cropped unmasked images, binary mask
- A quality-control (QC) panel
- Five analysis graphs (colorblind-safe palette, readable axes at any slice count)
- Three CSV metric files and a text summary

---

## Installation

### As a Python package

```bash
pip install -e .
```

This exposes both the importable API and a `cellfocus` command-line entry point.

### Inside Amira

Copy the three files to your Amira installation (adapt the root path):

| File | Destination |
|---|---|
| `src/cellfocus/cell_mask_core.py` | `<Amira>/share/python_modules/cellfocus/cell_mask_core.py` |
| `amira/script_objects/cellfocus.pyscro` | `<Amira>/share/python_script_objects/cellfocus.pyscro` |
| `amira/resources/cellfocus.rc` | `<Amira>/share/resources/cellfocus.rc` |

Restart Amira: a **"CellFocus - SBF-SEM Cell Mask"** button appears in the left
panel. Clicking it loads the module and shows its properties.

> Note: the path to `cellfocus.pyscro` is hard-coded in `cellfocus.rc`. If your
> Amira installation differs, update this path.

---

## Command-line usage

```bash
cellfocus --input /path/data --output /path/out --output-scale 2
```

or, without installing:

```bash
python src/cellfocus/cell_mask_core.py --input /path/data --output /path/out --output-scale 2
```

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | Input folder of TIFFs (e.g. denoised stack) |
| `--output` | required | Output folder (created if needed) |
| `--n-samples` | `0` | `0` = all slices; `>0` = test on N evenly-spaced slices |
| `--downsample` | `4` | Downsampling factor for mask computation |
| `--closing-radius` | `5` | Radius for morphological closing |
| `--padding` | `100` | Margin (pixels) around the maximum bounding box |
| `--output-scale` | `1` | `1` = full resolution, `2` = half, `4` = quarter |
| `--cell-dark` | off | Set when the cell is darker than the background |

### Python API

```python
from cellfocus import run_pipeline

result = run_pipeline("data/denoised", "data/out", output_scale=2)
print(result["reduction_pct"])
```

---

## Dependencies

Python 3.9+, numpy, tifffile, scikit-image, scipy, matplotlib, pandas.

```bash
pip install -r requirements.txt
```

## Tests

```bash
pip install -e ".[test]"
pytest
```

---

## Indicative performance

For 125 slices at 13664 x 13184 pixels, expect about 45 to 60 minutes of
processing. In Amira, the UI is frozen during processing (this is normal);
progress is shown in the status field and the console.

## Repository structure

```
amira-sbfsem-cell-mask/
├── README.md
├── LICENSE
├── requirements.txt
├── pyproject.toml
├── .gitignore
├── src/
│   └── cellfocus/
│       ├── __init__.py
│       └── cell_mask_core.py       # pipeline logic (CLI + importable)
├── amira/
│   ├── script_objects/
│   │   └── cellfocus.pyscro         # Amira UI module
│   └── resources/
│       └── cellfocus.rc             # Amira left-panel button
└── tests/
    └── test_pipeline.py
```

## Roadmap

- Multipage TIFF stack input (in addition to folders of individual slices)
- Automatic contrast-polarity detection (currently a manual `--cell-dark` switch)
- Performance: reuse the mask computed during scanning instead of recomputing it
- Documentation site (ReadTheDocs) and submission to the Journal of Open Source
  Software (JOSS)

---

## Context

M2 Bioinformatics internship, University of Montpellier.
Centre de Biologie Structurale (CBS), BSME team.
Supervisor: Patrick Bron.

## License

MIT, see [LICENSE](LICENSE).
