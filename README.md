# CellFocus

[![Documentation](https://img.shields.io/badge/docs-cellfocus.readthedocs.io-674EA7?style=for-the-badge&logo=readthedocs&logoColor=white)](https://cellfocus.readthedocs.io)

![Python](https://img.shields.io/badge/python-3.9+-3776AB?logo=python&logoColor=white)
![License](https://img.shields.io/badge/license-MIT-green)

**Isolate, crop and re-center the main object in large grayscale image stacks.**

CellFocus is a standalone Python tool that isolates a single main object (for
example a cell) from a large, mostly empty background across an image stack,
crops the images around it, re-centers them on a common canvas and reduces the
data size by up to around 90%. It runs three ways: from the command line, as a
Python library, or as an integrated module inside Amira-Avizo 3D.

It was developed for SBF-SEM organelle-segmentation preprocessing at the CBS, but
nothing in the core is specific to Amira or to electron microscopy.

Repository: `amira-sbfsem-cell-mask` · Package: `cellfocus`

---

## When it applies

CellFocus fits any grayscale image or stack where a single dominant object can be
separated from the background by intensity. It assumes one main object (it keeps
the largest connected region), enough contrast for a global Otsu threshold, and
supports either polarity: bright object on dark background (default) or dark
object on bright background (`--cell-dark`). The input is grayscale, not binary;
CellFocus does the thresholding itself. It is an intensity-plus-morphology
heuristic, not semantic segmentation, so it does not fit low-contrast or textured
backgrounds, or cases where several objects must be kept.

## Why this tool

Raw SBF-SEM volumes are large (> 20 GB) and made up of 60 to 95% empty
background. During segmentation training, random patch sampling then falls mostly
on this background and the model learns almost nothing useful. CellFocus solves
the problem upstream: it segments the object, re-centers it on a common canvas
and removes the useless border, which greatly increases the proportion of
informative pixels.

## What the pipeline does

For each slice of the stack:

1. Otsu thresholding on a downsampled version to separate the object from the
   background (both contrast polarities supported).
2. Morphological closing and hole filling to obtain a full, continuous mask.
3. Largest connected component kept.
4. Per-slice bounding box, then a common canvas (largest box plus a margin) to
   preserve the alignment of the 3D stack.
5. Re-centering of each object on this canvas.
6. Optional output downsampling.

## Outputs

- Three TIFF series: masked images, cropped unmasked images, binary mask
- A quality-control (QC) panel
- Five analysis graphs (colorblind-safe palette, readable axes at any slice count)
- Three CSV metric files and a text summary

---

## Documentation

Full documentation (installation, usage, a worked example on real data, Amira
integration): built with MkDocs in the `docs/` folder. A hosted version on Read
the Docs will be linked here once published.

Preview locally:

```bash
pip install mkdocs-material
mkdocs serve
```

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
panel.

> Note: the path to `cellfocus.pyscro` is hard-coded in `cellfocus.rc`. Update it
> if your Amira installation differs.

## Command-line usage

```bash
cellfocus --input /path/to/tiffs --output /path/to/out --n-samples 5
```

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | Input folder of TIFFs (e.g. denoised stack) |
| `--output` | required | Output folder (created if needed) |
| `--n-samples` | `0` | `0` = all slices; `>0` = test on N evenly-spaced slices |
| `--downsample` | `4` | Downsampling factor for mask computation |
| `--closing-radius` | `5` | Radius for morphological closing |
| `--padding` | `100` | Margin (pixels) around the maximum bounding box |
| `--output-scale` | `1` | `1` = full resolution, `2` = half, `4` = quarter |
| `--cell-dark` | off | Set when the object is darker than the background |

### Python API

```python
from cellfocus import run_pipeline

result = run_pipeline("data/denoised", "data/out", output_scale=2)
print(result["reduction_pct"])
```

## Tests

```bash
pip install -e ".[test]"
pytest
```

---

## Repository structure

```
amira-sbfsem-cell-mask/
├── README.md
├── LICENSE
├── requirements.txt
├── requirements-docs.txt
├── pyproject.toml
├── mkdocs.yml
├── .readthedocs.yaml
├── .gitignore
├── src/
│   └── cellfocus/
│       ├── __init__.py
│       └── cell_mask_core.py
├── amira/
│   ├── script_objects/
│   │   └── cellfocus.pyscro
│   └── resources/
│       └── cellfocus.rc
├── docs/
│   └── ... (MkDocs pages and images)
└── tests/
    └── test_pipeline.py
```

## Roadmap

- Multipage TIFF stack input (in addition to folders of individual slices)
- Automatic contrast-polarity detection (currently the manual `--cell-dark` switch)
- Performance: reuse the mask computed during scanning instead of recomputing it
- Hosted documentation on Read the Docs and submission to the Journal of Open
  Source Software (JOSS)

## Context

M2 Bioinformatics internship, University of Montpellier. Centre de Biologie
Structurale (CBS), BSME team. Supervisor: Patrick Bron.

## License

MIT, see [LICENSE](LICENSE).
