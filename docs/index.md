# CellFocus

**Isolate, crop and re-center the main object in large grayscale image stacks.**

CellFocus is a standalone Python tool that isolates a single main object (for
example a cell) from a large, mostly empty background across an image stack,
crops the images around it, re-centers them on a common canvas and reduces the
data size by up to around 90%. It runs three ways: from the command line, as a
Python library, or as an integrated module inside Amira-Avizo 3D.

It was developed for SBF-SEM organelle-segmentation preprocessing at the CBS, but
nothing in the core is specific to Amira or to electron microscopy. Amira is only
one of the ways to run it.

## When it applies

CellFocus fits any grayscale image or stack where **a single dominant object can
be separated from the background by intensity**. Concretely, it assumes:

- **One main object** of interest. It keeps the largest connected region, so it
  is not meant to preserve many scattered objects at once.
- **Enough contrast** for a single global threshold (Otsu) to separate object
  from background, that is a roughly bimodal histogram.
- Either polarity: bright object on dark background (default), or dark object on
  bright background (`--cell-dark`).

The input is a grayscale image, not a binary one: CellFocus does the
thresholding itself. It is an intensity-plus-morphology heuristic, not semantic
segmentation, so it will not fit low-contrast or textured backgrounds, or cases
where several objects must be kept.

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
   background.
2. Morphological closing and hole filling to obtain a full, continuous mask.
3. Largest connected component kept.
4. Per-slice bounding box, then a common canvas to preserve the alignment of the
   3D stack.
5. Re-centering of each object on this canvas.
6. Optional output downsampling.

Developed during an M2 Bioinformatics internship at the Centre de Biologie
Structurale (CBS), BSME team, University of Montpellier.
