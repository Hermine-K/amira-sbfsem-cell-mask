# CellFocus

**Cell-focused preprocessing for SBF-SEM organelle segmentation.**

CellFocus isolates the cell on each slice of an SBF-SEM volume, re-centers it on
a common canvas, removes the resin background and reduces the data size by around
90% while preserving the biological structures. It runs both as a command-line
tool and as an integrated module inside Amira-Avizo 3D 2025.1.

## Why this tool

Raw SBF-SEM volumes are large (> 20 GB) and made up of 60 to 95% empty
background. During segmentation training, random patch sampling then falls mostly
on this background and the model learns almost nothing useful. CellFocus solves
the problem upstream: it segments the cell, re-centers it on a common canvas and
removes the useless border, which greatly increases the proportion of informative
pixels.

## What the pipeline does

For each slice of the volume:

1. Otsu thresholding on a downsampled version to separate the cell from the
   background. Both contrast polarities are supported: bright cell on dark
   background (default) or dark cell on bright background.
2. Morphological closing and hole filling to obtain a full, continuous mask.
3. Largest connected component kept (the cell of interest).
4. Per-slice bounding box, then a common canvas to preserve the alignment of the
   3D stack.
5. Re-centering of each cell on this canvas.
6. Optional output downsampling.

Developed during an M2 Bioinformatics internship at the Centre de Biologie
Structurale (CBS), BSME team, University of Montpellier.
