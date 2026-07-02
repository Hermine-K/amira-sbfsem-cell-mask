"""CellFocus: cell-focused preprocessing for SBF-SEM organelle segmentation.

Automatic per-slice cell masking, common-canvas re-centering, background
removal and lossless data reduction, usable from the command line or inside
Amira-Avizo 3D.
"""
from .cell_mask_core import run_pipeline, compute_cell_mask, get_bbox, __version__

__all__ = ["run_pipeline", "compute_cell_mask", "get_bbox", "__version__"]
