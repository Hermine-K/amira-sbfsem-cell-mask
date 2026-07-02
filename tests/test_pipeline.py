"""End-to-end and unit tests for the CellFocus pipeline."""
import os
import numpy as np
import tifffile
import pytest

from cellfocus import run_pipeline, compute_cell_mask, get_bbox


def _make_slice(h=200, w=200, cx=100, cy=100, radius=45, dark=False, seed=0):
    """Synthetic SBF-SEM-like slice: a roundish 'cell' over a noisy background."""
    rng = np.random.default_rng(seed)
    if dark:
        img = np.full((h, w), 200, np.uint8)      # bright background
        cell_val = 40                              # dark cell
    else:
        img = np.full((h, w), 30, np.uint8)        # dark background
        cell_val = 190                             # bright cell
    yy, xx = np.ogrid[:h, :w]
    disk = (yy - cy) ** 2 + (xx - cx) ** 2 <= radius ** 2
    img[disk] = cell_val
    img = np.clip(img.astype(int) + rng.integers(-8, 8, (h, w)), 0, 255).astype(np.uint8)
    return img


def _make_stack(folder, n=6, dark=False):
    os.makedirs(folder, exist_ok=True)
    for i in range(n):
        cx = 100 + int(6 * i)          # cell drifts across slices
        img = _make_slice(cx=cx, cy=100 + i, dark=dark, seed=i)
        tifffile.imwrite(os.path.join(folder, f"slice_{i:03d}.tif"), img)


def test_compute_cell_mask_bright():
    img = _make_slice()
    mask, thr = compute_cell_mask(img, downsample=2, closing_radius=3)
    assert mask.shape == img.shape
    assert mask.sum() > 0
    # the mask centroid should sit near the disk center
    ys, xs = np.where(mask)
    assert abs(xs.mean() - 100) < 20 and abs(ys.mean() - 100) < 20


def test_compute_cell_mask_dark_polarity():
    img = _make_slice(dark=True)
    # without cell_dark it should mostly grab the bright background (wrong target)
    mask_dark, _ = compute_cell_mask(img, downsample=2, closing_radius=3, cell_dark=True)
    assert mask_dark.sum() > 0
    ys, xs = np.where(mask_dark)
    assert abs(xs.mean() - 100) < 20 and abs(ys.mean() - 100) < 20


def test_get_bbox():
    m = np.zeros((50, 50), bool)
    m[10:20, 5:15] = True
    assert get_bbox(m) == (10, 20, 5, 15)
    assert get_bbox(np.zeros((5, 5), bool)) is None


def test_run_pipeline_end_to_end(tmp_path):
    inp = tmp_path / "in"
    out = tmp_path / "out"
    _make_stack(str(inp), n=6)
    res = run_pipeline(str(inp), str(out), downsample=2, closing_radius=3,
                       padding=10, output_scale=1)
    assert res["n_files"] == 6
    # three TIFF series, one file each per slice
    for sub in ("masked", "unmasked", "masks"):
        files = [f for f in os.listdir(out / sub) if f.endswith(".tif")]
        assert len(files) == 6, f"{sub}: {len(files)}"
    # stats + graphs + qc
    assert (out / "stats" / "stats_per_slice.csv").exists()
    assert (out / "stats" / "summary.txt").exists()
    for i in range(1, 6):
        assert (out / "qc" / f"graph_{i}_" ).parent.exists()
    graphs = [f for f in os.listdir(out / "qc") if f.endswith(".png")]
    assert len(graphs) >= 6  # 5 analysis graphs + qc_panel
    assert isinstance(res["reduction_pct"], float)
