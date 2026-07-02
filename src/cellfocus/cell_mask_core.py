"""
=============================================================================
 CELLFOCUS - cell_mask_core.py
 Cell-focused preprocessing for SBF-SEM segmentation - Core logic
=============================================================================

 This module can be used in TWO ways :

 1. As a library, imported from the Amira .pyscro module :
        from cell_mask_core import run_pipeline
        run_pipeline(input_dir, output_dir, ..., progress_callback=cb)

 2. As a standalone script from the command line :
        python cell_mask_core.py --input ~/data --output ~/out --output-scale 2

 The pipeline segments the cell in each slice (Otsu + morphology), computes
 a per-slice bounding box, determines a common canvas size, and places each
 cell centered in the canvas. Three TIFF series are output (masked, unmasked,
 binary mask), plus a QC panel, 5 analysis graphs, and 3 CSV files.

 AUTHOR : Hermine Kiossou - M2 Bioinformatics - 2026
=============================================================================
"""

import os
import sys
import glob
import time
import warnings

warnings.filterwarnings('ignore', category=FutureWarning)

import numpy as np
import tifffile

from skimage.filters import threshold_otsu
from skimage.morphology import binary_closing, disk
from skimage.measure import label, regionprops
from scipy.ndimage import binary_fill_holes, zoom, center_of_mass

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


# =============================================================================
#  MASK COMPUTATION
# =============================================================================

def compute_cell_mask(image, downsample=4, closing_radius=5,
                      min_area_fraction=0.01, cell_dark=False):
    """
    Compute a binary mask of the main cell in a 2D image.

    Pipeline:
    1. Downsample image for speed.
    2. Otsu thresholding -> binary image (cell vs background).
    3. Morphological closing -> fill small gaps in membranes.
    4. Fill holes -> close internal cavities.
    5. Keep largest connected component -> remove debris.
    6. Upscale mask back to original resolution.

    Parameters
    ----------
    image : np.ndarray (2D)
    downsample : int
        Factor by which to downsample image before computing mask.
    closing_radius : int
        Radius of the disk structuring element for morphological closing.
    min_area_fraction : float
        Minimum area (relative to total) for a component to be kept.

    Returns
    -------
    mask : np.ndarray (2D, bool)
    threshold_value : int
    """
    small = image[::downsample, ::downsample]

    try:
        thr = threshold_otsu(small)
    except ValueError:
        return np.zeros_like(image, dtype=bool), 0

    # Contrast polarity: cell brighter than background (default) or darker.
    binary = (small < thr) if cell_dark else (small > thr)
    binary = binary_closing(binary, disk(closing_radius))
    binary = binary_fill_holes(binary)

    labeled = label(binary)
    if labeled.max() == 0:
        return np.zeros_like(image, dtype=bool), thr

    regions = regionprops(labeled)
    total_area = binary.size
    valid = [r for r in regions if r.area >= min_area_fraction * total_area]
    if not valid:
        return np.zeros_like(image, dtype=bool), thr

    biggest = max(valid, key=lambda r: r.area)
    mask_small = (labeled == biggest.label)

    # Upscale
    zoom_factors = (image.shape[0] / mask_small.shape[0],
                    image.shape[1] / mask_small.shape[1])
    mask_full = zoom(mask_small.astype(np.uint8), zoom=zoom_factors, order=0)
    mask_full = mask_full[:image.shape[0], :image.shape[1]].astype(bool)

    return mask_full, int(thr)


def get_bbox(mask):
    """Tight bounding box. Returns (rmin, rmax, cmin, cmax) or None."""
    if not mask.any():
        return None
    rows = np.any(mask, axis=1)
    cols = np.any(mask, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]
    return int(rmin), int(rmax) + 1, int(cmin), int(cmax) + 1


# =============================================================================
#  PIPELINE STEPS
# =============================================================================

def scan_slices(files, downsample, closing_radius, cell_dark=False,
                progress_callback=None):
    """
    First pass : compute mask, bbox, centroid for each slice.
    """
    records = []
    for i, fp in enumerate(files):
        if progress_callback:
            progress_callback(i, len(files),
                              f"Scan [{i+1}/{len(files)}] {os.path.basename(fp)}")

        img = tifffile.imread(fp)
        if img.ndim != 2:
            continue

        mask, thr = compute_cell_mask(img, downsample, closing_radius,
                                      cell_dark=cell_dark)
        bbox = get_bbox(mask)
        if bbox is None:
            continue

        cy, cx = center_of_mass(mask)
        records.append({
            'filepath': fp,
            'filename': os.path.basename(fp),
            'img_shape': img.shape,
            'dtype': str(img.dtype),
            'bbox': bbox,
            'bbox_h': bbox[1] - bbox[0],
            'bbox_w': bbox[3] - bbox[2],
            'centroid_y': float(cy),
            'centroid_x': float(cx),
            'threshold': thr,
            'cell_pixels': int(mask.sum()),
        })

    return records


def compute_canvas_size(records, padding=100):
    """Choose a canvas that holds every cell once centered, rounded to 16."""
    max_h = max(r['bbox_h'] for r in records)
    max_w = max(r['bbox_w'] for r in records)
    canvas_h = ((max_h + 2*padding + 15) // 16) * 16
    canvas_w = ((max_w + 2*padding + 15) // 16) * 16
    return canvas_h, canvas_w


def process_and_save(record, canvas_h, canvas_w, output_dirs,
                     downsample, closing_radius, output_scale=1, cell_dark=False):
    """
    Re-load slice, compute mask, center cell in canvas, save 3 TIFFs.
    """
    fp = record['filepath']
    img = tifffile.imread(fp)
    mask, _ = compute_cell_mask(img, downsample, closing_radius, cell_dark=cell_dark)

    cy, cx = center_of_mass(mask)
    canvas_cy, canvas_cx = canvas_h // 2, canvas_w // 2
    dy = canvas_cy - int(round(cy))
    dx = canvas_cx - int(round(cx))

    canvas_masked = np.zeros((canvas_h, canvas_w), dtype=img.dtype)
    canvas_unmasked = np.zeros((canvas_h, canvas_w), dtype=img.dtype)
    canvas_mask = np.zeros((canvas_h, canvas_w), dtype=np.uint8)

    H, W = img.shape
    y0_dst = max(0, dy);      x0_dst = max(0, dx)
    y1_dst = min(canvas_h, dy + H); x1_dst = min(canvas_w, dx + W)
    y0_src = max(0, -dy);     x0_src = max(0, -dx)
    y1_src = y0_src + (y1_dst - y0_dst)
    x1_src = x0_src + (x1_dst - x0_dst)

    img_region = img[y0_src:y1_src, x0_src:x1_src]
    mask_region = mask[y0_src:y1_src, x0_src:x1_src]

    canvas_unmasked[y0_dst:y1_dst, x0_dst:x1_dst] = img_region
    masked_region = img_region.copy()
    masked_region[~mask_region] = 0
    canvas_masked[y0_dst:y1_dst, x0_dst:x1_dst] = masked_region
    canvas_mask[y0_dst:y1_dst, x0_dst:x1_dst] = mask_region.astype(np.uint8) * 255

    # Optional output downsampling
    if output_scale > 1:
        canvas_masked = canvas_masked[::output_scale, ::output_scale]
        canvas_unmasked = canvas_unmasked[::output_scale, ::output_scale]
        canvas_mask = canvas_mask[::output_scale, ::output_scale]

    base = os.path.splitext(record['filename'])[0]
    tifffile.imwrite(os.path.join(output_dirs['masked'], f"{base}.tif"),
                     canvas_masked, compression='zlib')
    tifffile.imwrite(os.path.join(output_dirs['unmasked'], f"{base}.tif"),
                     canvas_unmasked, compression='zlib')
    tifffile.imwrite(os.path.join(output_dirs['masks'], f"{base}.tif"),
                     canvas_mask, compression='zlib')

    cell_vals = img_region[mask_region]
    mean_cell = float(cell_vals.mean()) if cell_vals.size else float('nan')
    std_cell = float(cell_vals.std()) if cell_vals.size else float('nan')

    size_before = os.path.getsize(fp) / (1024*1024)
    size_after = os.path.getsize(
        os.path.join(output_dirs['masked'], f"{base}.tif")
    ) / (1024*1024)

    return {
        'dy': dy, 'dx': dx,
        'size_mo_before': round(size_before, 2),
        'size_mo_after_masked': round(size_after, 2),
        'reduction_pct': round((1 - size_after / size_before) * 100, 1) if size_before > 0 else 0,
        'mean_cell': round(mean_cell, 2),
        'std_cell': round(std_cell, 2),
    }


# =============================================================================
#  ANALYSIS GRAPHS
# =============================================================================

# Colorblind-safe palette (Okabe-Ito)
CB_BLUE   = '#0072B2'   # retained series 1
CB_ORANGE = '#E69F00'   # retained series 2
CB_GRAY   = '#999999'   # context / non-retained
CB_BLACK  = '#000000'


def _thin_ticks(ax, positions, labels, max_ticks=12):
    """Show at most ~max_ticks evenly spaced x-ticks, so the axis stays
    readable whatever the number of slices."""
    positions = np.asarray(positions)
    labels = np.asarray(labels)
    if len(positions) > max_ticks:
        idx = np.linspace(0, len(positions) - 1, max_ticks).astype(int)
        positions = positions[idx]
        labels = labels[idx]
    ax.set_xticks(positions)
    ax.set_xticklabels([str(int(l)) for l in labels], fontsize=9)


def generate_analysis_graphs(df, output_dir, canvas_h, canvas_w):
    """Generate 5 analysis graphs for the report (colorblind-safe, readable axes)."""
    n = len(df)
    x = np.arange(1, n + 1)          # slice number, 1-based
    xp = np.arange(n)                # bar positions, 0-based

    # 1. Bounding box per slice
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x, df['bbox_h'], 'o-', label='Hauteur (bbox_h)', color=CB_BLUE, markersize=4)
    ax.plot(x, df['bbox_w'], 's-', label='Largeur (bbox_w)', color=CB_ORANGE, markersize=4)
    ax.set_xlabel('Numero de slice'); ax.set_ylabel('Dimensions (px)')
    ax.set_title('Taille de la bounding box cellulaire par slice')
    ax.legend(); ax.grid(True, alpha=0.3)
    _thin_ticks(ax, x, x)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'graph_1_bbox_per_slice.png'), dpi=120)
    plt.close()

    # 2. Size reduction
    fig, ax = plt.subplots(figsize=(10, 5))
    width = 0.42
    ax.bar(xp - width/2, df['size_mo_before'], width, label='Avant (original)', color=CB_GRAY)
    ax.bar(xp + width/2, df['size_mo_after_masked'], width, label='Apres (masque)', color=CB_BLUE)
    red = (1 - df['size_mo_after_masked'].sum() / df['size_mo_before'].sum()) * 100
    ax.set_xlabel('Numero de slice'); ax.set_ylabel('Taille fichier (Mo)')
    ax.set_title(f'Reduction de taille par slice (moyenne : -{red:.1f}%)')
    ax.legend(); ax.grid(True, alpha=0.3, axis='y')
    _thin_ticks(ax, xp, xp + 1)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'graph_2_size_reduction.png'), dpi=120)
    plt.close()

    # 3. Centroid trajectory
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(df['centroid_x'], df['centroid_y'], '-', color=CB_BLUE, linewidth=1, alpha=0.7)
    ax.scatter(df['centroid_x'], df['centroid_y'], color=CB_BLUE, s=30, zorder=3)
    ax.scatter(df['centroid_x'].iloc[0], df['centroid_y'].iloc[0],
               color=CB_BLUE, s=120, edgecolors=CB_BLACK, linewidths=1.5, zorder=4)
    ax.scatter(df['centroid_x'].iloc[-1], df['centroid_y'].iloc[-1],
               color=CB_ORANGE, s=120, edgecolors=CB_BLACK, linewidths=1.5, zorder=4)
    ax.annotate('Debut (slice 1)', (df['centroid_x'].iloc[0], df['centroid_y'].iloc[0]),
                xytext=(10, 10), textcoords='offset points', color=CB_BLACK, fontweight='bold')
    ax.annotate(f'Fin (slice {n})', (df['centroid_x'].iloc[-1], df['centroid_y'].iloc[-1]),
                xytext=(10, 10), textcoords='offset points', color=CB_BLACK, fontweight='bold')
    ax.set_xlabel('Centroide X (px)'); ax.set_ylabel('Centroide Y (px)')
    ax.set_title("Derive spatiale de la cellule entre slices\n(justifie le recalage)")
    ax.grid(True, alpha=0.3); ax.invert_yaxis(); ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'graph_3_centroid_drift.png'), dpi=120)
    plt.close()

    # 4. Intensity stability
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.errorbar(x, df['mean_cell'], yerr=df['std_cell'],
                fmt='o-', capsize=4, color=CB_BLUE, ecolor=CB_GRAY, markersize=4,
                label='Moyenne +/- ecart-type')
    ax.set_xlabel('Numero de slice'); ax.set_ylabel('Intensite (niveaux de gris)')
    ax.set_title("Distribution des intensites a l'interieur de la cellule\n(preuve de la preservation du signal)")
    ax.legend(); ax.grid(True, alpha=0.3)
    _thin_ticks(ax, x, x)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'graph_4_intensity_stability.png'), dpi=120)
    plt.close()

    # 5. Outlier detection
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    all_sizes = np.concatenate([df['bbox_h'].values, df['bbox_w'].values])
    q1, q3 = np.percentile(all_sizes, [25, 75])
    upper_fence = q3 + 1.5 * (q3 - q1)

    is_out_h = df['bbox_h'] > upper_fence
    is_out_w = df['bbox_w'] > upper_fence
    ax1.scatter(x, df['bbox_h'], color=CB_BLUE, s=45, label='bbox_h')
    ax1.scatter(x, df['bbox_w'], color=CB_ORANGE, s=45, marker='s', label='bbox_w')
    if is_out_h.any():
        ax1.scatter(x[is_out_h.values], df['bbox_h'][is_out_h.values],
                    facecolors='none', edgecolors=CB_BLACK, linewidths=2, s=140,
                    label='Outlier')
    if is_out_w.any():
        ax1.scatter(x[is_out_w.values], df['bbox_w'][is_out_w.values],
                    facecolors='none', edgecolors=CB_BLACK, linewidths=2, s=140)
    ax1.axhline(upper_fence, color=CB_BLACK, linestyle='--', alpha=0.6,
                label=f'Seuil outlier ({upper_fence:.0f} px)')
    ax1.set_xlabel('Numero de slice'); ax1.set_ylabel('Dimension (px)')
    ax1.set_title(f'Detection des outliers (IQR)\nCercles noirs : slices > {upper_fence:.0f} px')
    ax1.legend(); ax1.grid(True, alpha=0.3)
    _thin_ticks(ax1, x, x)

    normal = df[(df['bbox_h'] <= upper_fence) & (df['bbox_w'] <= upper_fence)]
    if len(normal) > 0:
        canvas_all = df['bbox_h'].max() * df['bbox_w'].max() / 1e6
        canvas_nor = normal['bbox_h'].max() * normal['bbox_w'].max() / 1e6
        red = (1 - canvas_nor / canvas_all) * 100
        ax2.bar(['Canvas actuel\n(avec outliers)', 'Canvas possible\n(sans outliers)'],
                [canvas_all, canvas_nor], color=[CB_GRAY, CB_BLUE])
        ax2.set_ylabel('Taille canvas (Mpx)')
        ax2.set_title(f'Impact des outliers\nReduction possible : -{red:.0f}%')
        for i, v in enumerate([canvas_all, canvas_nor]):
            ax2.text(i, v + max(canvas_all, canvas_nor) * 0.01, f'{v:.1f} Mpx',
                     ha='center', fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'graph_5_outlier_analysis.png'), dpi=120)
    plt.close()


# =============================================================================
#  QC PANEL
# =============================================================================

def generate_qc_panel(records, output_dirs, canvas_h, canvas_w,
                      out_path, n_thumb=10, thumb_size=400, cell_dark=False):
    """Visual QC : [original | mask overlay | final masked] for N slices."""
    n = min(n_thumb, len(records))
    if n == 0:
        return
    indices = np.linspace(0, len(records) - 1, n).astype(int)

    fig, axes = plt.subplots(n, 3, figsize=(12, 3 * n))
    if n == 1:
        axes = axes.reshape(1, 3)

    for row, i in enumerate(indices):
        r = records[i]
        img = tifffile.imread(r['filepath'])
        mask, _ = compute_cell_mask(img, downsample=4, closing_radius=5,
                                    cell_dark=cell_dark)

        def to_display(arr):
            if arr.dtype != np.uint8:
                lo, hi = np.percentile(arr, [1, 99])
                if hi > lo:
                    arr = np.clip((arr - lo) * 255 / (hi - lo), 0, 255).astype(np.uint8)
                else:
                    arr = arr.astype(np.uint8)
            return arr

        ts = thumb_size
        img_t = to_display(img)[::max(1, img.shape[0]//ts), ::max(1, img.shape[1]//ts)]
        mask_t = mask[::max(1, mask.shape[0]//ts), ::max(1, mask.shape[1]//ts)]

        masked_fp = os.path.join(output_dirs['masked'], r['filename'])
        if os.path.exists(masked_fp):
            masked_full = tifffile.imread(masked_fp)
            masked_t = to_display(masked_full)[::max(1, masked_full.shape[0]//ts),
                                                 ::max(1, masked_full.shape[1]//ts)]
        else:
            masked_t = np.zeros_like(img_t)

        axes[row, 0].imshow(img_t, cmap='gray')
        axes[row, 0].set_title(
            f"{r['filename']}\nOriginal {r['img_shape'][0]}x{r['img_shape'][1]}",
            fontsize=9)
        axes[row, 0].axis('off')

        overlay = np.stack([img_t, img_t, img_t], axis=-1).astype(float) / 255
        overlay[mask_t > 0] = overlay[mask_t > 0] * 0.6 + np.array([0.0, 0.5, 0.0]) * 0.4
        axes[row, 1].imshow(overlay)
        axes[row, 1].set_title(f"Mask overlay\nCell = {100*mask_t.mean():.1f}%",
                               fontsize=9)
        axes[row, 1].axis('off')

        axes[row, 2].imshow(masked_t, cmap='gray')
        axes[row, 2].set_title(f"Final (centered)\n{masked_t.shape[0]*max(1,img.shape[0]//ts)}x"
                               f"{masked_t.shape[1]*max(1,img.shape[1]//ts)}", fontsize=9)
        axes[row, 2].axis('off')

    plt.suptitle("QC Panel - Cell Masking v2 Pipeline", fontsize=14, y=1.0)
    plt.tight_layout()
    plt.savefig(out_path, dpi=100, bbox_inches='tight')
    plt.close()


# =============================================================================
#  HIGH-LEVEL ENTRY POINT (used by both CLI and Amira module)
# =============================================================================

def run_pipeline(input_dir, output_dir, n_samples=0,
                 downsample=4, closing_radius=5, padding=100, output_scale=1,
                 cell_dark=False, progress_callback=None):
    """
    Run the full pipeline on a folder of TIFF files.

    Parameters
    ----------
    input_dir : str
        Folder containing .tif files.
    output_dir : str
        Root output folder (subfolders masked/unmasked/masks/stats/qc are created).
    n_samples : int
        If 0, process all files. If > 0, sample N evenly-spaced slices (test mode).
    downsample : int
        Downsampling factor for mask computation.
    closing_radius : int
        Radius for morphological closing.
    padding : int
        Extra padding around max bbox in canvas.
    output_scale : int
        Divide final canvas by this factor (1 = full resolution).
    progress_callback : callable or None
        Function (current, total, message) called during processing.
        Use to update GUI progress.

    Returns
    -------
    dict with summary :
        {'n_files': int, 'canvas_h': int, 'canvas_w': int,
         'size_mo_before': float, 'size_mo_after': float, 'reduction_pct': float,
         'output_dir': str, 'masked_dir': str}
    """
    output_dirs = {
        'masked':   os.path.join(output_dir, 'masked'),
        'unmasked': os.path.join(output_dir, 'unmasked'),
        'masks':    os.path.join(output_dir, 'masks'),
        'qc':       os.path.join(output_dir, 'qc'),
        'stats':    os.path.join(output_dir, 'stats'),
    }
    for d in output_dirs.values():
        os.makedirs(d, exist_ok=True)

    # Find files
    files = sorted(glob.glob(os.path.join(input_dir, "*.tif")) +
                   glob.glob(os.path.join(input_dir, "*.tiff")))
    if not files:
        raise FileNotFoundError(f"No .tif files in {input_dir}")

    if n_samples > 0 and n_samples < len(files):
        idx = np.linspace(0, len(files) - 1, n_samples).astype(int)
        files = [files[i] for i in idx]

    total_steps = 2 * len(files) + 3  # Scan + process + 3 final steps

    if progress_callback:
        progress_callback(0, total_steps, f"Scanning {len(files)} files ...")

    # Pass 1 : scan
    records = scan_slices(files, downsample, closing_radius, cell_dark=cell_dark,
                          progress_callback=lambda c, t, m:
                              progress_callback(c, total_steps, m) if progress_callback else None)
    if not records:
        raise ValueError("No valid slices found (all masks empty ?)")

    # Pass 2 : canvas
    canvas_h, canvas_w = compute_canvas_size(records, padding=padding)
    if progress_callback:
        progress_callback(len(files), total_steps,
                          f"Canvas size : {canvas_h} x {canvas_w} px")

    # Pass 3 : process and save
    all_stats = []
    for i, r in enumerate(records):
        if progress_callback:
            progress_callback(len(files) + i + 1, total_steps,
                              f"Saving [{i+1}/{len(records)}] {r['filename']}")
        result = process_and_save(r, canvas_h, canvas_w, output_dirs,
                                  downsample, closing_radius, output_scale,
                                  cell_dark=cell_dark)
        stats = {**r, **result}
        stats.pop('filepath')
        stats.pop('img_shape')
        stats['bbox'] = str(stats['bbox'])
        all_stats.append(stats)

    # Stats CSV
    if progress_callback:
        progress_callback(2*len(files) + 1, total_steps, "Writing stats CSV ...")
    df = pd.DataFrame(all_stats)
    df.to_csv(os.path.join(output_dirs['stats'], 'stats_per_slice.csv'), index=False)

    df_trans = df[['filename', 'bbox', 'centroid_y', 'centroid_x', 'dy', 'dx']].copy()
    df_trans['canvas_h'] = canvas_h
    df_trans['canvas_w'] = canvas_w
    df_trans.to_csv(os.path.join(output_dirs['stats'], 'transformations.csv'), index=False)

    # Summary
    total_before = df['size_mo_before'].sum()
    total_after = df['size_mo_after_masked'].sum()
    reduction = (1 - total_after / total_before) * 100 if total_before > 0 else 0

    with open(os.path.join(output_dirs['stats'], 'summary.txt'), 'w', encoding='utf-8') as f:
        f.write("SBF-SEM Preprocessing v2 - Run summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Input folder     : {input_dir}\n")
        f.write(f"Output folder    : {output_dir}\n")
        f.write(f"Files processed  : {len(records)}\n\n")
        f.write(f"Canvas size      : {canvas_h} x {canvas_w} px\n")
        f.write(f"Output scale     : {output_scale} (final = canvas / scale)\n")
        f.write(f"Canvas padding   : {padding} px\n\n")
        f.write(f"Size before total : {total_before:.1f} Mo ({total_before/1024:.2f} Go)\n")
        f.write(f"Size after total  : {total_after:.1f} Mo ({total_after/1024:.2f} Go)\n")
        f.write(f"Reduction         : {reduction:.1f}%\n\n")
        f.write(f"Mean bbox per slice : "
                f"{df['bbox_h'].mean():.0f} x {df['bbox_w'].mean():.0f} px\n")
        f.write(f"Mean cell intensity : {df['mean_cell'].mean():.1f} +/- "
                f"{df['std_cell'].mean():.1f}\n")

    # QC panel
    if progress_callback:
        progress_callback(2*len(files) + 2, total_steps, "Generating QC panel ...")
    generate_qc_panel(records, output_dirs, canvas_h, canvas_w,
                      os.path.join(output_dirs['qc'], 'qc_panel.png'),
                      cell_dark=cell_dark)

    # Analysis graphs
    if progress_callback:
        progress_callback(2*len(files) + 3, total_steps, "Generating analysis graphs ...")
    generate_analysis_graphs(df, output_dirs['qc'], canvas_h, canvas_w)

    if progress_callback:
        progress_callback(total_steps, total_steps, "Done.")

    return {
        'n_files': len(records),
        'canvas_h': canvas_h,
        'canvas_w': canvas_w,
        'output_scale': output_scale,
        'size_mo_before': total_before,
        'size_mo_after': total_after,
        'reduction_pct': reduction,
        'output_dir': output_dir,
        'masked_dir': output_dirs['masked'],
    }


# =============================================================================
#  COMMAND-LINE ENTRY POINT
# =============================================================================

def _cli_progress(current, total, message):
    print(f"  [{current}/{total}] {message}")


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="SBF-SEM Cell Masking Pipeline v2")
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--n-samples', type=int, default=0)
    parser.add_argument('--downsample', type=int, default=4)
    parser.add_argument('--closing-radius', type=int, default=5)
    parser.add_argument('--padding', type=int, default=100)
    parser.add_argument('--output-scale', type=int, default=1)
    parser.add_argument('--cell-dark', action='store_true',
                        help='Cell is darker than background (invert threshold).')
    args = parser.parse_args()

    t0 = time.time()
    result = run_pipeline(
        input_dir=args.input,
        output_dir=args.output,
        n_samples=args.n_samples,
        downsample=args.downsample,
        closing_radius=args.closing_radius,
        padding=args.padding,
        output_scale=args.output_scale,
        cell_dark=args.cell_dark,
        progress_callback=_cli_progress,
    )
    print(f"\nDone in {time.time()-t0:.1f}s")
    print(f"  Files        : {result['n_files']}")
    print(f"  Canvas       : {result['canvas_h']} x {result['canvas_w']} px")
    print(f"  Reduction    : {result['reduction_pct']:.1f}%")
    print(f"  Output       : {result['output_dir']}")


if __name__ == '__main__':
    main()
