# Usage

## Command line

```bash
cellfocus --input /path/to/tiffs --output /path/to/out --n-samples 5
```

Without installing the entry point:

```bash
python src/cellfocus/cell_mask_core.py --input /path/to/tiffs --output /path/to/out --n-samples 5
```

!!! tip "Test quickly"
    `--n-samples 5` processes only 5 evenly spaced slices, a few seconds instead
    of the full volume. Remove it for a complete run.

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

## Python API

```python
from cellfocus import run_pipeline

result = run_pipeline("data/denoised", "data/out", output_scale=2)
print(result["reduction_pct"])
```
