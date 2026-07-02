# Installation

## As a Python package

From the repository root:

```bash
pip install -e .
```

This exposes both the importable API and the `cellfocus` command-line entry
point.

## Dependencies only

```bash
pip install -r requirements.txt
```

Requirements: Python 3.9+, numpy, tifffile, scikit-image, scipy, matplotlib,
pandas.

## Isolated environment (recommended)

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
pip install -e .
```
