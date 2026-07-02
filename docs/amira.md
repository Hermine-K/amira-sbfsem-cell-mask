# Amira integration

Copy the three files to your Amira installation (adapt the root path):

| File | Destination |
|---|---|
| `src/cellfocus/cell_mask_core.py` | `<Amira>/share/python_modules/cellfocus/cell_mask_core.py` |
| `amira/script_objects/cellfocus.pyscro` | `<Amira>/share/python_script_objects/cellfocus.pyscro` |
| `amira/resources/cellfocus.rc` | `<Amira>/share/resources/cellfocus.rc` |

Restart Amira. A **"CellFocus - SBF-SEM Cell Mask"** button appears in the left
panel. Clicking it loads the module and shows its properties.

!!! note
    The path to `cellfocus.pyscro` is hard-coded in `cellfocus.rc`. If your Amira
    installation differs, update this path.

## Performance

For 125 slices at 13664 x 13184 pixels, expect about 45 to 60 minutes of
processing. The Amira UI is frozen during processing (this is normal); progress
is shown in the status field and the console.
