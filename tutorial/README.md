# Tutorial cases

The tutorials are teaching-oriented entry points built on the public
`tsi_denoising` package. They are kept next to the case data layout and
generated `processed/` products.

## Quick start

From the repository root:

```bash
python -m pip install ".[tutorial]"
python retrieve_datasets.py
jupyter lab tutorial/MARS_DAS/run_example_MARS_DAS.ipynb
```

The notebooks locate their case directory when Jupyter is started from the
repository root, `tutorial/`, or an individual case directory. Downloaded
inputs and generated products are ignored by Git.

## Cases

| Case | Purpose | Public input | Expected scale |
|---|---|---|---|
| [MARS DAS](MARS_DAS/README.md) | Dominant Scholte-wave TSI demonstration | `MARS_DAS/input_public/RR/` | 1,128 pairs, 48 stations |
| [RR Array](RR_Array/README.md) | Four-component Rayleigh-wave separation and TSI | `RR_Array/input_public/{ZZ,ZR,RZ,RR}/` | 253 pairs/component, 23 stations |

See [`data_manifest.yml`](data_manifest.yml) for the machine-readable subset
description.

## Re-running and cache invalidation

Tutorial caches are parameter-dependent. If the frequency band, velocity
window, reference curve, distance threshold, or TSI settings change, remove
only the affected files under the case's `processed/` directory and rerun the
corresponding step.
