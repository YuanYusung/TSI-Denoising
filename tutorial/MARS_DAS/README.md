# MARS DAS tutorial

This case demonstrates a single-component Scholte-wave workflow:

1. read and preprocess the ANC SAC pairs;
2. inspect MASW energy in the selected frequency band;
3. inspect one representative three-station-interferometry result;
4. run iterative denoising and compare the original and final wavefields.

## Run

From the repository root:

```bash
python -m pip install ".[tutorial]"
python retrieve_datasets.py
jupyter lab tutorial/MARS_DAS/run_example_MARS_DAS.ipynb
```

The notebook is non-interactive by default. It reads `input_public/RR/` and
writes caches such as `processed/rr_wavefield.npz`,
`processed/rr_masw.npz`, and the paired `Denoised_*` result files.

The public subset currently contains 1,128 SAC pairs from 48 stations, with
2,501 samples at 25 Hz. These values document the teaching layout, not a
scientific quality guarantee.

## Interpretation boundary

The case is a teaching example for a selected approximately linear cable
segment. It does not claim that the full cable is one-dimensional or that all
Scholte-wave modes have been separated.
