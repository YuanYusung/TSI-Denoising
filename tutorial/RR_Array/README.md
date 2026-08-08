# RR Array tutorial

This case demonstrates the four-component Rayleigh-wave workflow:

1. read and preprocess `ZZ`, `ZR`, `RZ`, and `RR` ANC data;
2. compare component MASW diagnostics;
3. separate retrograde and prograde candidates;
4. apply reference-curve phase matching;
5. denoise the selected M0 and M1 candidates independently.

## Run

From the repository root:

```bash
python -m pip install ".[tutorial]"
python retrieve_datasets.py
jupyter lab tutorial/RR_Array/run_example_RR_Array.ipynb
```

The deterministic reference curves in the notebook are the default path and
do not require a GUI. An optional cell documents how to pick replacement
curves with `%matplotlib qt`; it should be run only when a local interactive
desktop is available.

Each component currently contains 253 SAC pairs from the same 23-station
public subset, with 4,001 samples at 50 Hz. The four components must remain
geometrically aligned before polarization separation.

## Interpretation boundary

Polarization separation produces candidate wavefields; it is not by itself a
proof of modal separation. Inspect dispersion continuity and spatial
stability before interpreting the M0/M1 TSI outputs.
