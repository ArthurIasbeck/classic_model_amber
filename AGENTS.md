# Repository Guide

## Scope and Layout

- This is a flat Python/MATLAB research workspace, not an installable package: there is no `pyproject.toml`, package `__init__.py`, CI, or task runner.
- `src/dataset.py` loads and plots 13-column experiment data, computes correlations, and plots FFTs; `src/experimental_frequency_response.py` contains the class-based MIMO frequency-response experiment; `src/experimental_frequency_response_siso.py` contains filtering, clipping, averaging, and SISO frequency-response analysis.
- `src/experimental_open_loop_frequency_response_siso.py` computes SISO open-loop frequency response from saved closed-loop data; `src/ol_freq_resp_siso.py` runs the `V13` and `W13` analyses and generates controller, open-loop, and delay-removal plots.
- `src/freq_resp.py` and `src/freq_resp_siso.py` are data-analysis entrypoints for real experiments; `tests/00_compute_exp_freq_resp.py` and `tests/01_read_csv.py` are scripts, not assertion-based automated tests.
- `src/build_mimo_model_from_frequency_response.py` fits a reduced ROSS/control MIMO model to saved open-loop responses; `src/utils/amber_models.py` builds rotor models, compares meshes, computes numerical FRFs, and compares them with experimental data.
- `src/utils.py` contains interpolation helpers. The `src/utils/` directory contains MATLAB conversion scripts; its `txt_data/` subdirectory is ignored because it contains multi-gigabyte raw exports.
- The MATLAB scripts in `src/utils/` export columns as time, four inputs, four disturbances, then four outputs, matching `Dataset.load()`.
- Check `git status` before cleanup or broad edits. Data files are local and ignored; generated plots and selected research artifacts are versioned when they document an analysis result.

## Commands

- Install the pinned environment: `python -m pip install -r requirements.txt`.
- Check formatting: `python -m black --check src tests`; apply it with `python -m black src tests`. No separate lint or typecheck command is configured.
- Run the dataset example from `src`, because its input path is working-directory-relative: `(cd src && MPLBACKEND=Agg python dataset.py)`.
- Run the class-based frequency-response example from `src`, because its output directory is working-directory-relative: `(cd src && MPLBACKEND=Agg python experimental_frequency_response.py)`.
- Run the real-data MIMO analysis from `src` only when the required files exist in `data/`: `(cd src && MPLBACKEND=Agg python freq_resp.py)`.
- Run the real-data SISO analysis from `src` only when the required files exist in `data/`: `(cd src && MPLBACKEND=Agg python freq_resp_siso.py)`.
- Run both labeled SISO open-loop analyses from `src` when the corresponding closed-loop frequency-response files exist in `data/`: `(cd src && MPLBACKEND=Agg python ol_freq_resp_siso.py)`.
- Run the ROSS rotor/model workflows only when their local experimental FRF files and optional model dependencies are available: `(cd src && MPLBACKEND=Agg python build_mimo_model_from_frequency_response.py)`.
- The only test-like check is `MPLBACKEND=Agg python tests/00_compute_exp_freq_resp.py` from the repository root. It generates plots but has no assertions and is not pytest-discoverable; do not report it as an automated test pass.
- The CSV reader can be exercised with `python tests/01_read_csv.py`, but this is also a script rather than an automated test.

## Data and Verification Hazards

- `Dataset.load()` prefers a sibling `.npz` even when passed a `.txt` path. If no `.npz` exists, it loads the entire text file and writes an uncompressed `.npz` beside it.
- Data files are very large (`chirp.txt` is about 2.6 GB and `chirp.npz` about 1.5 GB). Avoid regenerating caches, copying data, or using the full dataset for routine verification.
- Both frequency-response implementations construct a dense direct-DFT kernel. Their 20,000-sample examples allocate roughly 200 million complex values (about 3.2 GB for the kernel alone), so do not run them as a routine smoke test on a memory-constrained machine.
- Python entrypoints call `plt.show()`; use `MPLBACKEND=Agg` for non-interactive runs. They overwrite SVGs in `plots/`.
- Plot interpolation expands signals to 1,000,000 points, so plotting can use substantial memory even when the input data is smaller.
- `src/freq_resp.py` and `src/freq_resp_siso.py` require local experimental files under `data/`; those files are ignored and are not reproducible from the repository alone.
- In `ExperimentalFrequencyResponse`, `angular_frequencies` and `max_freq` are radians/second even though the plot x-axis says Hz. Preserve units deliberately when changing this code.
- The rotor-model workflows use ROSS, `control`, and Plotly in addition to the analysis dependencies; keep their relative data paths and configurable mesh/suspension settings documented when changing them.
