# Agent Guide

## Repository shape

- This is a research-analysis workspace, not an installable Python package. There is no `pyproject.toml`, CI workflow, task runner, or automated test suite.
- `src/dataset.py` owns TXT/CSV/NPZ loading and plotting; `src/experimental_frequency_response.py` and `src/experimental_frequency_response_siso.py` implement MIMO/SISO frequency-response analysis; `src/experimental_open_loop_frequency_response_siso.py` derives open-loop SISO results from saved closed-loop data.
- `src/freq_resp_siso.py` and `src/ol_freq_resp_siso.py` are real-data entrypoints. `src/build_mimo_model_from_frequency_response.py` is the ROSS rotor/model workflow and requires local frequency-response data.

## Setup and commands

- Install the pinned environment with `python -m pip install -r requirements.txt`.
- Run scripts from `src/`, not the repository root: their bare imports and default `../data`/`../plots` paths assume that working directory. Use `MPLBACKEND=Agg` for headless execution, for example `cd src && MPLBACKEND=Agg python ol_freq_resp_siso.py`.
- Useful focused commands are `cd src && MPLBACKEND=Agg python dataset.py`, `cd src && MPLBACKEND=Agg python experimental_frequency_response.py`, `cd src && MPLBACKEND=Agg python freq_resp_siso.py`, and `cd src && MPLBACKEND=Agg python build_mimo_model_from_frequency_response.py` when their required data exists.
- Format Python with `python -m black src main.py`; check with `python -m black --check src main.py`. Check syntax without running analyses with `python -m compileall src main.py`.

## Data and analysis hazards

- `data/` is ignored and must be supplied locally. `Dataset.load()` stores signals as `(n_signals, n_samples)` and, for TXT/CSV input, creates a sibling `.npz` cache; do not regenerate caches or copy raw data casually.
- Analysis scripts write SVG/PNG results into `plots/`, and most entrypoints call `plt.show()`. Treat generated plot changes as intentional research artifacts and use `MPLBACKEND=Agg` in automation.
- Frequency-response examples use dense direct-DFT matrices and can require multiple gigabytes of RAM; do not use the full-data workflows as routine smoke tests.
- `angular_frequencies` and frequency-response `max_freq` values are in radians/second, even where plot labels may say Hz. Preserve units when changing analysis code.
- There are no current assertion-based tests. A successful Black check, compile check, or script run is not a test-suite result; real-data scripts also depend on untracked files under `data/`.
