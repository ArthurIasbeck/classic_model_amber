# Repository Guide

## Scope and Layout

- This is a flat Python/MATLAB research workspace, not an installable package: there is no `pyproject.toml`, package `__init__.py`, CI, or task runner.
- `src/dataset.py` loads and plots 13-column experiment data; `src/experimental_frequency_response.py` contains the class-based direct-DFT experiment; `tests/00_compute_exp_freq_resp.py` is a separate visual comparison script, not an assertion-based test.
- The MATLAB scripts in `src/utils/` export columns as time, four inputs, four disturbances, then four outputs, matching `Dataset.load()`.
- Check `git status` before cleanup or broad edits. At the time this guide was written, only `README.md` was tracked; the implementation, tests, data, plots, requirements, and this guide were untracked.

## Commands

- Install the pinned environment: `python -m pip install -r requirements.txt`.
- Check formatting: `python -m black --check src tests`; apply it with `python -m black src tests`. No separate lint or typecheck command is configured.
- Run the dataset example from `src`, because its input path is working-directory-relative: `(cd src && MPLBACKEND=Agg python dataset.py)`.
- Run the class-based frequency-response example from `src`, because its output directory is working-directory-relative: `(cd src && MPLBACKEND=Agg python experimental_frequency_response.py)`.
- The only test-like check is `MPLBACKEND=Agg python tests/00_compute_exp_freq_resp.py` from the repository root. It generates plots but has no assertions and is not pytest-discoverable; do not report it as an automated test pass.

## Data and Verification Hazards

- `Dataset.load()` prefers a sibling `.npz` even when passed a `.txt` path. If no `.npz` exists, it loads the entire text file and writes an uncompressed `.npz` beside it.
- Data files are very large (`chirp.txt` is about 2.6 GB and `chirp.npz` about 1.5 GB). Avoid regenerating caches, copying data, or using the full dataset for routine verification.
- Both frequency-response implementations construct a dense direct-DFT kernel. Their 20,000-sample examples allocate roughly 200 million complex values (about 3.2 GB for the kernel alone), so do not run them as a routine smoke test on a memory-constrained machine.
- Python entrypoints call `plt.show()`; use `MPLBACKEND=Agg` for non-interactive runs. They overwrite SVGs in `plots/`.
- In `ExperimentalFrequencyResponse`, `angular_frequencies` and `max_freq` are radians/second even though the plot x-axis says Hz. Preserve units deliberately when changing this code.
