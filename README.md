# riskcontroltriage

This repository implements a training and evaluation workflow to study whether **short molecular dynamics (MD) windows** carry enough mechanistic signal to support **early termination decisions** under an explicitly controlled **false positive rate (FPR) budget**. Rather than treating early MD analysis as a purely predictive task, the pipeline converts model outputs into an operational screening policy calibrated to a user-specified risk threshold.

---

## Core idea

For each protein–ligand pocket system, a full 20 ns MD trajectory serves as the reference simulation. Early windows (e.g., 2 ns, 3 ns, 5 ns) are converted into mechanistic feature vectors and used to predict whether the pocket will remain stable at the long-horizon endpoint.

A calibrated threshold is selected to satisfy a user-defined FPR budget, turning probabilistic predictions into a practical decision rule:

- **Terminate early** for complexes predicted to be unstable
- **Continue simulation** for complexes not confidently classified as unstable

This framing is suited for real screening workflows where risk must be controlled explicitly.

---

## Project layout

```
riskcontroltriage/
├── README.md
└── pipeline/
    └── training/
        ├── compute_ground_truth_labels.py
        ├── generate_figure_both_cv_ligand_update4_patched.py
        ├── generate_fingerprint_summary.py
        ├── generate_timescale_fingerprints.py
        ├── label_cutoff.py
        └── outputs/
```

---

## Main scripts

### `generate_timescale_fingerprints.py`
Builds mechanistic feature representations from early MD trajectory windows. Reads trajectory-derived measurements, aggregates window-specific descriptors, and exports feature matrices indexed by protein / pocket / ligand / horizon.

### `compute_ground_truth_labels.py`
Computes endpoint labels from full-horizon MD trajectories. Reads long-horizon stability metrics, applies stability criteria, and generates binary labels for supervised training.

### `label_cutoff.py`
Centralizes cutoff logic for converting continuous endpoint metrics into stable/unstable labels. Ensures reproducible and consistent label generation across datasets, ligands, and runs.

### `generate_fingerprint_summary.py`
Produces summary statistics and reporting tables for the mechanistic fingerprints, including feature coverage and grouped component summaries.

### `generate_figure_both_cv_ligand_update4_patched.py`
Generates the main cross-validation and ligand-transfer figures used in the manuscript, including horizon-wise performance curves and ligand-level comparison panels.

---

## Reproducibility workflow

```bash
# 1. Generate endpoint labels
python pipeline/training/compute_ground_truth_labels.py

# 2. Build early-window fingerprints
python pipeline/training/generate_timescale_fingerprints.py

# 3. Summarize the feature set
python pipeline/training/generate_fingerprint_summary.py

# 4. Run training / ablation analysis
#    Save outputs to pipeline/training/outputs/

# 5. Generate manuscript figures
python pipeline/training/generate_figure_both_cv_ligand_update4_patched.py
```

---

## Expected inputs

This repository assumes access to preprocessed MD-derived data for each protein–ligand pocket system, which may include:

- Trajectory-derived time-series measurements
- Per-pocket or per-replica summary statistics
- Long-horizon endpoint metrics
- Metadata linking protein, ligand, pocket, replica, and simulation horizon

MD simulation files are typically large and stored outside the repository, referenced through local paths.

---

## Output naming convention

Output filenames encode three dimensions:

| Component | Meaning | Examples |
|-----------|---------|---------|
| Timescale | Early MD window used | `t2ns`, `t3ns`, `t5ns` |
| FPR budget | False positive rate constraint | `fpr10`, `fpr20` |
| Date stamp | Experiment run date | `20260213` |

Both `.csv` (analysis) and `.tex` (direct manuscript inclusion) formats are provided for ablation tables.

---

## Environment

**Python 3.10+** with the following dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install numpy pandas scipy scikit-learn matplotlib
```

Pin package versions when reproducing manuscript results exactly.

---

## Notes

- Some scripts may contain local path assumptions from the original experiment environment. Check input data directories, output save paths, and figure export locations before running.
- The `outputs/` directory contains cached result files from manuscript experiments to support figure and table regeneration without rerunning upstream steps.

---
