# Cosinor_power_calculation

## Overview

`multi_gene_cosinor_power.py` performs batch cosinor-based power and sample-size estimation for multiple genes using rhythmic expression data. It fits single-component cosinor models, estimates technical and biological variance, and simulates detection power across a range of sample sizes. Optionally, it can perform transcriptome-level discovery power simulations with multiple testing correction (Benjamini–Hochberg FDR).

---

## Method Summary

Each gene’s rhythmic expression is modeled as:

\[
Y(t) = M + A \cos(\omega t + \phi) + \epsilon
\]

where:

- \( M \): MESOR (mean level)  
- \( A \): amplitude  
- \( \phi \): phase (in hours)  
- \( \omega = 2\pi/24 \): angular frequency (assuming a 24-hour period)  
- \( \epsilon \): residual noise  

The model is fit via ordinary least squares, and power is estimated through Monte Carlo simulation. Each iteration simulates datasets given the fitted parameters and tests for rhythmicity (non-zero amplitude) using an F-test.

---

## Script Description

**Filename:** `multi_gene_cosinor_power.py`

This script:

1. Reads all `.tsv` files in the working directory.
2. Extracts gene names from filenames (e.g., `Bmal1.tsv` → `Bmal1`).
3. Fits a cosinor model for each gene.
4. Estimates technical and biological variance components using housekeeping gene data.
5. Runs power simulations over sample sizes \( n = 2 \) to `MAX_N`.
6. Determines the minimum \( n \) required to achieve target power (default 80%).
7. Optionally runs a transcriptome-scale simulation to estimate discovery power under multiple testing.

### Outputs

- Power curves for all genes (`.png`)
- Per-gene summary of cosinor parameters and sample-size estimates (`.csv`)
- Transcriptome-level discovery power curves and tables (optional)

---

## Input Data Format

Each gene is represented by a single `.tsv` file named:

Expected columns:

| mean | sd | cv | mean_hk | sd_hk | cv_hk | ZT | expr_ratio |
|------|----|----|----------|-------|-------|----|-------------|

- **ZT:** Zeitgeber time (e.g., 0, 4, 8, 12, 16, 20)
- **expr_ratio:** Expression ratio of gene of interest relative to housekeeping gene (e.g., GOI/HPRT1)
- **mean_hk**, **sd_hk**, **cv_hk:** Used to estimate technical variation
- **NA values:** Allowed and automatically handled during analysis

Example:
Bmal1.tsv; 
Per2.tsv; 
Nr1d1.tsv

---

## Configuration Parameters

The following parameters are defined in the 'CONFIG' section of the script and can be edited as needed:

```python
PERIOD = 24.0                 # assumed circadian period
TIMES_DEFAULT = [0, 4, 8, 12, 16, 20]
TARGET_POWER = 0.80
ALPHA = 0.05
N_SIM_POWER = 800             # Monte Carlo iterations per gene
MAX_N = 100                   # maximum n (per group) to test
N_SIM_TRANS = 200             # iterations per n for transcriptome simulation
G_TOTAL = 20000               # total genes in transcriptome simulation
PI1 = 0.40                    # proportion of rhythmic genes
ALPHA_FDR = 0.05              # FDR level for transcriptome simulation
OUT_PREFIX = "cosinor_batch"  # output file prefix
```

---

## Dependencies

Requires Python 3.9+ and the following packages:

```python
numpy
pandas
scipy
statsmodels
matplotlib
```

Install with:

```bash
pip3 install numpy pandas scipy statsmodels matplotlib
```

If using `NumPy ≥2.0 causes issues, install version 1.x:

```bash
pip install "numpy<2"
```
---

## Usage

1. Place all gene .tsv files in the same directory as the script.
2. Open a terminal in that directory.
3. Run:

```bash
python3 multi_gene_cosinor_power.py
```
The script will automatically:

* Process all .tsv files
* Generate power vs. sample-size curves for each gene
* Create a summary table of fitted parameters and estimated sample sizes
* **Optionally** simulate transcriptome-level discovery power (requires ≥5 genes)

--- 

## Outputs

| File                                    | Description                                                     |
| --------------------------------------- | --------------------------------------------------------------- |
| `cosinor_batch_power_curves.png`        | Power curves for all genes                                      |
| `cosinor_batch_per_gene_summary.csv`    | Summary table of fitted parameters and required n for 80% power |
| `cosinor_batch_transcriptome_power.csv` | Transcriptome discovery power estimates                         |
| `cosinor_batch_transcriptome_power.png` | Transcriptome discovery power plot                              |

| gene  | mesor | amp  | phase_h | sigma_total | sigma_tech | sigma_bio | req_n_80 |
| ----- | ----- | ---- | ------- | ----------- | ---------- | --------- | -------- |
| Bmal1 | 1.25  | 0.84 | 16.0    | 0.24        | 0.10       | 0.22      | 14       |
| Per2  | 1.05  | 0.30 | 8.0     | 0.28        | 0.09       | 0.27      | 42       |

---

## Runtime Estimates

| Task                                              | Approx. Runtime |
| ------------------------------------------------- | --------------- |
| Single gene (N_SIM_POWER=800)                     | 1–2 min         |
| 10 genes                                          | 10–15 min       |
| Transcriptome simulation (G=20k, N_SIM_TRANS=200) | 10–20 min       |

**Note**. Reduce `N_SIM_POWER` or `N_SIM_TRANS` for quicker exploratory runs.

---

## Interpretation of Results

* **amp**: Estimated amplitude of rhythmic oscillation
* **phase_h**: Estimated acrophase (hours of peak expression)
* **sigma_bio**: Biological variance estimate
* **req_n_80**: Estimated minimum per-group sample size to achieve 80% power

Larger amplitude and lower variance reduce the required sample size.
Transcriptome-level simulations summarize expected discovery power after multiple testing correction.

---

## Example Workflow

1. Prepare your time-course expression `.tsv` file(s) (one per gene).
2. Run `multi_gene_cosinor_power.py` in the directory which contains your `.tsv` file(s).
3. Review `cosinor_batch_power_curves.png` for visual power trends.
4. Consult `cosinor_batch_per_gene_summary.csv` for required sample sizes.
5. Review transcriptome simulation outputs if applicable.

---

## Notes

* Assumes a single-component 24-hour cosinor model.
* Handles missing data (NA).
* Power estimation is fully empirical via Monte Carlo simulation.
* Transcriptome simulation requires ≥5 genes to build empirical amplitude/variance distributions.

---

## References

* Nelson W., Tong Y.L., Lee J.K., Halberg F. (1979). Methods for cosinor-rhythmometry. Chronobiologia, 6(4), 305–323.
* Refinetti R., Cornélissen G., Halberg F. (2007). Procedures for numerical analysis of circadian rhythms. Biological Rhythm Research, 38(4), 275–325.

---


