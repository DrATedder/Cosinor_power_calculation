#!/usr/bin/env python3
"""
multi_gene_cosinor_power.py

Batch cosinor power and sample-size estimation for multiple genes (TSV files).
- Reads all files matching "*.tsv" in the working directory.
- Each file should contain tab-separated columns (no header required but must be in this order):
    mean, sd, cv, mean_hk, sd_hk, cv_hk, ZT, expr_ratio
  or the header names may match those labels.
- Fits a single-component cosinor to the per-ZT means to estimate mesor, amplitude, phase.
- Estimates technical and biological noise from housekeeping statistics.
- Produces power-vs-n sample-size curves for each gene (saved as a single PNG).
- Outputs per-gene summary CSV with fitted parameters and estimated required n for 80% power.
- Optional transcriptome-level simulation using BH-FDR to estimate discovery power across many genes.

Requirements:
    numpy, pandas, scipy, statsmodels, matplotlib

Usage: edit CONFIG or run as-is in the directory containing Genename.tsv files.
"""

import os, glob, io, math, time, argparse
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

# ---------------------- CONFIG ----------------------
PERIOD = 24.0
TIMES_DEFAULT = np.array([0,4,8,12,16,20])   # assumed timepoints; script will use per-file ZTs if present
TARGET_POWER = 0.80
ALPHA = 0.05
N_SIM_POWER = 800         # Monte-Carlo iterations for per-gene power curve
MAX_N = 100               # max n to search for per-gene sample-size (per-group)
N_SIM_TRANS = 200         # Monte-Carlo iterations per n for transcriptome simulation (keep smaller for speed)
G_TOTAL = 20000           # total genes for transcriptome simulation
PI1 = 0.40                # fraction of truly rhythmic genes in transcriptome sim
ALPHA_FDR = 0.05          # target FDR for transcriptome
OUT_PREFIX = "cosinor_batch"   # output file prefix
# ---------------------------------------------------

def read_gene_file(path):
    # Try to read with header; fall back to fixed columns order if no header
    try:
        df = pd.read_csv(path, sep=None, engine='python')  # auto-detect sep
    except Exception:
        df = pd.read_csv(path, sep='\t', header=None)
    # normalize column names (lowercase, remove spaces)
    df.columns = [str(c).strip().lower().replace(' ', '_') for c in df.columns]
    # try to find expression column
    expr_col = None
    for c in df.columns:
        if c in ('expr_ratio','expression_ratio','expression_ratio_(goi/hk)','expr_ratio(goi/hk)'):
            expr_col = c
            break
        if 'expr' in c and ('ratio' in c or 'hk' in c):
            expr_col = c
            break
        if 'expr' in c and 'ratio' in c:
            expr_col = c
            break
    if expr_col is None:
        # fallback: assume last column is expression ratio
        expr_col = df.columns[-1]
    # ensure ZT column exists
    zt_col = None
    for c in df.columns:
        if c in ('zt','time') or 'zt' in c:
            zt_col = c
            break
    if zt_col is None:
        raise ValueError(f"No ZT column found in {path}; expected column named 'ZT' or similar.")
    df2 = pd.DataFrame({'ZT': df[zt_col].astype(float).astype(int), 'expr_ratio': df[expr_col].astype(float)})
    # attach mean_hk and sd_hk if present
    if 'mean_hk' in df.columns:
        df2['mean_hk'] = df['mean_hk']
    else:
        df2['mean_hk'] = np.nan
    if 'sd_hk' in df.columns:
        df2['sd_hk'] = df['sd_hk']
    else:
        df2['sd_hk'] = np.nan
    return df2

def fit_cosinor_perZT(df):
    summary = df.groupby('ZT').agg(mean_expr=('expr_ratio','mean'), sd_expr=('expr_ratio','std')).reset_index().sort_values('ZT')
    times = summary['ZT'].values
    y_mean = summary['mean_expr'].values
    omega = 2*math.pi/PERIOD
    X = np.column_stack([np.cos(omega*times), np.sin(omega*times)])
    X = sm.add_constant(X)
    model = sm.OLS(y_mean, X).fit()
    C, b_cos, b_sin = model.params[0], model.params[1], model.params[2]
    amp = math.sqrt(b_cos**2 + b_sin**2)
    phi_rad = math.atan2(-b_sin, b_cos)
    phi_hours = (phi_rad * (24.0/(2*math.pi))) % 24
    return {'mesor': C, 'amp': amp, 'phase_hours': phi_hours, 'perZT_summary': summary}

def estimate_sigmas(df):
    summary = df.groupby('ZT').agg(sd_expr=('expr_ratio','std')).reset_index()
    within_z_sd2 = summary['sd_expr'].fillna(0).apply(lambda x: x**2).values
    sigma_total = math.sqrt(np.mean(within_z_sd2))
    # technical sigma approx from average hk CV * mean expr ratio
    if df['mean_hk'].notna().sum() > 0:
        df_temp = df.copy()
        df_temp['tech_cv'] = df_temp['sd_hk'] / (df_temp['mean_hk'].replace(0, np.nan))
        tech_cv_mean = df_temp['tech_cv'].dropna().mean()
        sigma_tech = tech_cv_mean * df['expr_ratio'].mean()
    else:
        sigma_tech = 0.0
    sigma_bio = math.sqrt(max(0, sigma_total**2 - sigma_tech**2))
    return {'sigma_total': sigma_total, 'sigma_tech': sigma_tech, 'sigma_bio': sigma_bio}

def simulate_once_pooled(mesor, amp, phase_hours, sigma, subj_sd, n_subj, times):
    omega = 2*math.pi/PERIOD
    t = np.tile(times, n_subj)
    subj = np.repeat(np.arange(n_subj), len(times))
    signal = mesor + amp * np.cos(omega*(t - phase_hours))
    bis = np.random.normal(0, subj_sd, n_subj)
    y = signal + bis[subj] + np.random.normal(0, sigma, len(t))
    cos_t = np.cos(omega*(t - phase_hours))
    sin_t = np.sin(omega*(t - phase_hours))
    X = sm.add_constant(np.column_stack([cos_t, sin_t]))
    model = sm.OLS(y, X).fit()
    R = np.zeros((2, len(model.params)))
    R[0,1] = 1
    R[1,2] = 1
    ftest = model.f_test(R)
    return float(ftest.pvalue)

def estimate_power(mesor, amp, phase_hours, sigma, subj_sd, n_subj, times, n_sim=N_SIM_POWER, alpha=ALPHA):
    pvals = []
    for i in range(n_sim):
        p = simulate_once_pooled(mesor, amp, phase_hours, sigma, subj_sd, n_subj, times)
        pvals.append(p)
    return np.mean(np.array(pvals) < alpha)

def power_curve_for_gene(mesor, amp, phase_hours, sigma_bio, subj_sd_assumption, times, n_min=2, n_max=MAX_N, n_sim=N_SIM_POWER, alpha=ALPHA):
    ns = list(range(n_min, n_max+1))
    powers = []
    for n in ns:
        pwr = estimate_power(mesor, amp, phase_hours, sigma_bio, subj_sd_assumption, n, times, n_sim=n_sim, alpha=alpha)
        powers.append(pwr)
    return np.array(ns), np.array(powers)

def find_n_for_target(ns, powers, target_power=TARGET_POWER):
    idx = np.where(powers >= target_power)[0]
    if len(idx)==0:
        return None
    return int(ns[idx[0]])

def transcriptome_simulation(empirical_params, G=G_TOTAL, pi1=PI1, times=TIMES_DEFAULT, alpha_fdr=ALPHA_FDR, n_sim=N_SIM_TRANS):
    amps = np.array([p['amp'] for p in empirical_params])
    sigmas = np.array([p['sigma_bio'] for p in empirical_params])
    def sample_params(G):
        amps_s = np.random.choice(amps, size=G, replace=True)
        sigs_s = np.random.choice(sigmas, size=G, replace=True)
        return amps_s, sigs_s
    ns = list(range(2, min(201, MAX_N+1), 2))
    avg_power_per_n = []
    for n in ns:
        power_runs = []
        for sim in range(n_sim):
            amps_s, sigs_s = sample_params(G)
            is_rhyth = np.random.rand(G) < pi1
            amps_true = amps_s * is_rhyth
            pvals = np.empty(G)
            for i in range(G):
                p = simulate_once_pooled(mesor=1.0, amp=amps_true[i], phase_hours=0.0, sigma=sigs_s[i], subj_sd=0.0, n_subj=n, times=times)
                pvals[i] = p
            rejected, pvals_adj, _, _ = multipletests(pvals, alpha=alpha_fdr, method='fdr_bh')
            if is_rhyth.sum() == 0:
                power_runs.append(0.0)
            else:
                tp = np.logical_and(rejected, is_rhyth).sum()
                power_runs.append(tp / is_rhyth.sum())
        avg_power_per_n.append(np.mean(power_runs))
        print(f"Transcriptome sim: n={n} -> discovery power (mean) = {avg_power_per_n[-1]:.3f}")
    return ns, avg_power_per_n

def main_process(folder='.', pattern='*.tsv', output_prefix=OUT_PREFIX, run_transcriptome=True):
    files = sorted(glob.glob(os.path.join(folder, pattern)))
    if len(files) == 0:
        print("No TSV files found in folder. Place Genename.tsv files in the folder and re-run.")
        return
    per_gene_results = []
    empirical_params = []
    plt.figure(figsize=(10,6))
    for path in files:
        gene = os.path.splitext(os.path.basename(path))[0]
        print(f"Processing {gene} ...")
        try:
            df = read_gene_file(path)
        except Exception as e:
            print(f"  Failed to read {path}: {e}")
            continue
        fit = fit_cosinor_perZT(df)
        sigs = estimate_sigmas(df)
        times = np.array(sorted(df['ZT'].unique()))
        if len(times)==0:
            times = TIMES_DEFAULT
        subj_sd_assump = sigs['sigma_bio'] / 2.0
        ns, powers = power_curve_for_gene(fit['mesor'], fit['amp'], fit['phase_hours'], sigs['sigma_bio'], subj_sd_assump, times, n_min=2, n_max=MAX_N, n_sim=N_SIM_POWER, alpha=ALPHA)
        req_n = find_n_for_target(ns, powers, TARGET_POWER)
        per_gene_results.append({
            'gene': gene,
            'mesor': fit['mesor'],
            'amp': fit['amp'],
            'phase_h': fit['phase_hours'],
            'sigma_total': sigs['sigma_total'],
            'sigma_tech': sigs['sigma_tech'],
            'sigma_bio': sigs['sigma_bio'],
            'req_n_80': req_n
        })
        empirical_params.append({'amp': fit['amp'], 'sigma_bio': sigs['sigma_bio']})
        plt.plot(ns, powers, label=gene, marker='o', linewidth=1)
    plt.axhline(TARGET_POWER, color='k', linestyle='--', label=f"target {TARGET_POWER*100:.0f}%")
    plt.xlabel('n (subjects per group)')
    plt.ylabel('Power to detect non-zero amplitude')
    plt.title('Power curves per gene (pooled cosinor)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plot_path = f"{output_prefix}_power_curves.png"
    plt.savefig(plot_path, dpi=200)
    print(f"Saved power curves figure to {plot_path}")
    out_df = pd.DataFrame(per_gene_results).sort_values('req_n_80')
    csv_path = f"{output_prefix}_per_gene_summary.csv"
    out_df.to_csv(csv_path, index=False)
    print(f"Saved per-gene summary to {csv_path}")
    if run_transcriptome and len(empirical_params)>=5:
        print("\nRunning transcriptome-level simulations (this may take some time).\nAdjust G_TOTAL, PI1, N_SIM_TRANS in the script for speed.)")
        ns_tx, power_tx = transcriptome_simulation(empirical_params, G=G_TOTAL, pi1=PI1, times=TIMES_DEFAULT, alpha_fdr=ALPHA_FDR, n_sim=N_SIM_TRANS)
        tx_df = pd.DataFrame({'n': ns_tx, 'disc_power': power_tx})
        tx_path = f"{output_prefix}_transcriptome_power.csv"
        tx_df.to_csv(tx_path, index=False)
        print(f"Saved transcriptome simulation results to {tx_path}")
        plt.figure(figsize=(6,4))
        plt.plot(ns_tx, power_tx, marker='o')
        plt.xlabel('n (subjects per group)')
        plt.ylabel('Mean discovery power (BH-FDR)')
        plt.title(f'Transcriptome discovery power (G={G_TOTAL}, pi1={PI1})')
        plt.grid(True)
        tx_plot = f"{output_prefix}_transcriptome_power.png"
        plt.savefig(tx_plot, dpi=200)
        print(f"Saved transcriptome plot to {tx_plot}")
    else:
        print("Not enough genes to run transcriptome sim (need >=5 empirical genes) or transcriptome sim disabled.")

if __name__ == '__main__':
    t0 = time.time()
    main_process(folder='.', pattern='*.tsv', output_prefix=OUT_PREFIX, run_transcriptome=True)
    print('Done. Elapsed time: {:.1f}s'.format(time.time()-t0))

