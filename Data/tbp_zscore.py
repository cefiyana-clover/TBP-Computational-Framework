# ==============================================================================
# STATISTICAL EVALUATOR - OMNICIDE V4 (FULL SCALE)
# PROTOCOL: CAUSALITY EXPANSION & POSITIVE IMPLANTATION (PHASE 4)
# METHOD: 1 GeV Precision Binning & Asymptotic Profile Likelihood Ratio
# TARGET: ULTRA-STRICT KINEMATIC MATRIX
# ==============================================================================

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("\n==================================================")
print("[SYSTEM] OMNICIDE: Z-SCORE STATISTICAL EVALUATOR V4 (ULTRA-STRICT)")
print("==================================================\n")

# ------------------------------------------------------------------------------
# 1. AKUISISI MATRIKS DATA
# ------------------------------------------------------------------------------
DATA_FILE = "met_all_distribution.npy"
try:
    met_all = np.load(DATA_FILE)
    print(f"[INFO] Matriks '{DATA_FILE}' berhasil dimuat ke dalam memori.")
except FileNotFoundError:
    print(f"[FATAL ERROR] File '{DATA_FILE}' tidak ditemukan. Pastikan OMNICIDE V4 telah tereksekusi hingga tuntas.")
    exit()

# ------------------------------------------------------------------------------
# 2. DEFINISI KOORDINAT KINEMATIK
# ------------------------------------------------------------------------------
sr_min = 215.0
sr_max = 216.0
obs_events = np.sum((met_all >= sr_min) & (met_all <= sr_max))

fit_min = 150.0
fit_max = 250.0

# ------------------------------------------------------------------------------
# 3. PREPARASI HISTOGRAM (RESOLUSI PRESISI 1.0 GeV)
# ------------------------------------------------------------------------------
# 101 titik batas (edges) menghasilkan tepat 100 interval bin selebar 1 GeV
bins = np.linspace(fit_min, fit_max, int(fit_max - fit_min) + 1)
bin_centers = (bins[:-1] + bins[1:]) / 2
counts, _ = np.histogram(met_all, bins=bins)

# Masking: Mengisolasi area sinyal (215-216 GeV) dari kalkulasi background
mask_cr = (bin_centers < sr_min) | (bin_centers > sr_max)
x_cr = bin_centers[mask_cr]
y_cr = counts[mask_cr]

# ------------------------------------------------------------------------------
# 4. PEMODELAN BACKGROUND (EKSPONENSIAL DECAY)
# ------------------------------------------------------------------------------
def bkg_model(x, a, b):
    return a * np.exp(-b * x)

try:
    popt, pcov = curve_fit(bkg_model, x_cr, y_cr, p0=(1e6, 0.05), maxfev=10000)
    fit_success = True
except Exception as e:
    print(f"[ERROR] Kegagalan pemodelan fisis latar belakang: {e}")
    popt = [0, 0]
    fit_success = False

# ------------------------------------------------------------------------------
# 5. EKSEKUSI STATISTIK (PROFILE LIKELIHOOD RATIO)
# ------------------------------------------------------------------------------
x_sr = bin_centers[~mask_cr]
exp_bkg = np.sum(bkg_model(x_sr, *popt))

O = obs_events
B = exp_bkg

# Kalkulasi Asimtotik LHC (Cowan et al.)
if B > 0 and O > B:
    z_score = np.sqrt(2 * (O * np.log(O / B) - (O - B)))
    p_value = norm.sf(z_score)
else:
    z_score = 0.0
    p_value = 1.0

# ------------------------------------------------------------------------------
# 6. OUTPUT MATEMATIS TERMINAL
# ------------------------------------------------------------------------------
print(f"Total Populasi MET (>100 GeV) : {len(met_all):,} events")
print(f"Observed Events (O)           : {O} anomali (Area 215-216 GeV)")
print(f"Expected Background (B)       : {B:.2f} events (Estimasi Fit Matematis)")
print("-" * 50)
print(f"Z-Score (Significance)        : {z_score:.5f} Sigma")
print(f"P-Value                       : {p_value:.5e}")

if z_score >= 5.0:
    print("[STATUS] KLAIM DISCOVERY VALID (≥ 5.0 SIGMA). CAUSALITY EXPANSION CONFIRMED.")
elif z_score >= 3.0:
    print("[STATUS] EVIDENCE KUAT (3.0 - 4.9 SIGMA). POSITIVE IMPLANTATION IN PROGRESS.")
else:
    print("[STATUS] FLUKTUASI KONSISTEN DENGAN BACKGROUND MODEL.")
print("==================================================\n")

# ------------------------------------------------------------------------------
# 7. RENDER GRAFIK DIAGNOSTIK
# ------------------------------------------------------------------------------
if fit_success:
    plt.figure(figsize=(10, 6))
    plt.hist(met_all[(met_all >= fit_min) & (met_all <= fit_max)], bins=bins, color='#333333', alpha=0.8, label='Observed Background (Ultra-Strict Cut)')
    
    x_plot = np.linspace(fit_min, fit_max, 500)
    plt.plot(x_plot, bkg_model(x_plot, *popt), color='red', linestyle='dashed', linewidth=2.5, label='Background Fit Model')
    
    plt.axvspan(sr_min, sr_max, color='yellow', alpha=0.4, label='Signal Region (Anomalies)')
    
    plt.yscale('log')
    plt.xlim(fit_min, fit_max)
    plt.xlabel(r'Missing Transverse Energy $E_T^{miss}$ [GeV]', fontsize=12)
    plt.ylabel('Events / 1.0 GeV', fontsize=12)
    plt.title('Ultra-Strict Full Scale Fit Model vs Target Anomalies', fontweight='bold')
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.2)
    
    OUTPUT_IMG = "Z_Score_Verification_UltraStrict_FullScale.png"
    plt.savefig(OUTPUT_IMG, dpi=300, bbox_inches='tight')
    print(f"[SUCCESS] Grafik verifikasi fitting disimpan sebagai '{OUTPUT_IMG}'.", flush=True)
