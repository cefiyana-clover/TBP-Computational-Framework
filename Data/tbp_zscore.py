# ==============================================================================
# STATISTICAL EVALUATOR - OMNICIDE V6 (DEFICIT LIKELIHOOD)
# PROTOCOL: CAUSALITY EXPANSION & DIMENSIONAL LEAKAGE (PHASE 4)
# METHOD: 1 GeV Precision Binning & Inverse Asymptotic Profile Likelihood
# TARGET: KINEMATIC ENDPOINT / HEAVISIDE CUTOFF AT 215.11111111111 GeV
# ==============================================================================

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("\n==================================================")
print("[SYSTEM] OMNICIDE: Z-SCORE STATISTICAL EVALUATOR (DEFICIT LIKELIHOOD)")
print("==================================================\n")

# ------------------------------------------------------------------------------
# 1. AKUISISI MATRIKS DATA MASSA TRANSVERSAL
# ------------------------------------------------------------------------------
DATA_FILE = "mt_all_distribution.npy"
try:
    mt_all = np.load(DATA_FILE)
    print(f"[INFO] Matriks '{DATA_FILE}' berhasil dimuat ke dalam memori.")
except FileNotFoundError:
    print(f"[FATAL ERROR] File '{DATA_FILE}' tidak ditemukan. Pastikan OMNICIDE V6 telah tereksekusi hingga tuntas.")
    exit()

# ------------------------------------------------------------------------------
# 2. DEFINISI KOORDINAT KINEMATIK (TEBING KEBOCORAN DIMENSI)
# ------------------------------------------------------------------------------
sr_min = 215.0
sr_max = 216.0
obs_events = np.sum((mt_all >= sr_min) & (mt_all <= sr_max))

fit_min = 150.0
fit_max = 250.0

# ------------------------------------------------------------------------------
# 3. PREPARASI HISTOGRAM (RESOLUSI PRESISI 1.0 GeV)
# ------------------------------------------------------------------------------
# 101 titik batas (edges) menghasilkan tepat 100 interval bin selebar 1 GeV
bins = np.linspace(fit_min, fit_max, int(fit_max - fit_min) + 1)
bin_centers = (bins[:-1] + bins[1:]) / 2
counts, _ = np.histogram(mt_all, bins=bins)

# Masking: Mengisolasi area sinyal (215-216 GeV) dari kalkulasi background
mask_cr = (bin_centers < sr_min) | (bin_centers > sr_max)
x_cr = bin_centers[mask_cr]
y_cr = counts[mask_cr]

# ------------------------------------------------------------------------------
# 4. PEMODELAN BACKGROUND TRANSVERSAL (EKSPONENSIAL DECAY)
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
# 5. EKSEKUSI STATISTIK KINEMATIC ENDPOINT (DEFICIT LIKELIHOOD)
# ------------------------------------------------------------------------------
x_sr = bin_centers[~mask_cr]
exp_bkg = np.sum(bkg_model(x_sr, *popt))

O = obs_events
B = exp_bkg

# Mengukur Signifikansi Defisit (Dimensional Leakage)
if B > 0 and O < B:
    # Menggunakan Poisson Likelihood untuk defisit (O < B)
    z_score = np.sqrt(2 * (B - O + O * np.log(O / B)))
    p_value = norm.sf(z_score)
    leakage_confirmed = True
elif B > 0 and O > B:
    # Surplus konvensional (Bukan prediksi TBP)
    z_score = np.sqrt(2 * (O * np.log(O / B) - (O - B)))
    p_value = norm.sf(z_score)
    leakage_confirmed = False
else:
    z_score = 0.0
    p_value = 1.0
    leakage_confirmed = False

# ------------------------------------------------------------------------------
# 6. OUTPUT MATEMATIS TERMINAL
# ------------------------------------------------------------------------------
print(f"Total Populasi M_T (>100 GeV) : {len(mt_all):,} events")
print(f"Observed Events (O)           : {O} events (Area 215-216 GeV)")
print(f"Expected Background (B)       : {B:.2f} events (Estimasi Fit Matematis)")
print("-" * 50)
print(f"Deviation (O - B)             : {O - B:.2f} events")
print(f"Z-Score (Significance)        : {z_score:.5f} Sigma")
print(f"P-Value                       : {p_value:.5e}")

if leakage_confirmed and z_score >= 5.0:
    print("[STATUS] KLAIM DISCOVERY VALID (≥ 5.0 SIGMA). HEAVISIDE CUTOFF / LEAKAGE CONFIRMED.")
elif leakage_confirmed and z_score >= 3.0:
    print("[STATUS] EVIDENCE KUAT (3.0 - 4.9 SIGMA). DIMENSIONAL LEAKAGE IN PROGRESS.")
elif not leakage_confirmed and z_score >= 3.0:
    print("[STATUS] PERINGATAN: SURPLUS DITEMUKAN. KONSISTEN DENGAN RESONANSI PARTIKEL, BUKAN LEAKAGE.")
else:
    print("[STATUS] FLUKTUASI KONSISTEN DENGAN BACKGROUND MODEL M_T.")
print("==================================================\n")

# ------------------------------------------------------------------------------
# 7. RENDER GRAFIK DIAGNOSTIK KINEMATIK
# ------------------------------------------------------------------------------
if fit_success:
    plt.figure(figsize=(10, 6))
    plt.hist(mt_all[(mt_all >= fit_min) & (mt_all <= fit_max)], bins=bins, color='#1e3d59', alpha=0.8, label='Observed M_T Distribution')
    
    x_plot = np.linspace(fit_min, fit_max, 500)
    plt.plot(x_plot, bkg_model(x_plot, *popt), color='red', linestyle='dashed', linewidth=2.5, label='Background Fit Model')
    
    plt.axvspan(sr_min, sr_max, color='yellow', alpha=0.4, label='Heaviside Cutoff Boundary (215.111 GeV)')
    
    plt.yscale('log')
    plt.xlim(fit_min, fit_max)
    plt.xlabel(r'Full Transverse Mass $M_T$ [GeV]', fontsize=12)
    plt.ylabel('Events / 1.0 GeV', fontsize=12)
    plt.title('Kinematic Endpoint Analysis: M_T Deficit vs Heaviside Cutoff', fontweight='bold')
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.2)
    
    OUTPUT_IMG = "Z_Score_Verification_Deficit_FullScale.png"
    plt.savefig(OUTPUT_IMG, dpi=300, bbox_inches='tight')
    print(f"[SUCCESS] Grafik verifikasi fitting M_T disimpan sebagai '{OUTPUT_IMG}'.", flush=True)
