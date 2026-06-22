# ==============================================================================
# STATISTICAL SIGNIFICANCE CALCULATOR - SIDEBAND ESTIMATION
# METHOD: Model-Independent Symmetric Sideband Estimation
# TARGET: 215.11111111111 GeV
# ==============================================================================

import os
import numpy as np

# ------------------------------------------------------------------------------
# 1. KONTROL PARAMETER STATISTIK
# ------------------------------------------------------------------------------
TBP_MASS = 215.11111111111
SIGNAL_WINDOW = 0.5  

SR_LOW, SR_HIGH = TBP_MASS - SIGNAL_WINDOW, TBP_MASS + SIGNAL_WINDOW
LSB_LOW, LSB_HIGH = 210.0, 214.0  
RSB_LOW, RSB_HIGH = 216.5, 220.5  

print("==================================================")
print("[SYSTEM] STATISTICAL SIGNIFICANCE EVALUATOR")
print("==================================================\n")

# ------------------------------------------------------------------------------
# 2. INGESTI DATA
# ------------------------------------------------------------------------------
DATA_FILE = "met_all_distribution.npy"

if os.path.exists(DATA_FILE):
    print(f"[INFO] File '{DATA_FILE}' terdeteksi.")
    try:
        met_data = np.load(DATA_FILE)
        print(f"[STATUS] {len(met_data)} sampel MET dimuat ke memori.")
        
        n_obs = np.sum((met_data >= SR_LOW) & (met_data <= SR_HIGH))
        n_lsb = np.sum((met_data >= LSB_LOW) & (met_data <= LSB_HIGH))
        n_rsb = np.sum((met_data >= RSB_LOW) & (met_data <= RSB_HIGH))
        
    except Exception as e:
        print(f"[ERROR] Gagal membaca file npy: {e}. Inisialisasi input manual.")
        n_obs = int(input(f"Input n_obs ({SR_LOW:.2f}-{SR_HIGH:.2f} GeV): "))
        n_lsb = int(input(f"Input n_lsb ({LSB_LOW:.2f}-{LSB_HIGH:.2f} GeV): "))
        n_rsb = int(input(f"Input n_rsb ({RSB_LOW:.2f}-{RSB_HIGH:.2f} GeV): "))
else:
    print(f"[WARNING] File '{DATA_FILE}' tidak ditemukan. Inisialisasi input manual.")
    n_obs = int(input(f"Input n_obs ({SR_LOW:.2f}-{SR_HIGH:.2f} GeV): "))
    n_lsb = int(input(f"Input n_lsb ({LSB_LOW:.2f}-{LSB_HIGH:.2f} GeV): "))
    n_rsb = int(input(f"Input n_rsb ({RSB_LOW:.2f}-{RSB_HIGH:.2f} GeV): "))

# ------------------------------------------------------------------------------
# 3. MATRIKS BACKGROUND (MODEL-INDEPENDENT)
# ------------------------------------------------------------------------------
w_sr = SR_HIGH - SR_LOW    
w_sb = (LSB_HIGH - LSB_LOW) + (RSB_HIGH - RSB_LOW)  

alpha = w_sr / w_sb  

n_sideband = n_lsb + n_rsb
n_bkg = alpha * n_sideband

print("\n--------------------------------------------------")
print("[DATA] EVENT CONTINGENCY MATRIX")
print("--------------------------------------------------")
print(f" Signal Region (SR)   [{SR_LOW:.2f} - {SR_HIGH:.2f} GeV] : {n_obs} events (N_obs)")
print(f" Left Sideband (LSB)  [{LSB_LOW:.2f} - {LSB_HIGH:.2f} GeV] : {n_lsb} events")
print(f" Right Sideband (RSB) [{RSB_LOW:.2f} - {RSB_HIGH:.2f} GeV] : {n_rsb} events")
print(f" Total Sideband Control                      : {n_sideband} events")
print(f" Estimated Background (N_bkg)                : {n_bkg:.4f} events")

# ------------------------------------------------------------------------------
# 4. KALKULASI Z-SCORE
# ------------------------------------------------------------------------------
if n_bkg <= 0:
    print("\n[ERROR] Background absolut tidak mencukupi untuk kalkulasi statistik.")
else:
    z_standard = (n_obs - n_bkg) / np.sqrt(n_bkg)
    
    if n_obs > 0:
        z_lhc = np.sqrt(2 * (n_obs * np.log(n_obs / n_bkg) - (n_obs - n_bkg)))
        if n_obs < n_bkg:
            z_lhc = -z_lhc
    else:
        z_lhc = 0.0

    print("\n--------------------------------------------------")
    print("[METRICS] SIGNIFICANCE EVALUATION (Z-SCORE)")
    print("--------------------------------------------------")
    print(f" Standard Significance (S/√B)   : {z_standard:.4f} Sigma")
    print(f" Profile Likelihood Ratio (LHC) : {z_lhc:.4f} Sigma")
    print("--------------------------------------------------")
    
    print("\n[STATUS] STATISTICAL CONCLUSION:")
    if z_lhc >= 5.0:
        print(" [RESULT] SIGNIFICANCE >= 5.0 SIGMA (DISCOVERY THRESHOLD MET).")
        print(" Data memenuhi kriteria observasi fisis anomali pada 215.11 GeV.")
    elif z_lhc >= 3.0:
        print(" [RESULT] SIGNIFICANCE >= 3.0 SIGMA (EVIDENCE THRESHOLD MET).")
        print(" Dibutuhkan perluasan dataset untuk mencapai ambang batas 5-Sigma.")
    else:
        print(" [RESULT] SIGNIFICANCE < 3.0 SIGMA.")
        print(" Fluktuasi konsisten dengan model background standar.")

print("\n==================================================")
