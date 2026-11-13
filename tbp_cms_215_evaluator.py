# ==============================================================================
# TBP-OMNI-EVALUATOR: UNIFIED STATISTICAL INFERENCE ENGINE (215.111 GeV CEILING)
# PROTOCOL: KINEMATIC PLR, MET DIVERGENCE, AND LUMINOSITY SCALING
# REPOSITORY: https://github.com/cefiyana-clover/TBP-Computational-Framework
# ==============================================================================
import os
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm, chi2, ks_2samp, ttest_ind
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

print("==================================================")
print(" [INIT] OMNI-EVALUATOR STATISTICAL INFERENCE ENGINE ")
print("==================================================")

# ------------------------------------------------------------------------------
# 0. AUTOMATED DATA RESOLUTION (KAGGLE & LOCAL REPOSITORY COMPATIBLE)
# ------------------------------------------------------------------------------
DEFAULT_DIR = "/kaggle/" if os.path.exists("/kaggle/") else "./"
BASE_DIR = os.environ.get("TBP_DATA_DIR", DEFAULT_DIR)
print(f"[INFO] Scanning for data matrices within {BASE_DIR}...")

npy_files = {}
for root, dirs, files in os.walk(BASE_DIR):
    for f in files:
        if f.endswith('.npy'):
            npy_files[f] = os.path.join(root, f)

print(f"[INFO] Located {len(npy_files)} .npy matrices. Commencing evaluation...\n")

# ------------------------------------------------------------------------------
# OPERATIONAL CONFIGURATION
# ------------------------------------------------------------------------------
WINDOW_LOW = 214.0
WINDOW_HIGH = 216.0

# 1. Z-SCORE EVALUATION CHANNELS (10 Sectors)
independent_channels = {
    "Cross-Flavor Di-Boson (MuonEG)": "muoneg_mass.npy",
    "Monophoton (Single-Photon)": "single_photon_mass.npy",
    "Single-Electron": "single_ele_mass.npy",
    "Single-Muon": "single_mu_mass.npy",
    "Tau": "single_tau_mass.npy",
    "Hadronic (JetHT)": "jetht_mass.npy",
    "Heavy Flavor (B-Jets)": "bjet_mass.npy",
    "W-Boson (Leptonic MT)": "wboson_mt.npy"
}

joint_channels = {
    "Di-Photon": "di_photon_mass.npy",
    "4-Lepton": "four_lepton_mass.npy"
}

# 2. KINEMATIC DIVERGENCE CHANNELS (Mass, MET_PT)
met_channels = {
    "Di-Photon": ("di_photon_mass.npy", "di_photon_met_pt.npy"),
    "4-Lepton": ("four_lepton_mass.npy", "four_lepton_met_pt.npy"),
    "MuonEG": ("muoneg_mass.npy", "muoneg_met_pt.npy"),
    "Single-Photon": ("single_photon_mass.npy", "single_photon_met_pt.npy"),
    "Single-Electron": ("single_ele_mass.npy", "single_ele_met_pt.npy"),
    "Single-Muon": ("single_mu_mass.npy", "single_mu_met_pt.npy"),
    "Tau": ("single_tau_mass.npy", "single_tau_met_pt.npy"),
    "JetHT": ("jetht_mass.npy", "jetht_met_pt.npy"),
    "Heavy Flavor (B-Jets)": ("bjet_mass.npy", "bjet_met_pt.npy"),
    "W-Boson (Leptonic MT)": ("wboson_mt.npy", "wboson_met_pt.npy")
}

# ==============================================================================
# MODULE 1: PROFILE LIKELIHOOD RATIO (PLR) EVALUATION
# ==============================================================================
print("\n==================================================")
print(" [MODULE 1] PROFILE LIKELIHOOD RATIO (PLR) EVALUATION ")
print("==================================================\n")

def background_model(m, c0, c1, c2):
    x = m / 13000.0
    x = np.clip(x, 1e-5, 1.0)
    return c0 * ((1.0 - x)**c1) / (x**c2)

def poisson_nll(expected, observed):
    if expected <= 0:
        return np.inf
    return expected - observed * np.log(expected)

def evaluate_plr_independent(mass_array, bins=100):
    mask_range = (mass_array >= 150.0) & (mass_array <= 250.0)
    m_fit = mass_array[mask_range]

    mask_sidebands = (m_fit < WINDOW_LOW) | (m_fit > WINDOW_HIGH)
    m_sidebands = m_fit[mask_sidebands]

    counts, bin_edges = np.histogram(m_sidebands, bins=bins, range=(150, 250))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

    valid_bins = (counts > 0) & ((bin_centers < WINDOW_LOW) | (bin_centers > WINDOW_HIGH))
    x_data = bin_centers[valid_bins]
    y_data = counts[valid_bins]

    if len(x_data) < 5:
        return np.sum((m_fit >= WINDOW_LOW) & (m_fit <= WINDOW_HIGH)), 0.0, 0.0

    def loss_func(params):
        c0, c1, c2 = params
        if c0 <= 0:
            return 1e12
        y_pred = background_model(x_data, c0, c1, c2)
        y_pred = np.clip(y_pred, 1e-10, 1e10)
        return np.sum(y_pred - y_data * np.log(y_pred))

    avg_count = np.mean(y_data)
    c0_guess = avg_count * ((200.0 / 13000.0)**5) if avg_count > 0 else 1e-5

    res = minimize(loss_func, x0=[c0_guess, 5.0, 5.0], method='Nelder-Mead', options={'maxiter': 10000})
    c0_fit, c1_fit, c2_fit = res.x

    obs_signal = np.sum((m_fit >= WINDOW_LOW) & (m_fit <= WINDOW_HIGH))
    bin_width = (250.0 - 150.0) / float(bins)
    m_signal_axis = np.linspace(WINDOW_LOW, WINDOW_HIGH, bins)

    bkg_density = background_model(m_signal_axis, c0_fit, c1_fit, c2_fit)
    bkg_density = np.clip(bkg_density, 1e-10, 1e10)

    exp_bkg = (np.sum(bkg_density) / float(bins)) * (WINDOW_HIGH - WINDOW_LOW) / bin_width

    O = float(obs_signal)
    B = float(exp_bkg)

    if B < 1e-5:
        return O, B, 0.0

    nll_null = poisson_nll(B, O)
    nll_alt = poisson_nll(O, O)

    q0 = 2.0 * (nll_null - nll_alt)
    z_score = np.sqrt(abs(q0)) if O > B else -np.sqrt(abs(q0))

    return O, B, z_score

if __name__ == '__main__':
    for channel_name, mass_file in independent_channels.items():
        if mass_file not in npy_files:
            print(f"[WARNING] Matrix {mass_file} not found. Skipping channel.")
            continue

        mass_array = np.load(npy_files[mass_file], mmap_mode='r')
        print("--------------------------------------------------")
        print(f" [EVALUATION] {channel_name.upper()} ")
        print(f" Analyzed Events       : {len(mass_array):,}")

        O, B, z_score = evaluate_plr_independent(mass_array, bins=100)

        print(f" Observed Events (O)   : {O:,.0f}")
        print(f" Expected Standard (B) : {B:,.2f}")
        print(f" STATISTICAL Z-SCORE   : {z_score:.2f} \u03c3\n")

    print("==================================================")
    print("[INFO] EXECUTING SIMULTANEOUS PLR FIT (JOINT CHANNELS)...\n")

    if joint_channels["Di-Photon"] in npy_files and joint_channels["4-Lepton"] in npy_files:
        dipho_arr = np.load(npy_files[joint_channels["Di-Photon"]], mmap_mode='r')
        fourl_arr = np.load(npy_files[joint_channels["4-Lepton"]], mmap_mode='r')

        O1, B1, _ = evaluate_plr_independent(dipho_arr, bins=50)
        O2, B2, _ = evaluate_plr_independent(fourl_arr, bins=50)
        print("[KINEMATIC PROFILE]")
        print(f" -> Di-Photon (Channel A) : Obs = {O1:,.0f} | Exp = {B1:.2f}")
        print(f" -> 4-Lepton  (Channel B) : Obs = {O2:,.0f} | Exp = {B2:.2f}")

        if B1 > 0 and B2 > 0:
            def joint_nll_alt(params):
                mu_leak, eps_transfer = params
                Exp1 = B1 - mu_leak
                Exp2 = B2 + (eps_transfer * mu_leak)
                if Exp1 <= 0 or Exp2 <= 0 or mu_leak < 0 or eps_transfer < 0 or eps_transfer > 1.0:
                    return np.inf
                return (Exp1 - O1 * np.log(Exp1)) + (Exp2 - O2 * np.log(Exp2))

            def nll_null():
                return (B1 - O1 * np.log(B1)) + (B2 - O2 * np.log(B2))

            res = minimize(joint_nll_alt, x0=[10.0, 0.1], method='Nelder-Mead', options={'maxiter': 5000})
            mu_leak_fit, eps_fit = res.x
            min_nll_alt = res.fun
            min_nll_null = nll_null()
            q0 = 2.0 * (min_nll_null - min_nll_alt)
            p_val_joint = chi2.sf(abs(q0), df=2)
            p_val_joint = np.clip(p_val_joint, 1e-15, 1.0)
            z_score_joint = norm.ppf(1.0 - (p_val_joint / 2.0))
            print("\n[RESULTS] JOINT PROFILE LIKELIHOOD RATIO")
            print(f" Absolute Event Deficit (\u03bc) : {mu_leak_fit:.2f} events")
            print(f" Transfer Efficiency (\u03b5)    : {eps_fit*100:.2f}%")
            print(f" STATISTICAL Z-SCORE          : {z_score_joint:.2f} \u03c3\n")
        else:
            print("[WARNING] Insufficient expected background to perform Joint PLR.")

    # ==========================================================================
    # MODULE 2 & 3: MET DIVERGENCE & LUMINOSITY SCALING
    # ==========================================================================
    print("==================================================")
    print(" [MODULE 2 & 3] MET DIVERGENCE & LUMINOSITY SCALING ")
    print("==================================================\n")
    FRACTIONS = [0.20, 0.50, 1.00]

    for channel_name, (mass_file, met_file) in met_channels.items():
        if mass_file not in npy_files or met_file not in npy_files:
            print(f"[WARNING] Missing matrices for {channel_name}. Skipping evaluation.")
            continue

        mass_full = np.load(npy_files[mass_file], mmap_mode='r')
        met_full = np.load(npy_files[met_file], mmap_mode='r')
        total_events = len(mass_full)

        print("==================================================")
        print(f" [EVALUATION] {channel_name.upper()} ")
        print("==================================================")

        mask_sr = (mass_full >= WINDOW_LOW) & (mass_full <= WINDOW_HIGH)
        mask_cr = ((mass_full >= 150.0) & (mass_full < WINDOW_LOW)) | ((mass_full > WINDOW_HIGH) & (mass_full <= 250.0))
        met_sr = met_full[mask_sr]
        met_cr = met_full[mask_cr]

        if len(met_sr) > 2 and len(met_cr) > 2:
            ks_stat_full, ks_pvalue_full = ks_2samp(met_sr, met_cr)
            t_stat_full, t_pvalue_full = ttest_ind(met_sr, met_cr, equal_var=False)
            mean_cr = np.mean(met_cr)
            mean_sr = np.mean(met_sr)
            divergence = mean_sr - mean_cr
            print(f" -> MET Control Region : {mean_cr:.2f} GeV")
            print(f" -> MET Signal Region  : {mean_sr:.2f} GeV (Divergence: {divergence:+.2f} GeV)")
            print(f" -> KS-Test Statistic  : {ks_stat_full:.4f} (p: {ks_pvalue_full:.2e})")
            print(f" -> Welch's T-Statistic: {t_stat_full:.2f} (p: {t_pvalue_full:.2e})\n")

            plt.figure(figsize=(10, 6))
            plt.style.use('dark_background')
            max_met = 200 if "JetHT" in channel_name else 150
            plt.hist(met_cr, bins=60, range=(0, max_met), density=True, color='gray', alpha=0.5, label='Control Region')
            plt.hist(met_sr, bins=60, range=(0, max_met), density=True, color='cyan', alpha=0.7, histtype='step', linewidth=2, label='Signal Region (214-216 GeV)')
            plt.axvline(mean_cr, color='gray', linestyle='dashed', linewidth=1)
            plt.axvline(mean_sr, color='cyan', linestyle='dashed', linewidth=2)
            plt.title(f'Kinematic Divergence of E_T^miss ({channel_name})', color='white', pad=15)
            plt.xlabel('Missing Transverse Energy (GeV)', color='white')
            plt.ylabel('Event Density', color='white')
            plt.legend(loc='upper right', framealpha=0.2)
            plt.tight_layout()
            plt.savefig(f"MET_Divergence_{channel_name}_FINAL.png", dpi=300)
            plt.close()

        np.random.seed(42)
        indices = np.arange(total_events)
        np.random.shuffle(indices)
        mass_shuffled = mass_full[indices]
        met_shuffled = met_full[indices]
        history_t, history_ks, labels = [], [], []

        for frac in FRACTIONS:
            limit = int(total_events * frac)
            m_sub = mass_shuffled[:limit]
            met_sub = met_shuffled[:limit]

            m_sr = (m_sub >= WINDOW_LOW) & (m_sub <= WINDOW_HIGH)
            m_cr = ((m_sub >= 150.0) & (m_sub < WINDOW_LOW)) | ((m_sub > WINDOW_HIGH) & (m_sub <= 250.0))

            if np.sum(m_sr) < 2 or np.sum(m_cr) < 2:
                continue

            ks_s, _ = ks_2samp(met_sub[m_sr], met_sub[m_cr])
            t_s, _ = ttest_ind(met_sub[m_sr], met_sub[m_cr], equal_var=False)

            history_t.append(t_s)
            history_ks.append(ks_s)
            labels.append(f"{int(frac*100)}%")
            print(f" [LUMI {int(frac*100):3}%] -> T-Stat: {t_s:5.2f} | KS-Stat: {ks_s:.4f}")

        if len(history_t) == 3:
            fig, ax1 = plt.subplots(figsize=(8, 6))
            plt.style.use('dark_background')
            ax1.set_xlabel('Luminosity Fraction', color='white')
            ax1.set_ylabel("Welch's T-Statistic", color='cyan')
            l1 = ax1.plot(labels, history_t, marker='o', color='cyan', lw=2, label='T-Stat')
            ax1.tick_params(axis='y', labelcolor='cyan')

            ax2 = ax1.twinx()
            ax2.set_ylabel('Kolmogorov-Smirnov Statistic', color='magenta')
            l2 = ax2.plot(labels, history_ks, marker='s', color='magenta', lw=2, ls='--', label='KS-Stat')
            ax2.tick_params(axis='y', labelcolor='magenta')

            ax1.legend(l1 + l2, [l.get_label() for l in l1 + l2], loc='upper left', framealpha=0.2)
            plt.title(f'Statistical Scaling Coherence ({channel_name})', color='white')
            fig.tight_layout()
            plt.savefig(f"Scaling_Evolution_{channel_name}_FINAL.png", dpi=300)
            plt.close(fig)

    print("\n==================================================")
    print(" [SUCCESS] EVALUATOR EXECUTION COMPLETED. ")
    print("==================================================")
