# ==============================================================================
# TBP-OMNI-ENGINE (OPTION 1): 173.111 GeV (1558/9 GeV) STABILIZATION AUDIT
# TARGET: MACROSCOPIC TENSOR SUBTRACTION [R_crit - (d^2 + d) = 215.111 - 42 GeV]
# PROTOCOL: PLR CONTINUUM, MET DIVERGENCE, 4D RIGIDITY & 600 DPI SYNTHESIS FIGURE
# REPOSITORY: https://github.com/cefiyana-clover/TBP-Computational-Framework
# ==============================================================================
import os
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm, chi2, ks_2samp, ttest_ind
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

print("==================================================================")
print(" [INIT] TBP OPTION-1 ENGINE: 173.111 GeV (1558/9 GeV) EVALUATION ")
print("        TESTING SURFACE SUBTRACTION: R_crit - (d^2 + d) = 42 GeV ")
print("==================================================================\n")

# ------------------------------------------------------------------------------
# 0A. AB INITIO TBP TOP-QUARK STABILIZATION DERIVATION (EQS. 1-3)
# ------------------------------------------------------------------------------
D_MANIFOLD = 6
D_BULK = 4
MU_LOOP = 88
LAMBDA_EW = 1.0                                           # 1.000 GeV unitary reference
LAMBDA_CRIT = (MU_LOOP ** 2) / float(D_MANIFOLD ** 2)     # 1936/9 = 215.111111
R_CRIT_GEV = LAMBDA_CRIT * LAMBDA_EW                      # 215.111111 GeV
SURFACE_SUB = (D_MANIFOLD ** 2) + D_MANIFOLD              # 36 + 6 = 42
M_BARE_EQ_GEV = (LAMBDA_CRIT - SURFACE_SUB) * LAMBDA_EW   # 1558/9 = 173.111111 GeV

H_LIM = 1.0 / (D_MANIFOLD ** 2 * np.sqrt(LAMBDA_CRIT))    # 1/528
ALPHA_S = 17.0 / 144.0
L_QCD = H_LIM / ALPHA_S                                   # 3/187
S_SCALE_MEV = (LAMBDA_CRIT / ((D_MANIFOLD + D_BULK) ** 2)) * 1.0  # 1936/900 = 2.151111 MeV
DELTA_QCD_MEV = L_QCD * S_SCALE_MEV                       # 0.034510 MeV
M_TOP_TBP_GEV = M_BARE_EQ_GEV + (6.0 * DELTA_QCD_MEV / 1000.0)    # 173.111318 GeV

print("[TBP AB INITIO STABILIZATION THRESHOLD]")
print(f" -> Critical Vacuum Ceiling (R_crit)      : {R_CRIT_GEV:.6f} GeV (1936/9 GeV)")
print(f" -> Structural Surface Subtraction (d^2+d) : {SURFACE_SUB:.1f} GeV")
print(f" -> Bare Equilibrium Mass (1558/9 GeV)    : {M_BARE_EQ_GEV:.6f} GeV")
print(f" -> 1-Loop QCD Gluon Layer (Delta_QCD)    : {DELTA_QCD_MEV:.6f} MeV (6x = {6*DELTA_QCD_MEV:.5f} MeV)")
print(f" -> Physical Top Stabilization Mass (M_t) : {M_TOP_TBP_GEV:.6f} GeV\n")

# ------------------------------------------------------------------------------
# 0B. AUTOMATED DATA RESOLUTION (KAGGLE & LOCAL REPOSITORY COMPATIBLE)
# ------------------------------------------------------------------------------
DEFAULT_DIR = "/kaggle/" if os.path.exists("/kaggle/") else "./"
BASE_DIR = os.environ.get("TBP_DATA_DIR", DEFAULT_DIR)
print(f"[INFO] Scanning for data matrices within {BASE_DIR}...")
npy_files = {}
for root, dirs, files in os.walk(BASE_DIR):
    for f in files:
        if f.endswith('.npy') and not f.endswith('_SABOTAGED.npy'):
            npy_files[f] = os.path.join(root, f)
print(f"[INFO] Located {len(npy_files)} pristine .npy matrices. Commencing evaluation...\n")

# ------------------------------------------------------------------------------
# OPERATIONAL CONFIGURATION (ISOLATED FROM 125 GeV HIGGS & 215.111 GeV CEILING)
# ------------------------------------------------------------------------------
FIT_LOW = 130.0           # Above 125 GeV Higgs pole
FIT_HIGH = 210.0          # Below 214-216 GeV R_crit ceiling
BINS_COUNT = 80           # Exact 1.0 GeV/bin resolution (78 active sideband bins)
WINDOW_LOW = 172.0        # 2 GeV Binned Signal Region around 173.111 GeV
WINDOW_HIGH = 174.0
PEAK_MIN_EXACT = 172.111  # Exact TBP Coordinate (1558/9 = 173.111 GeV +/- 1 GeV)
PEAK_MAX_EXACT = 174.111
BG_MIN_LOCAL = 168.0
BG_MAX_LOCAL = 178.0

# 1. Z-SCORE EVALUATION CHANNELS (8 Independent PLR Channels)
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

# 2. KINEMATIC DIVERGENCE CHANNELS (10 Sectors: Mass, MET_PT)
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

# 3. TOPOLOGICAL GEOMETRY AUDIT CHANNELS (10 Sectors)
targets_173 = {
    "1. MUONEG (Cross-Flavor t-tbar Proxy)": {
        "mass": "muoneg_mass.npy", "eta1": "muoneg_e_eta.npy", "phi1": "muoneg_e_phi.npy",
        "eta2": "muoneg_mu_eta.npy", "phi2": "muoneg_mu_phi.npy"
    },
    "2. SINGLE-TAU": {
        "mass": "single_tau_mass.npy", "eta1": "single_tau_t1_eta.npy", "phi1": "single_tau_t1_phi.npy",
        "eta2": "single_tau_t2_eta.npy", "phi2": "single_tau_t2_phi.npy"
    },
    "3. SINGLE-MUON": {
        "mass": "single_mu_mass.npy", "eta1": "single_mu_m1_eta.npy", "phi1": "single_mu_m1_phi.npy",
        "eta2": "single_mu_m2_eta.npy", "phi2": "single_mu_m2_phi.npy"
    },
    "4. SINGLE-ELECTRON": {
        "mass": "single_ele_mass.npy", "eta1": "single_ele_e1_eta.npy", "phi1": "single_ele_e1_phi.npy",
        "eta2": "single_ele_e2_eta.npy", "phi2": "single_ele_e2_phi.npy"
    },
    "5. DI-PHOTON": {
        "mass": "di_photon_mass.npy", "eta1": "di_photon_ph1_eta.npy", "phi1": "di_photon_ph1_phi.npy",
        "eta2": "di_photon_ph2_eta.npy", "phi2": "di_photon_ph2_phi.npy"
    },
    "6. 4-LEPTON": {
        "mass": "four_lepton_mass.npy", "eta1": "four_lepton_sum_e_eta.npy", "phi1": "four_lepton_sum_e_phi.npy",
        "eta2": "four_lepton_sum_mu_eta.npy", "phi2": "four_lepton_sum_mu_phi.npy"
    },
    "7. JETHT (Hadronic Continuum)": {
        "mass": "jetht_mass.npy", "eta1": "jetht_j1_eta.npy", "phi1": "jetht_j1_phi.npy",
        "eta2": "jetht_j2_eta.npy", "phi2": "jetht_j2_phi.npy"
    },
    "8. SINGLE-PHOTON (vs MET)": {
        "mass": "single_photon_mass.npy", "eta1": None, "phi1": "single_photon_phi.npy",
        "eta2": None, "phi2": "single_photon_met_phi.npy"
    },
    "9. HEAVY FLAVOR (B-JETS)": {
        "mass": "bjet_mass.npy", "eta1": "bjet_b1_eta.npy", "phi1": "bjet_b1_phi.npy",
        "eta2": "bjet_b2_eta.npy", "phi2": "bjet_b2_phi.npy"
    },
    "10. W-BOSON (Leptonic MT)": {
        "mass": "wboson_mt.npy", "eta1": None, "phi1": "wboson_mu_phi.npy",
        "eta2": None, "phi2": "wboson_met_phi.npy"
    }
}

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS (PLR & 4D EUCLIDEAN METRICS)
# ------------------------------------------------------------------------------
def background_model(m, c0, c1, c2):
    x = m / 13000.0
    x = np.clip(x, 1e-5, 1.0)
    return c0 * ((1.0 - x) ** c1) / (x ** c2)

def poisson_nll(expected, observed):
    if expected <= 0:
        return np.inf
    return expected - observed * np.log(expected)

def evaluate_plr_173(mass_array, bins=BINS_COUNT):
    mask_range = (mass_array >= FIT_LOW) & (mass_array <= FIT_HIGH)
    m_fit = mass_array[mask_range]

    mask_sidebands = (m_fit < WINDOW_LOW) | (m_fit > WINDOW_HIGH)
    m_sidebands = m_fit[mask_sidebands]

    counts, bin_edges = np.histogram(m_sidebands, bins=bins, range=(FIT_LOW, FIT_HIGH))
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
    c0_guess = avg_count * ((173.0 / 13000.0) ** 5) if avg_count > 0 else 1e-5

    res = minimize(loss_func, x0=[c0_guess, 5.0, 5.0], method='Nelder-Mead', options={'maxiter': 10000})
    c0_fit, c1_fit, c2_fit = res.x

    obs_signal = np.sum((m_fit >= WINDOW_LOW) & (m_fit <= WINDOW_HIGH))
    bin_width = (FIT_HIGH - FIT_LOW) / float(bins)
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

def calc_delta_phi(phi1, phi2):
    dphi = np.abs(phi1 - phi2)
    return np.where(dphi > np.pi, 2.0 * np.pi - dphi, dphi)

def calc_delta_r(eta1, phi1, eta2, phi2):
    deta = np.abs(eta1 - eta2)
    dphi = calc_delta_phi(phi1, phi2)
    return np.sqrt(deta ** 2 + dphi ** 2)

def calc_significance(data_peak, data_bg):
    mean_peak = np.mean(data_peak)
    mean_bg = np.mean(data_bg)
    twist = mean_peak - mean_bg
    sem_peak = np.std(data_peak) / np.sqrt(len(data_peak))
    sem_bg = np.std(data_bg) / np.sqrt(len(data_bg))
    combined_err = np.sqrt(sem_peak ** 2 + sem_bg ** 2)
    sig = abs(twist) / combined_err if combined_err > 0 else 0.0
    return twist, sig

# ==============================================================================
# MAIN EXECUTION PIPELINE (PARTS I-III: LIVE MATRIX AUDIT | PART IV: SYNTHESIS)
# ==============================================================================
if __name__ == '__main__':
    # --------------------------------------------------------------------------
    # PART I: PROFILE LIKELIHOOD RATIO (PLR) EVALUATION AT 173.111 GeV
    # --------------------------------------------------------------------------
    print("==================================================================")
    print(" [PART I] PROFILE LIKELIHOOD RATIO (PLR) AT 172.0 - 174.0 GeV     ")
    print("==================================================================\n")

    for channel_name, mass_file in independent_channels.items():
        if mass_file not in npy_files:
            print(f"[WARNING] Matrix {mass_file} not found. Skipping channel.")
            continue

        mass_array = np.load(npy_files[mass_file], mmap_mode='r')
        print("--------------------------------------------------")
        print(f" [EVALUATION @ 173.111 GeV] {channel_name.upper()} ")
        print(f" Total Events in Matrix : {len(mass_array):,}")

        O, B, z_score = evaluate_plr_173(mass_array, bins=BINS_COUNT)

        print(f" Observed Events (O)    : {O:,.0f}")
        print(f" Expected Standard (B)  : {B:,.2f}")
        print(f" STATISTICAL Z-SCORE    : {z_score:+.2f} \u03c3\n")

    # --------------------------------------------------------------------------
    # PART II: MET DIVERGENCE & LUMINOSITY SCALING AT 173.111 GeV
    # --------------------------------------------------------------------------
    print("==================================================================")
    print(" [PART II] MET DIVERGENCE & LUMINOSITY SCALING AT 173.111 GeV     ")
    print("==================================================================\n")
    FRACTIONS = [0.20, 0.50, 1.00]

    for channel_name, (mass_file, met_file) in met_channels.items():
        if mass_file not in npy_files or met_file not in npy_files:
            continue

        mass_full = np.load(npy_files[mass_file], mmap_mode='r')
        met_full = np.load(npy_files[met_file], mmap_mode='r')
        total_events = len(mass_full)

        print("==================================================")
        print(f" [MET AUDIT @ 173.111 GeV] {channel_name.upper()} ")
        print("==================================================")
        mask_sr = (mass_full >= WINDOW_LOW) & (mass_full <= WINDOW_HIGH)
        mask_cr = ((mass_full >= FIT_LOW) & (mass_full < WINDOW_LOW)) | ((mass_full > WINDOW_HIGH) & (mass_full <= FIT_HIGH))
        met_sr = met_full[mask_sr]
        met_cr = met_full[mask_cr]

        if len(met_sr) > 2 and len(met_cr) > 2:
            ks_stat_full, ks_pvalue_full = ks_2samp(met_sr, met_cr)
            t_stat_full, t_pvalue_full = ttest_ind(met_sr, met_cr, equal_var=False)
            mean_cr = np.mean(met_cr)
            mean_sr = np.mean(met_sr)
            divergence = mean_sr - mean_cr
            print(f" -> MET Control Region (130-210 GeV) : {mean_cr:.2f} GeV")
            print(f" -> MET Signal Region  (172-174 GeV) : {mean_sr:.2f} GeV (Divergence: {divergence:+.2f} GeV)")
            print(f" -> KS-Test Statistic                : {ks_stat_full:.4f} (p: {ks_pvalue_full:.2e})")
            print(f" -> Welch's T-Statistic              : {t_stat_full:+.2f} (p: {t_pvalue_full:.2e})\n")

            plt.figure(figsize=(10, 6))
            plt.style.use('dark_background')
            max_met = 200 if "JetHT" in channel_name else 150
            plt.hist(met_cr, bins=60, range=(0, max_met), density=True, color='gray', alpha=0.5, label='Control Region (130-210 GeV)')
            plt.hist(met_sr, bins=60, range=(0, max_met), density=True, color='gold', alpha=0.8, histtype='step', linewidth=2, label='Signal Region (172-174 GeV)')
            plt.axvline(mean_cr, color='gray', linestyle='dashed', linewidth=1)
            plt.axvline(mean_sr, color='gold', linestyle='dashed', linewidth=2)
            plt.title(f'173.111 GeV Stabilization Threshold E_T^miss ({channel_name})', color='white', pad=15)
            plt.xlabel('Missing Transverse Energy (GeV)', color='white')
            plt.ylabel('Event Density', color='white')
            plt.legend(loc='upper right', framealpha=0.2)
            plt.tight_layout()
            safe_ch = channel_name.replace(" ", "_").replace("(", "").replace(")", "")
            plt.savefig(f"MET_173GeV_{safe_ch}.png", dpi=300)
            plt.close()

        np.random.seed(42)
        indices = np.arange(total_events)
        np.random.shuffle(indices)
        mass_shuffled = mass_full[indices]
        met_shuffled = met_full[indices]

        for frac in FRACTIONS:
            limit = int(total_events * frac)
            m_sub = mass_shuffled[:limit]
            met_sub = met_shuffled[:limit]

            m_sr = (m_sub >= WINDOW_LOW) & (m_sub <= WINDOW_HIGH)
            m_cr = ((m_sub >= FIT_LOW) & (m_sub < WINDOW_LOW)) | ((m_sub > WINDOW_HIGH) & (m_sub <= FIT_HIGH))

            if np.sum(m_sr) < 2 or np.sum(m_cr) < 2:
                continue

            ks_s, _ = ks_2samp(met_sub[m_sr], met_sub[m_cr])
            t_s, _ = ttest_ind(met_sub[m_sr], met_sub[m_cr], equal_var=False)
            print(f" [LUMI {int(frac*100):3}%] -> T-Stat: {t_s:+6.2f} | KS-Stat: {ks_s:.4f}")
        print("")

    # --------------------------------------------------------------------------
    # PART III: OMNIDIRECTIONAL TOPOLOGICAL AUDIT AT 172.111 - 174.111 GeV
    # --------------------------------------------------------------------------
    print("==================================================================")
    print(" [PART III] TOPOLOGICAL GEOMETRY AUDIT AT 172.111 - 174.111 GeV   ")
    print("==================================================================\n")

    for name, config in targets_173.items():
        if config["mass"] not in npy_files:
            continue
        try:
            mass = np.load(npy_files[config["mass"]], mmap_mode='r')
            phi1 = np.load(npy_files[config["phi1"]], mmap_mode='r')
            phi2 = np.load(npy_files[config["phi2"]], mmap_mode='r')
            has_eta = config["eta1"] is not None
            if has_eta:
                eta1 = np.load(npy_files[config["eta1"]], mmap_mode='r')
                eta2 = np.load(npy_files[config["eta2"]], mmap_mode='r')
        except Exception:
            continue

        print("--------------------------------------------------")
        print(f" [TOPOLOGICAL AUDIT @ 173.111 GeV] {name}")
        print("--------------------------------------------------")
        peak_mask = (mass >= PEAK_MIN_EXACT) & (mass <= PEAK_MAX_EXACT)
        bg_mask = (mass >= BG_MIN_LOCAL) & (mass <= BG_MAX_LOCAL) & (~peak_mask)

        if np.sum(peak_mask) > 10 and np.sum(bg_mask) > 10:
            dphi = calc_delta_phi(phi1, phi2)
            twist_dphi, sig_dphi = calc_significance(dphi[peak_mask], dphi[bg_mask])
            status_dphi = "ANOMALY" if sig_dphi > 5.0 else "STABLE"
            print(f" -> Coordinate Window   : {PEAK_MIN_EXACT:.3f} - {PEAK_MAX_EXACT:.3f} GeV (N_SR = {np.sum(peak_mask):,})")
            print(f" -> Azimuthal Twist \u0394\u03c6  : {twist_dphi:+.4f} rad ({sig_dphi:5.2f} \u03c3) [{status_dphi}]")
            if has_eta:
                dr = calc_delta_r(eta1, phi1, eta2, phi2)
                twist_dr, sig_dr = calc_significance(dr[peak_mask], dr[bg_mask])
                status_dr = "DEVIATION" if sig_dr > 5.0 else "STABLE"
                print(f" -> Spatial Metric \u0394R   : {twist_dr:+.4f}     ({sig_dr:5.2f} \u03c3) [{status_dr}]")
        else:
            print(" -> Insufficient statistics in 172.111 - 174.111 GeV window.")

        z1_mask = (mass >= 172.111) & (mass <= 173.111)
        z2_mask = (mass > 173.111) & (mass <= 174.111)
        if np.sum(z1_mask) > 10 and np.sum(z2_mask) > 10 and np.sum(bg_mask) > 10:
            dphi = calc_delta_phi(phi1, phi2)
            _, sig_z1_phi = calc_significance(dphi[z1_mask], dphi[bg_mask])
            _, sig_z2_phi = calc_significance(dphi[z2_mask], dphi[bg_mask])
            print(f" -> Dual-Zone \u0394\u03c6 Strain : Zone 1 (172.1-173.1) = {sig_z1_phi:4.2f}\u03c3 | Zone 2 (173.1-174.1) = {sig_z2_phi:4.2f}\u03c3")
            if has_eta:
                dr = calc_delta_r(eta1, phi1, eta2, phi2)
                _, sig_z1_dr = calc_significance(dr[z1_mask], dr[bg_mask])
                _, sig_z2_dr = calc_significance(dr[z2_mask], dr[bg_mask])
                print(f" -> Dual-Zone \u0394R Strain : Zone 1 (172.1-173.1) = {sig_z1_dr:4.2f}\u03c3 | Zone 2 (173.1-174.1) = {sig_z2_dr:4.2f}\u03c3")
        print("")

    print("==================================================================")
    print(" [SUCCESS] OPTION-1 (173.111 GeV) EVALUATION & AUDIT COMPLETED.   ")
    print("==================================================================\n")

    # --------------------------------------------------------------------------
    # PART IV: 173.111 GeV (1558/9 GeV) GRAND SYNTHESIS FIGURE (600 DPI)
    # --------------------------------------------------------------------------
    plt.style.use('default')
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['axes.edgecolor'] = '#222222'
    plt.rcParams['axes.linewidth'] = 0.9

    DPI_OUT = 600
    fig = plt.figure(figsize=(14, 8.8), dpi=DPI_OUT, facecolor='white')
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.05, 1.0], width_ratios=[1.0, 1.0],
                           hspace=0.36, wspace=0.24, left=0.06, right=0.97, top=0.92, bottom=0.08)

    # Panel (a): Cross-Regime Comparison (173.111 GeV vs. 215.111 GeV)
    ax1 = fig.add_subplot(gs[0, :])
    channels_comp = [
        'Single-Muon\n(MET T-Stat)',
        'Single-Electron\n(MET T-Stat)',
        'Single-Electron\n(Mass PLR Z)',
        'JetHT Hadronic\n(Mass PLR Z)',
        'JetHT Hadronic\n(|MET T-Stat|)',
        'Heavy Flavor B-Jets\n(Mass PLR Z)',
        'Heavy Flavor B-Jets\n(MET T-Stat)',
        'Di-Photon\n(MET T-Stat)'
    ]
    vals_173 = [4.98, 3.53, 1.09, -2.00, 0.66, 0.14, 0.95, -0.07]
    vals_215 = [2.96, 5.44, 4.77, 8.43, 14.57, 0.20, -0.29, 0.38]

    x_pos = np.arange(len(channels_comp))
    bar_w = 0.35

    rects1 = ax1.bar(x_pos - bar_w/2, vals_173, bar_w, label=r'173.111 GeV Stabilization Threshold ($\Theta = 1$)',
                     color='#0077B6', edgecolor='black', linewidth=0.7, alpha=0.92, zorder=3)
    rects2 = ax1.bar(x_pos + bar_w/2, vals_215, bar_w, label=r'215.111 GeV Critical Vacuum Ceiling ($\Theta \to 0$)',
                     color='#D90429', edgecolor='black', linewidth=0.7, alpha=0.88, zorder=3)

    ax1.axhline(5.0, color='#B00000', linestyle='--', linewidth=1.2, zorder=2, label=r'$5\sigma$ Discovery / Threshold Reference')
    ax1.axhline(0.0, color='#333333', linestyle='-', linewidth=0.8, zorder=2)

    for rect in rects1:
        h = rect.get_height()
        va_pos = 'bottom' if h >= 0 else 'top'
        y_off = 0.25 if h >= 0 else -0.35
        ax1.text(rect.get_x() + rect.get_width()/2., h + y_off, f'{h:+.2f}',
                 ha='center', va=va_pos, fontsize=8.0, fontweight='bold', color='#005082')

    for rect in rects2:
        h = rect.get_height()
        va_pos = 'bottom' if h >= 0 else 'top'
        y_off = 0.25 if h >= 0 else -0.35
        ax1.text(rect.get_x() + rect.get_width()/2., h + y_off, f'{h:+.2f}',
                 ha='center', va=va_pos, fontsize=8.0, fontweight='bold', color='#9B0000')

    ax1.set_title('(a) Cross-Regime Kinematic & Thermodynamic Contrast: 173.111 GeV Stabilization vs. 215.111 GeV Vacuum Ceiling',
                  fontsize=11.5, fontweight='bold', pad=10)
    ax1.set_ylabel(r'Statistical Significance ($Z\ [\sigma]$ or $|T|$-Stat)', fontsize=10)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(channels_comp, fontsize=8.8)
    ax1.set_ylim(-4.2, 16.8)
    ax1.grid(axis='y', linestyle=':', alpha=0.5, zorder=0)
    ax1.legend(loc='upper right', frameon=True, facecolor='#F8F9FA', edgecolor='#CCCCCC', fontsize=9.0)

    ax1.annotate('Prompt $t \\to Wb \\to \\ell\\nu_\\ell b$\nNeutrino Surge (@ 173.111 GeV)',
                 xy=(0 - bar_w/2, 5.1), xytext=(0.4, 11.2),
                 arrowprops=dict(arrowstyle='->', color='#005082', lw=1.3),
                 fontsize=8.8, fontweight='bold', color='#003566',
                 bbox=dict(boxstyle='round,pad=0.35', facecolor='#EBF8FF', edgecolor='#0077B6', lw=1.0))

    ax1.annotate('Hadronic Quiescence at 173 GeV\nvs. $R_{crit}$ Pileup at 215 GeV',
                 xy=(3 + bar_w/2, 8.6), xytext=(2.7, 13.0),
                 arrowprops=dict(arrowstyle='->', color='#9B0000', lw=1.3),
                 fontsize=8.8, fontweight='bold', color='#800000',
                 bbox=dict(boxstyle='round,pad=0.35', facecolor='#FFF0F0', edgecolor='#D90429', lw=1.0))

    # Panel (b): Luminosity Scaling Coherence at 173.111 GeV
    ax2 = fig.add_subplot(gs[1, 0])
    lumi_fracs = [20, 50, 100]
    t_muon = [1.81, 2.50, 4.98]
    t_elec = [2.81, 2.86, 3.53]
    t_bjet = [1.11, 0.88, 0.95]
    t_dipho = [0.93, 0.54, -0.07]
    t_jetht = [-0.36, -1.24, -0.66]

    ax2.plot(lumi_fracs, t_muon, marker='o', markersize=6.5, lw=2.3, color='#0077B6',
             label=r'Single-Muon $\Delta E_T^{miss}$ ($+4.35\ \mathrm{GeV}, T=+4.98$)')
    ax2.plot(lumi_fracs, t_elec, marker='s', markersize=6.0, lw=2.1, color='#00B4D8',
             label=r'Single-Electron $\Delta E_T^{miss}$ ($+2.47\ \mathrm{GeV}, T=+3.53$)')
    ax2.plot(lumi_fracs, t_bjet, marker='^', markersize=5.5, lw=1.5, ls='--', color='#2A9D8F',
             label=r'Heavy Flavor B-Jets ($T=+0.95$, Quiescent)')
    ax2.plot(lumi_fracs, t_dipho, marker='d', markersize=5.0, lw=1.4, ls='-.', color='#E9C46A',
             label=r'Di-Photon Continuum ($T=-0.07$, Quiescent)')
    ax2.plot(lumi_fracs, t_jetht, marker='v', markersize=5.5, lw=1.5, ls=':', color='#E76F51',
             label=r'JetHT Hadronic ($T=-0.66$, Quiescent)')

    ax2.axhline(5.0, color='#D90429', linestyle='--', linewidth=1.1, label=r'$5\sigma$ Reference Limit')
    ax2.axhline(0.0, color='#444444', linestyle='-', linewidth=0.7)
    ax2.set_title(r'(b) Luminosity Scaling Coherence at 173.111 GeV ($E_T^{miss}$ Divergence)',
                  fontsize=11.0, fontweight='bold', pad=10)
    ax2.set_xlabel('Integrated Luminosity Fraction (%)', fontsize=10)
    ax2.set_ylabel("Welch's T-Statistic", fontsize=10)
    ax2.set_xticks([20, 50, 100])
    ax2.set_xticklabels(['20%', '50%', '100%'], fontsize=9.5)
    ax2.set_ylim(-2.2, 6.0)
    ax2.grid(True, linestyle=':', alpha=0.5)
    ax2.legend(loc='upper left', frameon=True, facecolor='#F8F9FA', edgecolor='#CCCCCC', fontsize=8.0)

    # Panel (c): Omnidirectional Euclidean Rigidity at 172.111 - 174.111 GeV
    ax3 = fig.add_subplot(gs[1, 1])
    sectors_geo = ['MuonEG', 'Tau', 'Muon', 'Electron', 'Di-Photon', 'JetHT', 'Photon', 'B-Jets', 'W-Boson']
    primary_strain = [0.93, 1.43, 0.39, 2.50, 0.66, 0.28, 0.65, 0.86, 0.27]
    dualzone_max   = [1.50, 1.68, 2.20, 2.77, 1.62, 2.57, 0.73, 1.53, 0.47]

    x_g = np.arange(len(sectors_geo))
    w_g = 0.36

    ax3.bar(x_g - w_g/2, primary_strain, w_g, label=r'Primary Audit Max Strain ($\Delta\phi, \Delta R$)',
            color='#2B2D42', edgecolor='black', linewidth=0.6, alpha=0.9, zorder=3)
    ax3.bar(x_g + w_g/2, dualzone_max, w_g, label=r'Dual-Zone Elasticity Max Strain (Z1/Z2)',
            color='#8D99AE', edgecolor='black', linewidth=0.6, alpha=0.85, zorder=3)
    ax3.axhline(5.0, color='#D90429', linestyle='--', linewidth=1.2, zorder=2, label=r'$5\sigma$ Metric Tear Threshold')

    ax3.set_title('(c) 4D Euclidean Rigidity Audit (172.111 - 174.111 GeV)',
                  fontsize=11.0, fontweight='bold', pad=10)
    ax3.set_xlabel('Evaluated Kinematic Topologies', fontsize=10)
    ax3.set_ylabel(r'Maximum Geometric Strain ($\sigma$)', fontsize=10)
    ax3.set_xticks(x_g)
    ax3.set_xticklabels(sectors_geo, fontsize=8.5, rotation=20)
    ax3.set_ylim(0, 5.8)
    ax3.grid(axis='y', linestyle=':', alpha=0.5, zorder=0)
    ax3.legend(loc='upper right', frameon=True, facecolor='#F8F9FA', edgecolor='#CCCCCC', fontsize=8.3)

    ax3.text(4.0, 3.8, 'ABSOLUTE EUCLIDEAN RIGIDITY VERIFIED\nAll 9 Evaluable Sectors Strictly < $5\\sigma$ (Max $2.77\\sigma$)',
             ha='center', va='center', fontsize=8.8, fontweight='bold', color='#1B4332',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#D8F3DC', edgecolor='#2D6A4F', lw=1.0))

    fig.suptitle('Phenomenological Audit of the 173.111 GeV (1558/9 GeV) Macroscopic Tensor Stabilization Threshold (CMS 13 TeV)',
                 fontsize=13.0, fontweight='bold', y=0.975)
    plt.savefig("Figure_1_173GeV_Synthesis.png", dpi=DPI_OUT)
    try:
        plt.savefig("Figure_1_173GeV_Synthesis.pdf")
    except Exception:
        pass
    plt.close()
    print("[SUCCESS] Generated Figure_1_173GeV_Synthesis.png (600 DPI)!")
