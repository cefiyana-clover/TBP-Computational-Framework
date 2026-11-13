# ==============================================================================
# TBP-OMNI-AUDITOR: UNIFIED TOPOLOGICAL INFERENCE ENGINE (215.111 GeV CEILING)
# PROTOCOL: KINEMATIC GEOMETRY, NULL-HYPOTHESIS, SYSTEMATIC CALIBRATION, & SCALING
# REPOSITORY: https://github.com/cefiyana-clover/TBP-Computational-Framework
# ==============================================================================
import os
import numpy as np

print("==================================================")
print(" [INIT] OMNI-AUDITOR: TOPOLOGICAL INFERENCE ENGINE ")
print("==================================================\n")

# ------------------------------------------------------------------------------
# 0. AUTOMATED DATA RESOLUTION (KAGGLE & LOCAL REPOSITORY COMPATIBLE)
# ------------------------------------------------------------------------------
DEFAULT_DIR = "/kaggle/" if os.path.exists("/kaggle/") else "./"
BASE_DIR = os.environ.get("TBP_DATA_DIR", DEFAULT_DIR)
print(f"[INFO] Scanning for data matrices within {BASE_DIR}...")

npy_files = {}
for root, dirs, files in os.walk(BASE_DIR):
    if files:
        for f in files:
            if f.endswith('.npy'):
                npy_files[f] = os.path.join(root, f)

print(f"[INFO] Located {len(npy_files)} .npy matrices. Commencing geometric audit...\n")

# ------------------------------------------------------------------------------
# 1. MATHEMATICAL TENSOR FUNCTIONS & STATISTICAL ERROR
# ------------------------------------------------------------------------------
def calc_delta_phi(phi1, phi2):
    dphi = np.abs(phi1 - phi2)
    dphi = np.where(dphi > np.pi, 2.0 * np.pi - dphi, dphi)
    return dphi

def calc_delta_r(eta1, phi1, eta2, phi2):
    deta = np.abs(eta1 - eta2)
    dphi = calc_delta_phi(phi1, phi2)
    return np.sqrt(deta**2 + dphi**2)

def calc_significance(data_peak, data_bg):
    mean_peak = np.mean(data_peak)
    mean_bg = np.mean(data_bg)
    twist = mean_peak - mean_bg

    sem_peak = np.std(data_peak) / np.sqrt(len(data_peak))
    sem_bg = np.std(data_bg) / np.sqrt(len(data_bg))

    combined_err = np.sqrt(sem_peak**2 + sem_bg**2)
    sig = abs(twist) / combined_err if combined_err > 0 else 0.0
    return twist, sig

# ------------------------------------------------------------------------------
# 2. OPERATIONAL CONFIGURATION (10 SECTORS) - STANDARDIZED TBP BOUNDARY
# ------------------------------------------------------------------------------
targets = {
    "1. MUONEG (Cross-Flavor)": {
        "mass": "muoneg_mass.npy", "eta1": "muoneg_e_eta.npy", "phi1": "muoneg_e_phi.npy",
        "eta2": "muoneg_mu_eta.npy", "phi2": "muoneg_mu_phi.npy",
        "peak_min": 214.111, "peak_max": 216.111
    },
    "2. SINGLE-TAU": {
        "mass": "single_tau_mass.npy", "eta1": "single_tau_t1_eta.npy", "phi1": "single_tau_t1_phi.npy",
        "eta2": "single_tau_t2_eta.npy", "phi2": "single_tau_t2_phi.npy",
        "peak_min": 214.111, "peak_max": 216.111
    },
    "3. SINGLE-MUON": {
        "mass": "single_mu_mass.npy", "eta1": "single_mu_m1_eta.npy", "phi1": "single_mu_m1_phi.npy",
        "eta2": "single_mu_m2_eta.npy", "phi2": "single_mu_m2_phi.npy",
        "peak_min": 214.111, "peak_max": 216.111
    },
    "4. SINGLE-ELECTRON": {
        "mass": "single_ele_mass.npy", "eta1": "single_ele_e1_eta.npy", "phi1": "single_ele_e1_phi.npy",
        "eta2": "single_ele_e2_eta.npy", "phi2": "single_ele_e2_phi.npy",
        "peak_min": 214.111, "peak_max": 216.111
    },
    "5. DI-PHOTON": {
        "mass": "di_photon_mass.npy", "eta1": "di_photon_ph1_eta.npy", "phi1": "di_photon_ph1_phi.npy",
        "eta2": "di_photon_ph2_eta.npy", "phi2": "di_photon_ph2_phi.npy",
        "peak_min": 214.111, "peak_max": 216.111
    },
    "6. 4-LEPTON": {
        "mass": "four_lepton_mass.npy", "eta1": "four_lepton_sum_e_eta.npy", "phi1": "four_lepton_sum_e_phi.npy",
        "eta2": "four_lepton_sum_mu_eta.npy", "phi2": "four_lepton_sum_mu_phi.npy",
        "peak_min": 214.111, "peak_max": 216.111
    },
    "7. JETHT (Hadronic)": {
        "mass": "jetht_mass.npy", "eta1": "jetht_j1_eta.npy", "phi1": "jetht_j1_phi.npy",
        "eta2": "jetht_j2_eta.npy", "phi2": "jetht_j2_phi.npy",
        "peak_min": 214.111, "peak_max": 216.111
    },
    "8. SINGLE-PHOTON (vs MET)": {
        "mass": "single_photon_mass.npy", "eta1": None, "phi1": "single_photon_phi.npy",
        "eta2": None, "phi2": "single_photon_met_phi.npy",
        "peak_min": 214.111, "peak_max": 216.111
    },
    "9. HEAVY FLAVOR (B-JETS)": {
        "mass": "bjet_mass.npy", "eta1": "bjet_b1_eta.npy", "phi1": "bjet_b1_phi.npy",
        "eta2": "bjet_b2_eta.npy", "phi2": "bjet_b2_phi.npy",
        "peak_min": 214.111, "peak_max": 216.111
    },
    "10. W-BOSON (Leptonic MT)": {
        "mass": "wboson_mt.npy", "eta1": None, "phi1": "wboson_mu_phi.npy",
        "eta2": None, "phi2": "wboson_met_phi.npy",
        "peak_min": 214.111, "peak_max": 216.111
    }
}

# ------------------------------------------------------------------------------
# 3. UNIFIED EXECUTION ENGINE
# ------------------------------------------------------------------------------
if __name__ == '__main__':
    for name, config in targets.items():
        if config["mass"] not in npy_files:
            print(f"[WARNING] Matrices for {name} missing. Skipping sector.")
            continue

        try:
            mass = np.load(npy_files[config["mass"]], mmap_mode='r')
            phi1 = np.load(npy_files[config["phi1"]], mmap_mode='r')
            phi2 = np.load(npy_files[config["phi2"]], mmap_mode='r')

            has_eta = config["eta1"] is not None
            if has_eta:
                eta1 = np.load(npy_files[config["eta1"]], mmap_mode='r')
                eta2 = np.load(npy_files[config["eta2"]], mmap_mode='r')
        except Exception as e:
            print(f"[ERROR] Matrix integrity failure on {name}: {e}")
            continue

        print("\n==================================================")
        print(f" [SECTOR] {name} ")
        print("==================================================")

        # MODULE A: PRIMARY ANGULAR TENSOR AUDIT (EPI-CENTER)
        print("\n[MODULE A] PRIMARY ANGULAR TENSOR AUDIT")
        peak_mask = (mass >= config["peak_min"]) & (mass <= config["peak_max"])
        bg_mask = (mass >= 210.0) & (mass <= 220.0) & (~peak_mask)

        if np.sum(peak_mask) > 10 and np.sum(bg_mask) > 10:
            dphi = calc_delta_phi(phi1, phi2)
            twist_dphi, sig_dphi = calc_significance(dphi[peak_mask], dphi[bg_mask])

            status_dphi = "ANOMALY" if sig_dphi > 5.0 else "STABLE"
            print(f" -> Vacuum Coordinates : {config['peak_min']:.3f} - {config['peak_max']:.3f} GeV")
            print(f" -> Azimuthal Twist    : {twist_dphi:+.4f} rad ({sig_dphi:5.2f} \u03c3) [{status_dphi}]")

            if has_eta:
                dr = calc_delta_r(eta1, phi1, eta2, phi2)
                twist_dr, sig_dr = calc_significance(dr[peak_mask], dr[bg_mask])
                status_dr = "DEVIATION" if sig_dr > 5.0 else "STABLE"
                print(f" -> Spatial Metric Tear: {twist_dr:+.4f}     ({sig_dr:5.2f} \u03c3) [{status_dr}]")
            else:
                print(" -> Spatial Metric Tear: [N/A - Transverse Kinematics]")
        else:
            print(" -> Insufficient statistics for Primary Audit.")

        # MODULE B: DUAL-ZONE ELASTICITY EVALUATION
        print("\n[MODULE B] DUAL-ZONE ELASTICITY EVALUATION")
        z1_mask = (mass >= 214.111) & (mass <= 215.111)
        z2_mask = (mass > 215.111) & (mass <= 216.111)
        bg_dual_mask = (mass >= 210.0) & (mass <= 220.0) & (~z1_mask) & (~z2_mask)

        if np.sum(z1_mask) > 10 and np.sum(z2_mask) > 10 and np.sum(bg_dual_mask) > 10:
            dphi = calc_delta_phi(phi1, phi2)
            twist_z1_phi, sig_z1_phi = calc_significance(dphi[z1_mask], dphi[bg_dual_mask])
            twist_z2_phi, sig_z2_phi = calc_significance(dphi[z2_mask], dphi[bg_dual_mask])

            print(f" -> ZONE 1 (Redshift) 214.1-215.1 GeV | \u0394\u03c6: {twist_z1_phi:+.4f} rad ({sig_z1_phi:5.2f} \u03c3)")
            print(f" -> ZONE 2 (Overflow) 215.1-216.1 GeV | \u0394\u03c6: {twist_z2_phi:+.4f} rad ({sig_z2_phi:5.2f} \u03c3)")

            if has_eta:
                dr = calc_delta_r(eta1, phi1, eta2, phi2)
                twist_z1_dr, sig_z1_dr = calc_significance(dr[z1_mask], dr[bg_dual_mask])
                twist_z2_dr, sig_z2_dr = calc_significance(dr[z2_mask], dr[bg_dual_mask])
                print(f" -> ZONE 1 (Redshift) 214.1-215.1 GeV | \u0394R: {twist_z1_dr:+.4f}     ({sig_z1_dr:5.2f} \u03c3)")
                print(f" -> ZONE 2 (Overflow) 215.1-216.1 GeV | \u0394R: {twist_z2_dr:+.4f}     ({sig_z2_dr:5.2f} \u03c3)")
        else:
            print(" -> Insufficient statistics for Dual-Zone Evaluation.")

        # MODULE C: NULL-HYPOTHESIS CALIBRATION (Z, HIGGS, HIGH-MASS)
        print("\n[MODULE C] NULL-HYPOTHESIS CALIBRATION")
        test_zones = {
            "Z-BOSON (91 GeV)": {"peak": (90.0, 92.0), "bg": (85.0, 97.0)},
            "HIGGS BOSON (125 GeV)": {"peak": (124.0, 126.0), "bg": (119.0, 131.0)},
            "HIGH-MASS (300 GeV)": {"peak": (299.0, 301.0), "bg": (294.0, 306.0)}
        }

        for zone_name, bounds in test_zones.items():
            nh_peak_mask = (mass >= bounds["peak"][0]) & (mass <= bounds["peak"][1])
            nh_bg_mask = (mass >= bounds["bg"][0]) & (mass <= bounds["bg"][1]) & (~nh_peak_mask)

            if np.sum(nh_peak_mask) < 10 or np.sum(nh_bg_mask) < 10:
                print(f" -> {zone_name:21} | Insufficient Data")
                continue

            dphi = calc_delta_phi(phi1, phi2)
            twist_nh_phi, sig_nh_phi = calc_significance(dphi[nh_peak_mask], dphi[nh_bg_mask])
            status_nh_phi = "ANOMALY" if sig_nh_phi > 5.0 else "STABLE"

            if has_eta:
                dr = calc_delta_r(eta1, phi1, eta2, phi2)
                twist_nh_dr, sig_nh_dr = calc_significance(dr[nh_peak_mask], dr[nh_bg_mask])
                status_nh_dr = "ANOMALY" if sig_nh_dr > 5.0 else "STABLE"
                print(f" -> {zone_name:21} | \u0394\u03c6: {twist_nh_phi:+.4f} ({sig_nh_phi:4.1f}\u03c3) [{status_nh_phi:7}] | \u0394R: {twist_nh_dr:+.4f} ({sig_nh_dr:4.1f}\u03c3) [{status_nh_dr}]")
            else:
                print(f" -> {zone_name:21} | \u0394\u03c6: {twist_nh_phi:+.4f} ({sig_nh_phi:4.1f}\u03c3) [{status_nh_phi:7}]")

        # MODULE D: SYSTEMATIC KINEMATIC PERTURBATION (JES +/- 5%)
        print("\n[MODULE D] SYSTEMATIC KINEMATIC PERTURBATION (\u00b1 5%)")
        shifts = {"NOMINAL (+0%)": 1.00, "SCALE UP (+5%)": 1.05, "SCALE DOWN (-5%)": 0.95}
        WIN_MIN, WIN_MAX = 214.111, 216.111

        for shift_name, factor in shifts.items():
            mass_shifted = mass * factor
            sys_bg_mask = (mass_shifted >= 210.0) & (mass_shifted <= 220.0) & ((mass_shifted < WIN_MIN) | (mass_shifted > WIN_MAX))
            sys_peak_mask = (mass_shifted >= WIN_MIN) & (mass_shifted <= WIN_MAX)

            if np.sum(sys_bg_mask) < 10 or np.sum(sys_peak_mask) < 10:
                print(f" -> {shift_name:17} | Insufficient data after shift.")
                continue

            dphi = calc_delta_phi(phi1, phi2)
            sys_twist_phi, sys_sig_phi = calc_significance(dphi[sys_peak_mask], dphi[sys_bg_mask])
            sys_status_phi = "ANOMALY" if sys_sig_phi > 5.0 else "STABLE"

            if has_eta:
                dr = calc_delta_r(eta1, phi1, eta2, phi2)
                sys_twist_dr, sys_sig_dr = calc_significance(dr[sys_peak_mask], dr[sys_bg_mask])
                sys_status_dr = "DEVIATION" if sys_sig_dr > 5.0 else "STABLE"
                print(f" -> {shift_name:17} | \u0394\u03c6: {sys_twist_phi:+.4f} ({sys_sig_phi:4.1f}\u03c3) [{sys_status_phi:7}] | \u0394R: {sys_twist_dr:+.4f} ({sys_sig_dr:4.1f}\u03c3) [{sys_status_dr}]")
            else:
                print(f" -> {shift_name:17} | \u0394\u03c6: {sys_twist_phi:+.4f} ({sys_sig_phi:4.1f}\u03c3) [{sys_status_phi:7}]")

        # MODULE E: LUMINOSITY SCALING EVOLUTION (PROJECTED SIGNIFICANCE)
        print("\n[MODULE E] LUMINOSITY SCALING EVOLUTION (SIGMA TRACKING)")
        FRACTIONS = [0.20, 0.50, 1.00]
        total_evts = len(mass)

        np.random.seed(42)
        indices = np.arange(total_evts)
        np.random.shuffle(indices)

        m_shuf = mass[indices]
        phi1_shuf = phi1[indices]
        phi2_shuf = phi2[indices]
        if has_eta:
            eta1_shuf = eta1[indices]
            eta2_shuf = eta2[indices]

        for frac in FRACTIONS:
            limit = int(total_evts * frac)
            m_sub = m_shuf[:limit]
            p1_sub = phi1_shuf[:limit]
            p2_sub = phi2_shuf[:limit]

            peak_mask_sub = (m_sub >= config["peak_min"]) & (m_sub <= config["peak_max"])
            bg_mask_sub = (m_sub >= 210.0) & (m_sub <= 220.0) & (~peak_mask_sub)

            if np.sum(peak_mask_sub) > 10 and np.sum(bg_mask_sub) > 10:
                dphi_sub = calc_delta_phi(p1_sub, p2_sub)
                _, sig_phi_sub = calc_significance(dphi_sub[peak_mask_sub], dphi_sub[bg_mask_sub])

                dr_str = ""
                if has_eta:
                    e1_sub = eta1_shuf[:limit]
                    e2_sub = eta2_shuf[:limit]
                    dr_sub = calc_delta_r(e1_sub, p1_sub, e2_sub, p2_sub)
                    _, sig_dr_sub = calc_significance(dr_sub[peak_mask_sub], dr_sub[bg_mask_sub])
                    dr_str = f"| \u0394R: {sig_dr_sub:4.2f} \u03c3"

                print(f" -> [LUMI {int(frac*100):3}%] \u0394\u03c6: {sig_phi_sub:4.2f} \u03c3 {dr_str}")
            else:
                print(f" -> [LUMI {int(frac*100):3}%] Insufficient data")

    print("\n==================================================")
    print(" [SUCCESS] UNIFIED AUDITOR EXECUTION COMPLETED. ")
    print("==================================================")
