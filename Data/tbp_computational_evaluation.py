# ==============================================================================
# THE CONNECTING THREAD THEORY (TBP): AB INITIO SYMBOLIC EVALUATION ENGINE
# PROTOCOL: EXACT RATIONAL & HIGH-PRECISION SYMPY VERIFICATION (PAPER I)
# SOURCE EQUATIONS: SSRN 6953858 / ZENODO 10.5281/zenodo.19978186
# REPOSITORY: https://github.com/cefiyana-clover/TBP-Computational-Framework
# ==============================================================================
import sympy as sp

print("==========================================================================")
print(" [INIT] TBP AB INITIO SYMBOLIC & FRACTIONAL COMPUTATIONAL EVALUATION      ")
print("        PURE FORMULA-DRIVEN DERIVATION FROM (d = 6, d_bulk = 4)           ")
print("==========================================================================\n")

# ------------------------------------------------------------------------------
# 0. FOUNDATIONAL AXIOMATIC INTEGERS & EXACT SYMPY OPERATORS (SECTIONS III-V)
# ------------------------------------------------------------------------------
d = sp.Integer(6)                      # Internal 6D manifold capacity (SU(3)_c -> R^6)
d_bulk = sp.Integer(4)                 # Macroscopic 4D spacetime bulk
d_int = d + d_bulk                     # Total 10D interaction phase space
chi = d_int ** 2                       # Phase-space tensor area (100)
vol_scale = d_int ** 3                 # Volumetric condensation ratio (1000)

# Eq. (3): Chiral Parity Bifurcation & Active Loop Constant
mu = (d ** 2) + (d_bulk ** 2) + (d ** 2)

# Eq. (1) & (5): Dimensionless Critical Eigenvalue & Vacuum Ceiling (Lambda_EW = 1 GeV)
lambda_crit = sp.Rational(mu ** 2, d ** 2)
R_crit_GeV = lambda_crit               # GeV (anchored at Lambda_EW = 1.000 GeV)
R_crit_MeV = R_crit_GeV * vol_scale    # MeV

# Eq. (2), (11), (14), (16): 1-Loop Friction, 2-Loop Backscattering & Damping
H_lim = sp.Integer(1) / ((d ** 2) * sp.sqrt(lambda_crit))
F1 = H_lim
F2 = (F1 ** 2) / (d ** 2)
Gamma_TBP = F1 - F2
contraction_ratio = H_lim / (d ** 2)

# Symmetry Node Fixation & Volumetric Dilution
active_nodes = (d ** 2) / 2 - 1        # 18 - 1 = 17
sym_ratio = sp.Rational(active_nodes, (d ** 2) // 2)   # 17 / 18
V_3D = sp.Rational(d ** 2, vol_scale)                  # 36 / 1000

print("[SECTION III-V] AXIOMATIC MATRIX INVARIANTS & DAMPED WAVE OPERATORS")
print(f" -> Internal Dimension (d)               : {d} (d^2 = {d**2})")
print(f" -> Bulk Dimension (d_bulk)              : {d_bulk} (d_bulk^2 = {d_bulk**2})")
print(f" -> Active Loop Constant (mu, Eq. 3)     : {mu}")
print(f" -> Critical Ceiling (R_crit, Eq. 1,5)   : {R_crit_GeV} = {float(R_crit_GeV.evalf(12)):.9f} GeV")
print(f" -> 1-Loop Friction (H_lim, Eq. 2)       : {H_lim} = {float(H_lim.evalf(12)):.10f}")
print(f" -> 2-Loop Backscattering (F^(2), Eq.75) : {F2} = {float(F2.evalf(12)):.10e}")
print(f" -> Net Damping (Gamma_TBP, Eq. 14)      : {Gamma_TBP} = {float(Gamma_TBP.evalf(12)):.10f}")
print(f" -> Perturbative Contraction (H_lim/d^2) : {contraction_ratio} = {float(contraction_ratio.evalf(10)):.8e}")
print(f" -> Symmetry-Breaking Ratio              : {sym_ratio} = {float(sym_ratio.evalf(10)):.8f}")
print(f" -> Volumetric Dilution (V_3D, Eq. 60)   : {V_3D} = {float(V_3D.evalf(10)):.6f}\n")

# ------------------------------------------------------------------------------
# 1. ELECTROMAGNETIC, STRONG & ELECTROWEAK SECTORS (SECTION VI)
# ------------------------------------------------------------------------------
print("==========================================================================")
print(" [SECTION VI.A] ELECTROMAGNETIC, STRONG & ELECTROWEAK SECTORS             ")
print("==========================================================================\n")

# Electromagnetic Sector (Eqs. 18-19)
Lambda_scalar = mu + (d ** 2) + 1
E_edges = 2 * d
alpha_ideal_inv = Lambda_scalar + E_edges + V_3D
alpha_bare_inv = alpha_ideal_inv + H_lim
alpha_obs_inv_codata = sp.Rational(137035999084, 1000000000)  # CODATA 2022 reference for L1 & C1

# Strong Interaction Sector (Eqs. 20-21)
container_cap = (d ** 2) * d_bulk
alpha_s = sp.Rational(active_nodes, container_cap)
L_QCD = H_lim / alpha_s

# Electroweak Sector (Eqs. 22-33)
sin2_theta_W = sp.Rational(2, 9)
cos2_theta_W = 1 - sin2_theta_W

# Unitary Node Energy E_node (Eqs. 23-24)
V4_hyp = (sp.pi ** 2) / 2
V6_hyp = (sp.pi ** 3) / 6
E_node_GeV = (V4_hyp / V6_hyp) + (sp.Rational(d, d_bulk) * H_lim)

# Z Boson (Eqs. 26-28)
M_Z_bare = mu + sp.pi
B_Z = sp.Rational(d_bulk, mu)
M_Z = M_Z_bare + B_Z

# W Boson (Eqs. 29-31 & Eq. 73)
M_W_bare = M_Z * sp.sqrt(cos2_theta_W)
dM_W = sp.Rational(mu, 4) * H_lim
M_W_1loop = M_W_bare - dM_W
F_W_2 = (dM_W ** 2) / (d ** 2)
M_W_2loop = M_W_1loop + F_W_2

# Higgs Boson (Eqs. 32-33)
M_H_bare = mu + (d ** 2) + 1
M_H_obs = M_H_bare + sin2_theta_W

print(f" -> alpha_ideal^-1 (Eq. 18)        : {alpha_ideal_inv} = {float(alpha_ideal_inv.evalf(12)):.9f}")
print(f" -> alpha_Bare^-1  (Eq. 19)        : {alpha_bare_inv} = {float(alpha_bare_inv.evalf(12)):.9f}")
print(f" -> alpha_s        (Eq. 20)        : {alpha_s} = {float(alpha_s.evalf(12)):.9f}")
print(f" -> L_QCD          (Eq. 21)        : {L_QCD} = {float(L_QCD.evalf(12)):.9f}")
print(f" -> sin^2(theta_W) (Eq. 22)        : {sin2_theta_W} = {float(sin2_theta_W.evalf(12)):.9f}")
print(f" -> E_node         (Eq. 23-24)     : {float(E_node_GeV.evalf(12)):.7f} GeV ({float((E_node_GeV*1000).evalf(12)):.4f} MeV)")
print(f" -> M_Z(Bare)      (Eq. 26)        : {float(M_Z_bare.evalf(12)):.6f} GeV")
print(f" -> B_Z            (Eq. 27)        : {B_Z} = {float(B_Z.evalf(12)):.6f} GeV")
print(f" -> M_Z(Obs)       (Eq. 28)        : {float(M_Z.evalf(12)):.6f} GeV")
print(f" -> M_W(Bare)      (Eq. 29)        : {float(M_W_bare.evalf(12)):.6f} GeV")
print(f" -> Delta M_W      (Eq. 30)        : {dM_W} = {float(dM_W.evalf(12)):.6f} GeV")
print(f" -> M_W(1-Loop)    (Eq. 31)        : {float(M_W_1loop.evalf(12)):.6f} GeV")
print(f" -> F_W^(2)        (Eq. 73)        : {F_W_2} = {float(F_W_2.evalf(12)):.6e} GeV")
print(f" -> M_W(2-Loop)    (Eq. 31+73)     : {float(M_W_2loop.evalf(12)):.6f} GeV")
print(f" -> M_H(Bare)      (Eq. 32)        : {M_H_bare} GeV")
print(f" -> M_H(Obs)       (Eq. 33)        : {M_H_obs} = {float(M_H_obs.evalf(12)):.6f} GeV\n")

# ------------------------------------------------------------------------------
# 2. FERMION MASS HIERARCHY, NEUTRINOS, QUARKS & PROTON (SECTION VI.1-4)
# ------------------------------------------------------------------------------
print("==========================================================================")
print(" [SECTION VI.B] FERMION MASS HIERARCHY, NEUTRINOS, QUARKS & PROTON        ")
print("==========================================================================\n")

# Metric Scale Operator S (Eq. 34, Lambda_cond = 1.000 MeV)
s_dimless = lambda_crit / chi
S_MeV = s_dimless

# 1. Lepton Sector (Eqs. 35-39)
L_1 = H_lim / alpha_obs_inv_codata
Omega_1 = 2 * sp.pi
Omega_2 = 2 * (d ** 2)
Omega_3 = (d_bulk * d) ** 2

delta_1 = Omega_1 * H_lim - 1 * L_1
delta_2 = Omega_2 * H_lim - 2 * L_1
delta_3 = Omega_3 * H_lim - 3 * L_1

LPhi_1 = 1 * sp.E * S_MeV
LPhi_2 = 2 * sp.E * S_MeV
LPhi_3 = 3 * sp.E * S_MeV

C_1 = (2 * d) / alpha_obs_inv_codata
C_2 = d + sp.pi
C_3 = 3 * (d ** 2)

M_e = (LPhi_1 - delta_1) * C_1
M_mu = (LPhi_2 - delta_2) * C_2
M_tau = (LPhi_3 - delta_3) * C_3

print("[1. LEPTON SECTOR (MeV)]")
print(f" -> Metric Scale Operator S (Eq. 34): {S_MeV} = {float(S_MeV.evalf(12)):.6f} MeV")
print(f" -> Vacuum Turbulence L_1   (Eq. 35): {float(L_1.evalf(12)):.8e}")
print(f" -> delta_Final(1,2,3)      (Eq. 35): {float(delta_1.evalf(12)):.8f}, {float(delta_2.evalf(12)):.8f}, {float(delta_3.evalf(12)):.8f}")
print(f" -> Electron Mass M_e       (Eq. 37): {float(M_e.evalf(12)):.6f} MeV (2*M_e = {float((2*M_e).evalf(12)):.6f} MeV)")
print(f" -> Muon Mass     M_mu      (Eq. 38): {float(M_mu.evalf(12)):.6f} MeV")
print(f" -> Tau Mass      M_tau     (Eq. 39): {float(M_tau.evalf(12)):.6f} MeV\n")

# 2. Neutrino Dissipation Sector (Eqs. 40-44)
R_crit_sq = lambda_crit ** 2
M_nu_e_eV   = (delta_1 / R_crit_sq) * sp.Integer(10 ** 6)
M_nu_mu_eV  = (delta_2 / R_crit_sq) * sp.Integer(10 ** 6)
M_nu_tau_eV = (delta_3 / R_crit_sq) * sp.Integer(10 ** 6)

print("[2. NEUTRINO DISSIPATION SECTOR (eV)]")
print(f" -> Universal Divisor R_crit^2 (Eq. 40): {R_crit_sq} = {float(R_crit_sq.evalf(12)):.6f} GeV^2")
print(f" -> Electron Neutrino M_nu_e   (Eq. 42): {float(M_nu_e_eV.evalf(12)):.6f} eV")
print(f" -> Muon Neutrino     M_nu_mu  (Eq. 43): {float(M_nu_mu_eV.evalf(12)):.6f} eV")
print(f" -> Tau Neutrino      M_nu_tau (Eq. 44): {float(M_nu_tau_eV.evalf(12)):.6f} eV\n")

# 3. Quark Sector (Eqs. 45-52)
Delta_QCD_MeV = L_QCD * S_MeV
M_u_MeV = S_MeV + 1 * Delta_QCD_MeV
M_d_MeV = (S_MeV ** 2) + 2 * Delta_QCD_MeV
M_s_MeV = (mu + d) + 3 * Delta_QCD_MeV
M_c_MeV = (R_crit_MeV / ((2 * d + 1) ** 2)) + 4 * Delta_QCD_MeV
M_b_MeV = (sp.Rational(chi, 24) * vol_scale) + 5 * Delta_QCD_MeV
M_t_bare_GeV = lambda_crit - ((d ** 2) + d)
M_t_GeV = M_t_bare_GeV + (6 * Delta_QCD_MeV / vol_scale)

print("[3. QUARK SECTOR (Geometric Regime)]")
print(f" -> Gluon Cloud Layer Delta_QCD (Eq. 45): {Delta_QCD_MeV} = {float(Delta_QCD_MeV.evalf(12)):.6f} MeV")
print(f" -> Up Quark      M_u (Eq. 46)          : {float(M_u_MeV.evalf(12)):.6f} MeV")
print(f" -> Down Quark    M_d (Eq. 47)          : {float(M_d_MeV.evalf(12)):.6f} MeV")
print(f" -> Strange Quark M_s (Eq. 48)          : {float(M_s_MeV.evalf(12)):.6f} MeV")
print(f" -> Charm Quark   M_c (Eq. 49)          : {float(M_c_MeV.evalf(12)):.6f} MeV")
print(f" -> Bottom Quark  M_b (Eq. 50)          : {float(M_b_MeV.evalf(12)):.6f} MeV")
print(f" -> Top Bare      M_t(Bare) (Eq. 51)    : {M_t_bare_GeV} = {float(M_t_bare_GeV.evalf(12)):.6f} GeV")
print(f" -> Top Physical  M_t (Eq. 52)          : {float(M_t_GeV.evalf(12)):.6f} GeV\n")

# 4. Proton Equilibrium (Eqs. 53-56 & Eq. 74)
E_Core_MeV = R_crit_MeV * L_QCD
M_p_bare = E_Core_MeV * (sp.E / sp.sqrt(chi))
dM_p = chi * H_lim
M_p_1loop = M_p_bare + dM_p
F_p_2 = (dM_p ** 2) / (d ** 2)
M_p_2loop = M_p_1loop - F_p_2

print("[4. PROTON EQUILIBRIUM (MeV)]")
print(f" -> Hadronic Core E_Core (Eq. 53)       : {E_Core_MeV} = {float(E_Core_MeV.evalf(12)):.6f} MeV")
print(f" -> Proton Bare M_p(Bare)(Eq. 54)       : {float(M_p_bare.evalf(12)):.6f} MeV")
print(f" -> Friction Shift dM_p  (Eq. 55)       : {dM_p} = {float(dM_p.evalf(12)):.6f} MeV")
print(f" -> Proton 1-Loop M_p    (Eq. 56)       : {float(M_p_1loop.evalf(12)):.6f} MeV")
print(f" -> 2-Loop Shift F_p^(2) (Eq. 74)       : {F_p_2} = {float(F_p_2.evalf(12)):.6e} MeV")
print(f" -> Proton 2-Loop M_p    (Eq. 56-74)    : {float(M_p_2loop.evalf(12)):.6f} MeV\n")

# ------------------------------------------------------------------------------
# 3. DETERMINISTIC COSMOLOGY & GRAVITATIONAL DILUTION (SECTION VII)
# ------------------------------------------------------------------------------
print("==========================================================================")
print(" [SECTION VII] DETERMINISTIC COSMOLOGY & GRAVITATIONAL DILUTION           ")
print("==========================================================================\n")

# Newton's Constant & Planck Scale (Eqs. 57-62) using exact CODATA 2022 SI constants
hbar_SI = sp.Float("1.054571817e-34", 20)
c_SI    = sp.Float("2.99792458e8", 20)
GeV_to_kg = sp.Float("1.7826619216279e-27", 20)

M_P_GeV = sp.exp(mu / 2)
M_P_kg = M_P_GeV * GeV_to_kg
G_Base = (hbar_SI * c_SI) / (M_P_kg ** 2)
V_G = V_3D + (1 + H_lim) / chi
G_Final = G_Base * (1 / sym_ratio) * (1 + V_G)

print("[1. NEWTON'S CONSTANT & PLANCK SCALE]")
print(f" -> Geometric Planck Mass M_P (Eq. 57)  : e^{mu//2} = {float(M_P_GeV.evalf(12)):.6e} GeV")
print(f" -> Planck Mass in SI         (Eq. 58)  : {float(M_P_kg.evalf(12)):.6e} kg")
print(f" -> Baseline Gravity G_Base   (Eq. 59)  : {float(G_Base.evalf(12)):.6e} m^3 kg^-1 s^-2")
print(f" -> Regge-Volumetric Factor V_G (Eq.61) : {V_G} = {float(V_G.evalf(12)):.8f}")
print(f" -> Macroscopic G_Final       (Eq. 62)  : {float(G_Final.evalf(12)):.6e} m^3 kg^-1 s^-2\n")

# Cosmic Budget (Eqs. 63-68)
OM_phys = sp.Rational(d_bulk, mu) * sym_ratio
OM_obs  = OM_phys * (1 + V_3D)

DM_phys = sp.Rational(d * d_bulk, mu) * sym_ratio
DM_obs  = DM_phys * (1 + V_3D)

DE_phys = sp.Rational(mu - d_bulk - (d * d_bulk), mu) + sp.Rational(d_bulk + (d * d_bulk), mu) * (1 - sym_ratio)
Omega_r = H_lim
DE_obs  = 1 - (OM_obs + DM_obs + Omega_r)
w_DE    = -1 + (H_lim / 3)

print("[2. DETERMINISTIC COSMIC BUDGET]")
print(f" -> Ordinary Matter (OM_obs, Eq. 64)    : {OM_obs} = {float((OM_obs*100).evalf(12)):.5f}%")
print(f" -> Dark Matter     (DM_obs, Eq. 66)    : {DM_obs} = {float((DM_obs*100).evalf(12)):.5f}%")
print(f" -> Relic Radiation (Omega_r = H_lim)   : {Omega_r} = {float((Omega_r*100).evalf(12)):.5f}%")
print(f" -> Dark Energy     (DE_phys, Eq. 67)   : {DE_phys} = {float((DE_phys*100).evalf(12)):.5f}%")
print(f" -> Dark Energy     (DE_obs,  Eq. 68)   : {DE_obs} = {float((DE_obs*100).evalf(12)):.5f}%")
print(f" -> Flatness Check  (Sum of Omegas)     : {OM_obs + DM_obs + Omega_r + DE_obs}")
print(f" -> Equation of State w_DE (-1+H_lim/3) : {w_DE} = {float(w_DE.evalf(12)):.6f}\n")

# ------------------------------------------------------------------------------
# 4. SENSITIVITY ANALYSIS, PERTURBATIVE STABILITY & OBSERVER BIAS (SECTION VIII)
# ------------------------------------------------------------------------------
print("==========================================================================")
print(" [SECTION VIII] SENSITIVITY ANALYSIS, 2-LOOP STABILITY & OBSERVER BIAS    ")
print("==========================================================================\n")

# Dimensional Stability Margin (Eqs. 69-70)
mu_5 = 2 * (5 ** 2) + 16
R_crit_5 = sp.Rational(mu_5 ** 2, 5 ** 2)
mu_7 = 2 * (7 ** 2) + 16
M_H_bare_7 = mu_7 + (7 ** 2) + 1

# 2-Loop Lepton Vacuum Turbulence Rescattering (Eq. 72)
F_e_2 = (L_1 ** 2) / (d ** 2)

# Local Observer Equation & Terrestrial Vacuum Microlensing (Eqs. 75-81)
Delta_alpha_emp = alpha_ideal_inv - alpha_obs_inv_codata
Phi_MW = (sp.Float("2.2e5") / sp.Float("2.99e8")) ** 2
Phi_Sun = sp.Float("9.87e-9")
Phi_Earth = sp.Float("6.96e-10")
Sigma_Phi_orb = Phi_MW + Phi_Sun + Phi_Earth
Phi_MW_full = Phi_MW * (1 + sp.Rational(24, 36)) * (1 + V_3D)
Sigma_Phi_full = Phi_MW_full + Phi_Sun + Phi_Earth

print(f" -> d = 5 Sensitivity (Eq. 69)          : mu_5 = {mu_5}, R_crit(5) = {R_crit_5} = {float(R_crit_5.evalf(10)):.2f} GeV")
print(f" -> d = 7 Sensitivity (Eq. 70)          : mu_7 = {mu_7}, M_H(Bare,7) = {M_H_bare_7} GeV")
print(f" -> 2-Loop Lepton F_e^(2) (Eq. 72)      : {float(F_e_2.evalf(12)):.4e} MeV")
print(f" -> Empirical CODATA Gap (Eq. 76)       : {float(Delta_alpha_emp.evalf(12)):.3e}")
print(f" -> Hypothetical 2-Loop F_a^(2) (Eq.75) : {float(F2.evalf(12)):.3e}")
print(f" -> Orbital Potential Lower Bound (Eq81): {float(Sigma_Phi_orb.evalf(12)):.3e}")
print(f" -> Full Halo Potential Upper Bound     : {float(Sigma_Phi_full.evalf(12)):.3e}")

print("\n==========================================================================")
print(" [SUCCESS] TBP AB INITIO SYMBOLIC EVALUATION COMPLETED.                   ")
print("==========================================================================")
