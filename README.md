# The Connecting Thread Theory (TBP): Computational & Multi-Messenger Phenomenological Auditing Framework

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg)](https://www.python.org/downloads/)
[![SymPy](https://img.shields.io/badge/SymPy-Symbolic_Verification-3B5526.svg)](https://www.sympy.org/)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0008--4324--9515-A6CE39.svg)](https://orcid.org/0009-0008-4324-9515)
[![DOI: Zenodo](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.19978186-blue.svg)](https://doi.org/10.5281/zenodo.19978186)
[![SSRN: Theory](https://img.shields.io/badge/SSRN-6953858-darkred.svg)](https://doi.org/10.2139/ssrn.6953858)
[![SSRN: Collider](https://img.shields.io/badge/SSRN-7075779-darkred.svg)](https://doi.org/10.2139/ssrn.7075779)

---

## 1. Overview & Axiomatic Architecture

This repository provides the official open-source symbolic computation scripts and empirical data-auditing pipelines for **The Connecting Thread Theory (TBP)**.

TBP formulates a discrete 6-dimensional topological matrix framework designed to derive fundamental particle masses, gauge couplings, leptonic mixing angles, and macroscopic gravitational/cosmological metrics *ab initio* without adjustable numerical fitting parameters. The mathematical architecture is anchored on six foundational postulates derived from three structural integers:

1. **Internal 6D Manifold Capacity ($d = 6$):** Derived from the real linear-algebraic isomorphism ($\mathbb{C}^3 \cong \mathbb{R}^6$) required to house the 3-dimensional complex color space of $SU(3)_c$ in Quantum Chromodynamics, defining a rank-2 spatial matrix capacity of $N_{\text{matrix}}(6) = d^2 = 36$ components ($21$ symmetric metric + $15$ antisymmetric spin-torsion components on a discrete causal graph $\mathcal{G}_6$).
2. **Macroscopic 4D Spacetime Bulk ($d_{bulk} = 4$):** Operational trajectories project onto a $3+1$ dimensional observable spacetime governed by a rank-2 operator matrix ($M_{bulk}$) containing $N_{\text{matrix}}(4) = d_{bulk}^2 = 16$ components ($10$ symmetric metric components $g_{(\mu\nu)}$ and $6$ antisymmetric spin-torsion components $g_{[\mu\nu]}$; canonically decomposing into $4$ diagonal trace nodes $\mu=\nu$ and $12$ off-diagonal shear nodes $\mu\neq\nu$).
3. **Chiral Parity Bifurcation & Active Loop Constant ($\mu = 88$):** Under the chiral projection operators $P_{L,R} = (1 \mp \gamma^5)/2$, weak parity violation bifurcates the $d=6$ manifold into orthogonal Left-Chiral ($M_{6D_L} = 36$ nodes) and Right-Chiral ($M_{6D_R} = 36$ nodes) boundaries sandwiching the central 4D bulk ($M_{bulk} = 16$ nodes), fixing the total topological loop load at:
   $$\mu = \dim(M_{6D_L}) + \dim(M_{bulk}) + \dim(M_{6D_R}) = 6^2 + 4^2 + 6^2 = 36 + 16 + 36 = 88$$
4. **Critical Vacuum Ceiling ($R_{crit}$) & Heaviside Truncation ($\Theta$):** Bounded by the ratio of the squared active loop constant to the internal spatial capacity anchored at the unitary electroweak metric scale $\Lambda_{EW} \equiv 1.000\text{ GeV}$:
   $$R_{crit} = \left(\frac{\mu^2}{d^2}\right)\Lambda_{EW} = \frac{88^2}{36}\text{ GeV} = \frac{1936}{9}\text{ GeV} \approx 215.111111\text{ GeV}$$
   The Grand Action $S_{TBP}$ incorporates the Lorentz-invariant Heaviside cutoff $\Theta(1 - \mathcal{M}_i^2/R_{crit}^2)$. As single-node invariant mass-energy $\mathcal{M}_i \to R_{crit}$, 4D action probabilities truncate ($\Theta \to 0$) and excess momentum is transferred orthogonally into the 2D internal boundary via the vielbein transition tensor $E_b^\nu$ (**Dimensional Leakage**) under a globally self-adjoint $2\times 2$ block Hamiltonian $\mathcal{H}_{6D} = \mathcal{H}_{6D}^\dagger$ ($S_{6D}^\dagger S_{6D} = \mathbb{I}$).
5. **1-Loop Geometric Friction Limit ($H_{lim}$) & Perturbative Backscattering ($F^{(n)}$):** Wave propagation across the discrete lattice encounters a deterministic first-order impedance limit and second-order backscattering:
   $$H_{lim} \equiv F^{(1)} = \frac{1}{d^2\sqrt{R_{crit}/\Lambda_{EW}}} = \frac{1}{36 \times (44/3)} = \frac{1}{528} \approx 0.00189394, \quad F^{(2)} = \frac{(F^{(1)})^2}{d^2} = \frac{1}{10036224} \approx 9.9639 \times 10^{-8}$$
6. **Symmetry Node Fixation ($17/18$) & Volumetric Phase-Space Dilution ($V_{3D}$):** Binary $SU(2)$ polarization divides the 36 internal nodes into 18 matter and 18 vacuum nodes; fixing 1 static reference node to prevent unconstrained rotational instability leaves $18 - 1 = 17$ active matter nodes (asymmetry ratio $17/18$). Mapping $d=6$ onto $d_{bulk}=4$ defines a 10D interaction phase space ($d_{int} = 10$), tensor area $\chi = d_{int}^2 = 100$, volumetric condensation ratio $\Lambda_{EW}/\Lambda_{cond} = d_{int}^3 = 1000$ ($\Lambda_{cond} = 1.000\text{ MeV}$), and 3D volumetric dilution factor $V_{3D} = d^2/\chi^{3/2} = 36/1000 = 0.036$.

---

## 2. Corpus Architecture & Pipeline Mapping (8-Paper Series)

This repository is organized around the complete 8-manuscript TBP research corpus (1 Foundational Theoretical Manuscript + 7 Companion Phenomenological Audits):

| Paper | Domain / Scope | Repository Script(s) | Primary Datasets Audited | Key Empirical Verdict |
| :--- | :--- | :--- | :--- | :--- |
| **Paper I** | **Foundational 6D Theory:** Unified Matrix Architecture & Symbolic Derivations | `tbp_computational_evaluation.py` | PDG (2024), CODATA (2022), Planck (2018), KATRIN (2024) | 24 fundamental constants derived with zero free fitting parameters; Lyapunov & $S_{6D}$ unitarity proven. |
| **Paper II** | **Collider Ceiling ($215.111\text{ GeV}$):** Topological Deformation & Kinematic Divergence | `tbp_cms_omni_extractor.py`<br>`tbp_cms_215_evaluator.py`<br>`tbp_cms_215_auditor.py` | CERN CMS 2016 Legacy $\sqrt{s}=13\text{ TeV}$ Open Data ($>1.5\times 10^8$ events, 10 channels) | $+8.43\sigma$ JetHT mass bottleneck; $+4.77\sigma$ electron reflection; $+6.05\text{ GeV}$ ($5.44\sigma$) $e^\pm$ $E_T^{miss}$ shift; 4D Euclidean rigidity ($\Delta\phi, \Delta R < 5\sigma$). |
| **Paper III** | **Top-Quark Threshold ($173.111\text{ GeV}$):** Sub-Ceiling Control & HLT Diagnostics | `tbp_cms_173_top_auditor.py` | CERN CMS 2016 Legacy $\sqrt{s}=13\text{ TeV}$ (93 Four-Vector `.npy` Matrices, $130\text{--}210\text{ GeV}$) | Hadronic pileup vanishes ($T=-0.66$); prompt $t\to Wb\to \ell\nu_\ell b$ neutrino surge in $\mu$ ($T=+4.98$) & $e$ ($T=+3.53$); HLT plateau verified. |
| **Paper IV** | **Gravitational Waves:** Topological Low-Pass Filtering & Jitter Suppression | `tbp_ligo_omni_auditor_v2.py` | LIGO-Virgo GWTC-1 Catalog (10 BBH Mergers, 20 H1/L1 Streams at $4096\text{ Hz}$) | Zero quantum jitter $>528\text{ Hz}$ ($\langle Z_{HF}\rangle = 1.24\sigma$, max $3.56\sigma$); High-SNR ($\ge 5.50$) Parseval residual $\mathcal{R}_{att} = 0.78\times 10^{-3} < H_{lim}$. |
| **Paper V** | **Precision Metrology & Quasars:** Vacuum Microlensing at $\alpha_{ideal}^{-1} = 137.036000$ | `tbp_metrology_quasar_auditor_v2.py` | CODATA 2014--2022, Cs-133/Rb-87 Recoil, VLT/Keck Quasars ($N=293$), Fermilab $(g-2)$ | 2-loop $F_\alpha^{(2)}$ rejected at $>24\sigma$; lab gap ($8.79\times 10^{-7}$) inside Milky Way bracket $[5.48, 9.45]\times 10^{-7}$; $\Delta a_\mu^{\text{TBP}} = 2.27\times 10^{-9}$ ($0.46\sigma$). |
| **Paper VI** | **Galactic Dynamics (SPARC):** $24/88$ Cross-Manifold Torque & Scale $a_{TBP}$ | `tbp_sparc_torque_auditor_v2.py` | SPARC Database ($N=280$ Radial Points across 10 Archetypal Galaxies) | $a_{TBP} = cH_0(289/1584) = 1.1947\times 10^{-10}\text{ m/s}^2$ ($99.56\%, 0.26\sigma$); RAR RMS $0.0183\text{ dex}$ ($r=0.99900$); halo ratio $5.76\pm 0.83$ ($0.29\sigma$ from $6.00$). |
| **Paper VII** | **Neutrino Sector:** $R_{crit}^2$ Mass Dilution, Chiral PMNS Mixing & IceCube HESE | `tbp_neutrino_icecube_auditor_v2.py` | KATRIN KNM1--5, PDG 2024 / NuFIT 5.3 PMNS, IceCube 7.5-Yr HESE Catalog | $M_{\nu_e} = 0.25687\text{ eV}$ ($<0.45\text{ eV}$ KATRIN); PMNS angles within $0.03\sigma\text{--}0.34\sigma$; IceCube index $\gamma_{TBP} = 2.8691$ ($<0.01\sigma, \chi_\nu^2 = 0.16$). |
| **Paper VIII** | **Black Hole Horizons:** Pre-Leakage Boundary States Across 4 Messengers | `tbp_blackhole_horizon_auditor_v2.py` | Fermi-LAT/Cherenkov Blazars, XRISM/XMM/NuSTAR $\text{Fe K}\alpha$, IceCube AGN, LIGO | Intrinsic blazar breaks cluster at $215.17\pm 5.41\text{ GeV}$ ($0.01\sigma$); 2PN disk microlensing $\chi_\nu^2 = 0.30$; coronal $\Phi_\nu/\Phi_\gamma = 1/L_{QCD} = 62.33\times$. |

---

## 3. Consolidated Summary of Ab Initio TBP Derivations (`tbp_computational_evaluation.py`)

The core symbolic evaluation engine verifies the dimensional homogeneity and exact fractional values of 24 fundamental parameters across the Electroweak, Strong, Fermionic, and Cosmological sectors:

| Fundamental Parameter | TBP Topological Formula | TBP 1-Loop Ideal | TBP Corrected (2-Loop / Bias) | Empirical Benchmark (PDG / CODATA / Planck) |
| :--- | :--- | :--- | :--- | :--- |
| **Critical Vacuum Ceiling ($R_{crit}$)** | $(\mu^2 / d^2)\Lambda_{EW}$ | $1936/9 \approx 215.111111\text{ GeV}$ | $215.111111\text{ GeV}$ | Tested at CMS / Blazar Jets ($215.17 \pm 5.41\text{ GeV}$) |
| **1-Loop Geometric Friction ($H_{lim}$)** | $1 / (d^2\sqrt{R_{crit}/\Lambda_{EW}})$ | $1/528 \approx 0.00189394$ | $0.00189384$ ($\Gamma_{TBP}$) | Tested via LIGO GWTC-1 ($\mathcal{R}_{att} < 1.894\times 10^{-3}$) |
| **Unitary Node Energy ($E_{node}$)** | $(V_4/V_6 + \frac{d}{d_{bulk}}H_{lim})\Lambda_{EW}$ | $(3/\pi + 1/352)\text{ GeV}$ | $957.7706\text{ MeV}$ | $957.78 \pm 0.06\text{ MeV}$ ($\eta'$ axial-anomaly meson) |
| **Inverse Fine-Structure ($\alpha^{-1}$)** | $(\mu + d^2 + 1) + 2d + V_{3D}$ | $137.036000$ | $137.0359991$ ($-\Sigma\Phi_{local}$) | $137.035999084 \pm 0.000000021$ (CODATA 2022) |
| **Strong Coupling Constant ($\alpha_s$)** | $17 / (d^2 \cdot d_{bulk})$ | $17/144 \approx 0.118056$ | $0.118056$ (Static) | $0.1179 \pm 0.0009$ (PDG 2024 at $M_Z$) |
| **1-Loop QCD Feedback ($L_{QCD}$)** | $H_{lim} / \alpha_s$ | $3/187 \approx 0.0160428$ | $0.0160428$ | Governs $\Delta_{QCD}$ & AGN $\Phi_\nu/\Phi_\gamma = 187/3$ |
| **Electroweak Mixing ($\sin^2\theta_W$)** | $2 / 9$ | $0.222222\dots$ | $0.222222$ | $0.2223$ (On-shell $1 - M_W^2/M_Z^2$) |
| **$Z$ Boson Mass ($M_Z$)** | $(\mu + \pi + d_{bulk}/\mu)\Lambda_{EW}$ | $91.187047\text{ GeV}$ | $91.187047\text{ GeV}$ | $91.1876 \pm 0.0021\text{ GeV}$ |
| **$W$ Boson Mass ($M_W$)** | $M_Z\sqrt{7/9} - (\mu/4)H_{lim}\Lambda_{EW}$ | $80.377751\text{ GeV}$ | $80.377798\text{ GeV}$ (2-Loop) | $80.377 \pm 0.012\text{ GeV}$ (PDG 2024) |
| **Higgs Boson Mass ($M_H$)** | $(\mu + d^2 + 1 + \sin^2\theta_W)\Lambda_{EW}$ | $125.222222\text{ GeV}$ | $125.222222\text{ GeV}$ | $125.25 \pm 0.17\text{ GeV}$ |
| **Metric Scale Operator ($S$)** | $(\lambda_{crit}/\chi)\Lambda_{cond}$ | $1936/900\text{ MeV}$ | $2.151111\text{ MeV}$ | Fermionic baseline ($s \equiv 1936/900$ dimensionless) |
| **Electron Mass ($M_e$)** | $(e\cdot S - \delta_{Final(1)})(12/\alpha_{obs}^{-1})$ | $0.510999\text{ MeV}$ | $0.510999\text{ MeV}$ | $0.51099895\text{ MeV}$ |
| **Muon Mass ($M_\mu$)** | $(2e\cdot S - \delta_{Final(2)})(d + \pi)$ | $105.661421\text{ MeV}$ | $105.661421\text{ MeV}$ | $105.6583755\text{ MeV}$ |
| **Tau Mass ($M_\tau$)** | $(3e\cdot S - \delta_{Final(3)})(3d^2)$ | $1776.719999\text{ MeV}$ | $1776.719999\text{ MeV}$ | $1776.86 \pm 0.12\text{ MeV}$ |
| **Up Quark Mass ($M_u$)** | $S + 1\cdot\Delta_{QCD}$ | $2.185622\text{ MeV}$ | $2.185622\text{ MeV}$ | $2.16^{+0.49}_{-0.26}\text{ MeV}$ |
| **Down Quark Mass ($M_d$)** | $(S^2/\Lambda_{cond}) + 2\cdot\Delta_{QCD}$ | $4.696301\text{ MeV}$ | $4.696301\text{ MeV}$ | $4.67^{+0.48}_{-0.17}\text{ MeV}$ |
| **Strange Quark Mass ($M_s$)** | $(\mu + d)\Lambda_{cond} + 3\cdot\Delta_{QCD}$ | $94.103534\text{ MeV}$ | $94.103534\text{ MeV}$ | $93.4^{+8.6}_{-3.4}\text{ MeV}$ |
| **Charm Quark Mass ($M_c$)** | $R_{crit}/(2d+1)^2 + 4\cdot\Delta_{QCD}$ | $1272.984857\text{ MeV}$ | $1272.984857\text{ MeV}$ | $1270 \pm 20\text{ MeV}$ |
| **Bottom Quark Mass ($M_b$)** | $(\chi/24)\Lambda_{EW} + 5\cdot\Delta_{QCD}$ | $4166.839223\text{ MeV}$ | $4166.839223\text{ MeV}$ | $4180^{+30}_{-20}\text{ MeV}$ |
| **Top Quark Mass ($M_t$)** | $[\lambda_{crit} - (d^2+d)]\Lambda_{EW} + 6\Delta_{QCD}$ | $173.111318\text{ GeV}$ | $173.111318\text{ GeV}$ | $172.5\text{--}173.0\text{ GeV}$ |
| **Proton Mass ($M_p$)** | $R_{crit}L_{QCD}(e/10) + 100H_{lim}\Lambda_{cond}$ | $938.2631\text{ MeV}$ | $938.2621\text{ MeV}$ (2-Loop) | $938.272088\text{ MeV}$ |
| **Electron Neutrino ($M_{\nu_e}$)** | $(\delta_{Final(1)} / R_{crit}^2)\times 10^6\text{ eV}$ | $0.25687\text{ eV}$ | $0.25687\text{ eV}$ (Ground) | $< 0.45\text{ eV}$ (KATRIN 2024, $90\%$ CL) |
| **Newton's Constant ($G_{Final}$)** | $\frac{\hbar c}{(e^{\mu/2}\text{ GeV})^2}\left(\frac{18}{17}\right)(1 + V_G)$ | $6.67126 \times 10^{-11}$ | $6.67126 \times 10^{-11}$ | $6.67430 \times 10^{-11}\text{ m}^3\text{kg}^{-1}\text{s}^{-2}$ |
| **Cosmic Budget ($\Omega_{\Lambda}, \Omega_c, \Omega_b$)** | Partitions of $\mu = 60 + 24 + 4$ | $68.68\%, 26.68\%, 4.45\%$ | $w_{DE} = -0.99937$ | $68.31\pm 0.61\%, 26.81\pm 0.41\%, 4.91\%$ |

---

## 4. Detailed Architecture of the 7 Phenomenological Auditing Pipelines

### Pipeline 1 & 2: CERN CMS $13\text{ TeV}$ Collider Audits at $215.111\text{ GeV}$ & $173.111\text{ GeV}$
* **Scripts:** `tbp_cms_omni_extractor.py`, `tbp_cms_215_evaluator.py`, `tbp_cms_215_auditor.py`, `tbp_cms_173_top_auditor.py`
* **Data Provenance:** Streams reconstructed AOD/MiniAOD/NanoAOD collision records from the **CERN Open Data Portal** (CMS 2016 Legacy $\sqrt{s}=13\text{ TeV}$ primary datasets), sanitized via Golden JSON (`Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt`), six hardware noise filters (`Flag_goodVertices`, `Flag_globalSuperTightHalo2016Filter`, `Flag_HBHENoiseFilter`, `Flag_HBHENoiseIsoFilter`, `Flag_EcalDeadCellTriggerPrimitiveFilter`, `Flag_BadPFMuonFilter`), explicit High-Level Triggers (HLT), and Tight identification cuts across 10 independent topologies (93 four-vector `.npy` matrices).
* **Key Findings:**
  1. **At the $215.111\text{ GeV}$ Ceiling ($\Theta \to 0$):** Isolates a $+8.43\sigma$ hadronic invariant mass condensation bottleneck in `JetHT` ($429,625$ observed vs. $424,125.26$ expected) with suppressed transverse leakage ($-0.93\text{ GeV}$), paired with a $+4.77\sigma$ single-electron boundary reflection pileup and a monotonically scaling $+6.05\text{ GeV}$ ($T = +5.44, p = 6.32\times 10^{-8}$) Missing Transverse Energy ($E_T^{miss}$) shift, while preserving 4D Euclidean spatial rigidity ($\Delta\phi, \Delta R < 5.0\sigma$) across nominal and $\pm 5\%$ JES/LES perturbed matrices.
  2. **At the $173.111\text{ GeV}$ Top-Quark Threshold ($\Theta = 1$ Negative Control):** Within the isolated $130\text{--}210\text{ GeV}$ window, hadronic pileup and compression vanish (`JetHT` $\Delta E_T^{miss} = -0.05\text{ GeV}, T = -0.66$), while on-shell top-quark cascades ($t \to Wb \to \ell\nu_\ell b$) produce a clean, monotonically scaling prompt neutrino $E_T^{miss}$ surge in `Single-Muon` ($+4.35\text{ GeV}, T = +4.98$) and `Single-Electron` ($+2.47\text{ GeV}, T = +3.53$). Unmasked PLR diagnostics transparently isolate the `HLT_Photon165/175` prescale turn-on step ($+45.73\sigma$) below $200\text{ GeV}$, proving that the $200\text{--}230\text{ GeV}$ window used for the $215.111\text{ GeV}$ audit operates on the fully efficient trigger plateau.

### Pipeline 3: LIGO-Virgo GWTC-1 Gravitational-Wave Low-Pass Filter Audit
* **Script:** `tbp_ligo_omni_auditor_v2.py`
* **Data Provenance:** Streams all 10 binary black hole (BBH) coalescences (20 Hanford H1 and Livingston L1 strain streams sampled at $4096\text{ Hz}$, $T = 32\text{ s}$) from the **Gravitational Wave Open Science Center (GWOSC)**.
* **Methodology & Findings:** Applies median Welch PSD whitening and zero-phase fourth-order Butterworth bandpass filtering to partition each stream into **Band A (Macroscopic IMR: $35\text{--}350\text{ Hz}$)** and **Band B (High-Frequency Lattice Test Band: $528\text{--}1600\text{ Hz}$, above the second quasi-normal mode harmonic $f_c \ge 2f_{QNM} \approx 528\text{ Hz}$)**.
  1. **Universal Absence of Quantum Lattice Jitter:** Across all 10 mergers, on-source ($\pm 100\text{ ms}$) high-frequency strain energy exhibits zero deviations above the $5.0\sigma$ threshold ($\langle Z_{HF}\rangle = 1.24\sigma \pm 0.94\sigma$, maximum $3.56\sigma$ in GW170814).
  2. **Compliance with the $H_{lim} = 1/528$ Parseval Bound:** In High-SNR Sub-Ensemble A ($\text{SNR}_{Net} \ge 5.50$: GW150914, GW170814, GW170104, GW170818, GW170809), the bandwidth-normalized coherent high-frequency spectral residual ratio is strictly bounded below $H_{lim} = 1/528 \approx 1.894 \times 10^{-3}$, averaging $\langle\mathcal{R}_{att}\rangle = (0.775 \pm 0.326) \times 10^{-3}$. Apparent excursions in Low-SNR Sub-Ensemble B ($\text{SNR}_{Net} < 5.50$) are quantitatively resolved by low-mass inspiral chirp overshoot above $350\text{ Hz}$ ($M_{tot} \approx 19\text{--}21.7\,M_\odot$ in GW170608 and GW151226) and high redshift ($z \approx 0.48$ in GW170729).

### Pipeline 4: Terrestrial Vacuum Microlensing, Quasar Dipoles & Muon $(g-2)$ Audit
* **Script:** `tbp_metrology_quasar_auditor_v2.py`
* **Data Provenance:** Evaluates six terrestrial $\alpha_{obs}^{-1}$ determinations (CODATA 2014/2018/2022, Cs-133 recoil, Rb-87 recoil, electron $g_e - 2$), 293 Keck/HIRES and VLT/UVES quasar absorption systems ($0.22 \le z \le 4.18$, CDS VizieR `J/MNRAS/422/3370`), HST/STIS white dwarf G191-B2B spectra, and the Fermilab/BNL muon $(g-2)$ consensus.
* **Key Findings:**
  1. **Rejection of 2-Loop Internal Perturbation:** Both the independent primary experimental mean ($\langle\Delta\alpha_{emp}^{-1}\rangle_{\text{primary}} = (8.607 \pm 0.679) \times 10^{-7}$) and full benchmark mean ($(8.792 \pm 0.547) \times 10^{-7}$) reject a hypothetical 2-loop internal canvas perturbation ($F_\alpha^{(2)} = 9.964 \times 10^{-8}$) by $24.6\sigma$ to $63.1\sigma$ while residing strictly inside the parameter-free Milky Way virial potential bracket $\Sigma\Phi_{local}^{\text{TBP}} \in [5.52 \times 10^{-7}, 9.45 \times 10^{-7}]$.
  2. **Quasar Dipole & Muon $(g-2)$ Resolution:** Resolves the $8.14\sigma$ Keck/VLT spatial dipole via differential cosmic-web gravitational lensing ($\Delta\Phi \sim 10^{-5}$) and predicts a vacuum magnetic permeability shift of $\Delta a_\mu^{\text{TBP}} = a_\mu^{(1)}(\pi/\sqrt{2})\langle\Delta\alpha_{emp}^{-1}\rangle = 2.27 \times 10^{-9}$, consistent within $0.46\sigma$ of the empirical Fermilab/BNL anomaly ($2.49 \pm 0.48 \times 10^{-9}$).

### Pipeline 5: SPARC Galactic Rotation Curves & $24/88$ Cross-Manifold Torque Audit
* **Script:** `tbp_sparc_torque_auditor_v2.py`
* **Data Provenance:** Evaluates $N = 280$ radially resolved rotation curve measurements across 10 benchmark galaxies from the **SPARC (Spitzer Photometry and Accurate Rotation Curves)** database (`J/AJ/152/157`), spanning HSB spirals (NGC2841, NGC7814) to gas-rich dwarfs (NGC3741, DDO154).
* **Key Findings:**
  1. **Ab Initio Radial Acceleration Scale ($a_{TBP}$):** Projecting the cosmic horizon acceleration ($cH_0$, Planck $H_0 = 67.4\text{ km/s/Mpc}$) through the active matter node fraction ($17/88$) and symmetry-breaking ratio ($17/18$) yields $a_{TBP} = cH_0(289/1584) = 1.1947 \times 10^{-10}\text{ m/s}^2$, converging to $99.56\%$ ($0.26\sigma$) of the empirical SPARC RAR scale $g_\dagger = (1.20 \pm 0.02) \times 10^{-10}\text{ m/s}^2$.
  2. **Zero-Parameter Rotation Curves & Cosmic Invariant Convergence:** Reduces the pure 4D Newtonian RAR discrepancy from $0.6528\text{ dex}$ to $0.0183\text{ dex}$ ($r = 0.99900$, mean velocity $\text{RMS} = 2.38\text{ km/s}$), while the ensemble mean outer-halo torque-to-baryon mass ratio $\langle\Delta M_{torque}/M_{bar}\rangle = 5.76 \pm 0.83$ converges within $0.29\sigma$ of the TBP cosmic partition $\Omega_{DM}/\Omega_{OM} = 24/4 = d = 6.00$.

### Pipeline 6: Neutrino Mass Dilution ($R_{crit}^2$), Chiral PMNS Mixing & IceCube HESE Audit
* **Script:** `tbp_neutrino_icecube_auditor_v2.py`
* **Data Provenance:** Audits direct tritium $\beta$-decay kinematic bounds (Mainz, Troitsk, KATRIN KNM1--5), PDG 2024 / NuFIT 5.3 global oscillation parameters, and the 7.5-year IceCube High-Energy Starting Events (HESE) catalog ($60\text{ TeV}\text{--}2.2\text{ PeV}$).
* **Key Findings:**
  1. **Sub-eV Mass Dilution & KATRIN Detection Window:** Diluting weak-vertex rotational friction $\delta_{Final(n)}$ by $R_{crit}^2 = (1936/9)^2 \approx 46272.7901\text{ GeV}^2$ establishes an electron-neutrino ground state at $M_{\nu_e} = 0.25687\text{ eV}$ (satisfying the KATRIN 2024 bound $<0.45\text{ eV}$ with a $+0.1931\text{ eV}$ margin and sitting inside KATRIN's final $0.20\text{ eV}$ sensitivity zone) alongside local vertex thresholds $M_{\nu_\mu}(0) = 2.94635\text{ eV}$ and $M_{\nu_\tau}(0) = 23.57471\text{ eV}$ that relax exponentially via $\Box_{dis}\Phi_k + (\Gamma_{TBP}/\lambda_{coh})u^\mu\nabla_\mu\Phi_k = 0$.
  2. **Zero-Parameter Chiral PMNS Angles:** Derives all three leptonic mixing angles *ab initio*: $\sin^2\theta_{13}^{\text{TBP}} = \frac{1}{2}\Omega_{OM}^{\text{TBP}} = \frac{17}{792}(1.036) \approx 0.02224$ ($0.34\sigma$ from PDG $0.02200 \pm 0.00069$), $\sin^2\theta_{12}^{\text{TBP}} = \frac{1}{3}(1 - \frac{4}{54}) = \frac{25}{81} \approx 0.30864$ ($0.13\sigma$ from $0.3070 \pm 0.0130$), and $\sin^2\theta_{23}^{\text{TBP}} = \frac{1}{2} + B_Z = \frac{6}{11} \approx 0.54545$ ($0.03\sigma$ from $0.5460 \pm 0.0210$).
  3. **IceCube HESE PeV Spectral Softening:** Non-bulk chiral impedance ($\Delta\gamma_{\text{TBP}} = \frac{84}{88}\frac{17}{18}(0.964) = +0.8691$) shifts the Fermi shock baseline ($\gamma_0 = 2.000$) to $\gamma_{\text{TBP}} = 2.8691$, matching the empirical IceCube HESE index ($2.87 \pm 0.20$, $<0.01\sigma$, reducing $\chi_\nu^2$ from $133.39$ to $0.16$) and yielding an Earth flavor ratio of $(0.322 : 0.341 : 0.337)$ at the orthogonal chiral phase $\delta_{CP} = 3\pi/2$.

### Pipeline 7: Multi-Messenger Black Hole Pre-Leakage Horizon Boundary Audit
* **Script:** `tbp_blackhole_horizon_auditor_v2.py`
* **Data Provenance:** Synthesizes observations from Fermi-LAT 4FGL-DR4 & Cherenkov observatories (MAGIC, VERITAS, H.E.S.S.), inner accretion disk X-ray spectroscopy (**XRISM Resolve 2024**, XMM-Newton, NuSTAR), IceCube AGN point sources (NGC 1068, TXS 0506+056), and LIGO GWTC-1 post-merger ringdown records.
* **Key Findings:**
  1. **Intrinsic VHE Blazar Jet Breaks at $R_{crit} = 215.111\text{ GeV}$:** Across 10 benchmark blazars ($0.031 \le z \le 0.536$), transforming EBL-deabsorbed break energies to the black hole rest frame ($E_{\text{brk}}^{\text{int}} = E_{\text{brk}}^{\text{obs}}[1+z]$) collapses the redshift drift onto a weighted mean of $215.17 \pm 5.41\text{ GeV}$ ($0.01\sigma$ from $R_{crit} = 215.11\text{ GeV}$) with mean spectral steepening $\langle\Delta\Gamma\rangle = +0.891 \pm 0.016$ (vs. $\Delta\gamma_{\text{TBP}} = +0.8691$).
  2. **Strong-Field Inner Disk Vacuum Microlensing ($\text{Fe K}\alpha$):** Evaluating $\alpha_{obs}^{-1}(r) = 137.036000 - \Sigma\Phi_{local}(r)$ with the exact 2PN effective potential $\Sigma\Phi_{local}(r) = 1/r + 1/(2r^2)$ ($r \equiv R_{\text{in}}/R_g$) accounts for the apparent supersolar iron abundance anomaly via $A_{\text{Fe}}^{\text{TBP}}(r) = 1 + 12\,\Sigma\Phi_{local}(r)$ across NGC 4151, Cygnus X-1, Mrk 335, GX 339-4, MCG--6-30-15, and 1H 0707--495 within $0.17\sigma\text{--}0.69\sigma$ ($\chi_\nu^2 = 0.30, r = 0.99861$).
  3. **Coronal Neutrino Conduit ($1/L_{QCD} = 62.33\times$) & Unitary Ringdown:** 1-loop QCD confinement ($L_{QCD} = 3/187$) suppresses hadronic gamma-ray escape above $R_{crit}$ while channeling momentum into neutrinos, yielding a coronal neutrino-to-gamma flux ratio of $1/L_{QCD} = 62.33\times$ (matching NGC 1068 at $\sim 50\text{--}65\times$) and $29.44\times$ in blazar jets (TXS 0506+056 at $\sim 25\text{--}35\times$), while post-merger LIGO residuals remain strictly below $H_{lim} = 1.894 \times 10^{-3}$.

---

## 5. Popperian Falsifiability Matrix

Because TBP commits to exact rational fractions with zero adjustable parameters, every sectoral derivation is exposed to decisive empirical refutation:

| Operational Sector | Rigid TBP Prediction | Target Facility / Experiment | Exact Empirical Falsification Condition |
| :--- | :--- | :--- | :--- |
| **1. Direct Neutrino Mass** | $M_{\nu_e} = 0.25687\text{ eV}$ | KATRIN Final ($0.20\text{ eV}$) & Project 8 | A $90\%$ CL upper bound of $M_{\nu_e} < 0.24\text{ eV}$ without a tritium $\beta$-decay endpoint distortion at $0.257\text{ eV}$. |
| **2. PMNS Mixing Angles** | $\sin^2\theta_{23} = 6/11 \approx 0.54545$<br>$\sin^2\theta_{12} = 25/81 \approx 0.30864$ | JUNO, DUNE & Hyper-Kamiokande | Confirmation of lower-octant $\theta_{23}$ ($\sin^2\theta_{23} < 0.50$) or deviation from $0.30864$ by $>5\sigma$. |
| **3. Strong-Field $\alpha^{-1}(r)$** | $\alpha_{obs}^{-1}(r) = 137.036000 - \Sigma\Phi_{local}(r)$ | XRISM Resolve (Gate-Valve-Open), Solar Probes & ELT | Failure to observe $+0.53\%$ to $+1.72\%$ $\text{Fe XXVI Ly}\alpha$ fine-structure doublet widening at $R_{\text{in}} \le 6\,R_g$. |
| **4. Vacuum Ceiling ($R_{crit}$)** | $\Theta(1 - \mathcal{M}_i^2/R_{crit}^2) \to 0$ at $215.111\text{ GeV}$ | LHC Run-3 / HL-LHC & CTA Blazar Array | Discovery of a fundamental 4th-generation fermion $>215.111\text{ GeV}$ or intrinsic blazar breaks deviating from $215.11\text{ GeV}$ by $>5\sigma$. |
| **5. GW Low-Pass Filter** | $\mathcal{R}_{att} < H_{lim} = 1/528 \approx 1.894\times 10^{-3}$ | LIGO-Virgo-KAGRA O4/O5 ($\text{SNR}_{Net} > 20$) | Observation of $>5\sigma$ lattice jitter above $528\text{ Hz}$ or post-merger echoes $\mathcal{R}_{\text{echo}} > 1.894\times 10^{-3}$ in heavy BBH mergers. |

---

## 6. Installation, Dependencies & Execution

### Prerequisites
Clone the repository and install the required scientific Python packages (Python 3.10+ recommended):

    git clone https://github.com/cefiyana-clover/TBP-Computational-Framework.git
    cd TBP-Computational-Framework
    pip install -r requirements.txt

### Standard Dependencies (`requirements.txt`)

    numpy>=1.24.0
    scipy>=1.10.0
    sympy>=1.12
    matplotlib>=3.7.0
    h5py>=3.8.0
    requests>=2.31.0

*(Note: Running the raw CERN CMS `.root` extraction pipeline `tbp_cms_omni_extractor.py` from scratch additionally requires `uproot>=5.0`, `awkward>=2.0`, and `coffea>=0.7`.)*

### Executing the Verification & Auditing Suite
Each script prints a structured statistical audit log to the console and exports a 600 DPI publication-grade synthesis figure:

    # 1. Foundational Symbolic & Rational Constant Verification (Paper I)
    python tbp_computational_evaluation.py

    # 2. CERN CMS 13 TeV Collider Audits at 215.111 GeV & 173.111 GeV (Papers II & III)
    python tbp_cms_215_auditor.py
    python tbp_cms_173_top_auditor.py

    # 3. LIGO GWTC-1 Gravitational-Wave Low-Pass Filter Audit (Paper IV)
    python tbp_ligo_omni_auditor_v2.py

    # 4. Terrestrial Vacuum Microlensing, Quasar Dipoles & Muon g-2 Audit (Paper V)
    python tbp_metrology_quasar_auditor_v2.py

    # 5. SPARC Galactic Rotation Curves & 24/88 Cross-Manifold Torque Audit (Paper VI)
    python tbp_sparc_torque_auditor_v2.py

    # 6. Neutrino Mass Dilution, PMNS Mixing & IceCube HESE Audit (Paper VII)
    python tbp_neutrino_icecube_auditor_v2.py

    # 7. Multi-Messenger Black Hole Pre-Leakage Horizon Boundary Audit (Paper VIII)
    python tbp_blackhole_horizon_auditor_v2.py

---

## 7. Citation & License

All custom Python source code in this repository is released under the **GNU General Public License v3.0 (GPLv3)**.

If you utilize this computational framework or reference the Connecting Thread Theory in your research, please cite the foundational theoretical manuscript and relevant companion phenomenological audits:

1. **Cefiyana, C.** (2026). *"The Connecting Thread Theory: An Ab Initio 6D Unified Geometric Architecture for Quantum and Gravitational Metrics."* SSRN Preprint, `doi:10.2139/ssrn.6953858` / Zenodo, `doi:10.5281/zenodo.19978186`.
2. **Cefiyana, C.** (2026). *"The Connecting Thread Theory: Multi-Sector Phenomenological Evaluation of Topological Deformation and Kinematic Divergence at the 215.111 GeV Critical Vacuum Ceiling."* SSRN Preprint, `doi:10.2139/ssrn.7075779`.
3. **Cefiyana, C.** (2026). *"Phenomenological Audit of the 173.111 GeV Macroscopic Tensor Stabilization Threshold: Lepton-Specific Neutrino Kinematics, Trigger Turn-On Diagnostics, and Sub-Ceiling Euclidean Rigidity in 13 TeV CMS Open Data."* SSRN Preprint.
4. **Cefiyana, C.** (2026). *"Phenomenological Audit of Gravitational-Wave Topological Low-Pass Filtering across the GWTC-1 Catalog: High-Frequency Lattice Jitter Suppression and Parseval Spectral Attenuation Bounds in LIGO Open Data."* SSRN Preprint.
5. **Cefiyana, C.** (2026). *"Phenomenological Audit of Terrestrial Vacuum Microlensing at the $\alpha_{ideal}^{-1} = 137.036000$ Electromagnetic Anchor: CODATA Metrology Gaps, VLT/Keck Quasar Dipoles, and Muon $(g-2)$ Permeability Scaling."* SSRN Preprint.
6. **Cefiyana, C.** (2026). *"Phenomenological Audit of Cross-Manifold Topological Torque ($24/88$ Nodes) across SPARC Galactic Rotation Curves: Ab Initio Derivation of the Radial Acceleration Scale $a_{TBP} = cH_0(289/1584)$ and Halo Torque Saturation."* SSRN Preprint.
7. **Cefiyana, C.** (2026). *"Phenomenological Audit of Neutrino Mass Dilution ($R_{crit}^2$), Chiral Boundary PMNS Mixing Angles, and IceCube PeV Spectral Attenuation within the 6D Connecting Thread Theory."* SSRN Preprint.
8. **Cefiyana, C.** (2026). *"Multi-Messenger Phenomenological Audit of Black Hole Pre-Leakage Boundary States: Relativistic Blazar Jet Spectral Breaks at 215.111 GeV, Inner-Disk Fe K-alpha Vacuum Microlensing, Coronal Neutrino Conduits, and Unitary Ringdown Bounds."* SSRN Preprint.

---

**Author:** Cefiyana  
**Affiliation:** Independent Interdisciplinary Researcher, Indonesia  
**Email:** [cefiyana@fiveleafclover.org](mailto:cefiyana@fiveleafclover.org)  
**ORCID iD:** [https://orcid.org/0009-0008-4324-9515](https://orcid.org/0009-0008-4324-9515)
