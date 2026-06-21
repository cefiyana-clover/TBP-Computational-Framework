# ==============================================================================
# AUTOMATED KINEMATIC PIPELINE (TBP FRAMEWORK) - HEADLESS CLI MODE
# TARGET: MULTI-CHANNEL EVALUATION AT 215.11111111111 GeV
# ==============================================================================

import os
import time
import uproot
import awkward as ak
import vector
import numpy as np
import matplotlib
# Set backend to 'Agg' for headless server execution (no GUI required)
matplotlib.use('Agg') 
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# 1. PARAMETER INITIALIZATION
# ------------------------------------------------------------------------------
TBP_THRESHOLD = 215.11111111111
CHUNK_SIZE = 100000

print("[INFO] Initializing TBP Kinematic Pipeline (CLI Headless Mode).")
print(f"[INFO] Target threshold defined at: {TBP_THRESHOLD} GeV.")

# ------------------------------------------------------------------------------
# 2. AUTOMATED DATA ACQUISITION (RUN 1 & RUN 2 INTEGRATION)
# ------------------------------------------------------------------------------
# Integrasi data Run 1 (8 TeV) dan Run 2 (13 TeV) CMS Open Data
DATA_SOURCES = {
    # --- CALIBRATION BASELINE (RUN 1 - 8 TeV) ---
    "Run2012B_DoubleElectron.root": "https://opendata.cern.ch/record/12367/files/Run2012B_DoubleElectron.root",
    "Run2012BC_DoubleMuParked.root": "https://eospublic.cern.ch/eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/Run2012BC_DoubleMuParked_Muons.root",
    
    # --- HIGH ENERGY TARGETS (RUN 2 - 13 TeV) ---
    # DoubleEG menggabungkan trigger Electron dan Photon (Gamma)
    "Run2015D_DoubleEG.root": "https://opendata.cern.ch/record/24360/files/Run2015D_DoubleEG.root",
    "Run2015D_DoubleMuon.root": "https://opendata.cern.ch/record/24361/files/Run2015D_DoubleMuon.root"
}

data_files = []

print("[INFO] Executing automated data acquisition protocol...")
for filename, url in DATA_SOURCES.items():
    if not os.path.exists(filename) or os.path.getsize(filename) < 1000000:
        print(f"[PROCESS] Downloading {filename}...")
        # '-c' ditambahkan agar wget melanjutkan unduhan jika terputus (resume)
        os.system(f"wget -c -q --show-progress --no-check-certificate {url} -O {filename}")
    
    if os.path.exists(filename) and os.path.getsize(filename) > 1000000:
        data_files.append(filename)
        print(f"[STATUS] {filename} acquired and verified.")
    else:
        print(f"[WARNING] Failed to acquire {filename}. Bypassing to next target.")

if not data_files:
    raise FileNotFoundError("[FATAL] No valid root files available for processing. Terminating execution.")

# ------------------------------------------------------------------------------
# 3. KINEMATIC EXTRACTION PIPELINE
# ------------------------------------------------------------------------------
results = {
    "di_photon_mass": [],
    "four_lepton_mass": [],
    "met_distribution": []
}

total_events = 0
start_time = time.time()

print(f"[INFO] Initiating iterative extraction on {len(data_files)} files (Chunk size: {CHUNK_SIZE})...")

for chunk in uproot.iterate([f + ":Events" for f in data_files], step_size=CHUNK_SIZE):
    total_events += len(chunk)
    
    # --- CHANNEL 1: DI-PHOTON INVARIANT MASS ---
    if "nPhoton" in chunk.fields:
        mask_2ph = chunk["nPhoton"] >= 2
        if ak.any(mask_2ph):
            g1 = vector.zip({"pt": chunk["Photon_pt"][mask_2ph, 0], "eta": chunk["Photon_eta"][mask_2ph, 0], 
                             "phi": chunk["Photon_phi"][mask_2ph, 0], "mass": 0})
            g2 = vector.zip({"pt": chunk["Photon_pt"][mask_2ph, 1], "eta": chunk["Photon_eta"][mask_2ph, 1], 
                             "phi": chunk["Photon_phi"][mask_2ph, 1], "mass": 0})
            
            ph_cut = (g1.pt > 40) & (g2.pt > 40) & (abs(g1.eta) < 2.5) & (abs(g2.eta) < 2.5) & (abs(g1.deltaphi(g2)) > 2.5)
            ph_mass = (g1[ph_cut] + g2[ph_cut]).mass
            results["di_photon_mass"].append(ak.to_numpy(ph_mass))

    # --- CHANNEL 2: FOUR-LEPTON SYSTEM ---
    if "nElectron" in chunk.fields and "nMuon" in chunk.fields:
        tot_lep = chunk["nElectron"] + chunk["nMuon"]
        mask_4l = tot_lep >= 4
        
        if ak.any(mask_4l):
            e_pt = ak.sum(chunk["Electron_pt"][mask_4l], axis=-1)
            mu_pt = ak.sum(chunk["Muon_pt"][mask_4l], axis=-1)
            total_4l_pt = e_pt + mu_pt
            results["four_lepton_mass"].append(ak.to_numpy(total_4l_pt[total_4l_pt > 150]))

    # --- CHANNEL 3: MISSING TRANSVERSE ENERGY (MET) ---
    if "MET_pt" in chunk.fields and "nElectron" in chunk.fields and "nMuon" in chunk.fields:
        mask_2l = (chunk["nElectron"] + chunk["nMuon"]) == 2
        
        if ak.any(mask_2l):
            met_data = chunk["MET_pt"][mask_2l]
            results["met_distribution"].append(ak.to_numpy(met_data))

    if total_events % 1000000 == 0:
        print(f"[PROCESS] {total_events} events processed...")

exec_time = time.time() - start_time
print(f"[INFO] Extraction complete. Total events: {total_events}. Execution time: {exec_time:.2f} seconds.")

# ------------------------------------------------------------------------------
# 4. DATA VISUALIZATION (EXPORT TO PNG)
# ------------------------------------------------------------------------------
print("[INFO] Compiling diagnostic plots to image file...")

fig, axs = plt.subplots(3, 1, figsize=(15, 18))
plt.subplots_adjust(hspace=0.4)

# Plot 1: Di-Photon
if results["di_photon_mass"]:
    ph_all = np.concatenate(results["di_photon_mass"])
    ph_plot = ph_all[(ph_all > 200) & (ph_all < 230)]
    counts_ph, bins_ph = np.histogram(ph_plot, bins=600, range=(200, 230))
    axs[0].bar((bins_ph[:-1]+bins_ph[1:])/2, counts_ph, width=0.05, color='cyan', alpha=0.9)
axs[0].axvline(x=TBP_THRESHOLD, color='red', linestyle='dashed', linewidth=2.5)
axs[0].set_title(r'Channel 1: Di-Photon Invariant Mass $M_{\gamma\gamma}$ (Bin: 0.05 GeV)', fontweight='bold', color='white')
axs[0].set_xlabel('Invariant Mass (GeV)', color='white')
axs[0].set_ylabel('Events / 0.05 GeV', color='white')
axs[0].set_xlim(210, 220)
axs[0].set_facecolor('black')

# Plot 2: Four-Lepton
if results["four_lepton_mass"]:
    l4_all = np.concatenate(results["four_lepton_mass"])
    l4_plot = l4_all[(l4_all > 200) & (l4_all < 230)]
    counts_l4, bins_l4 = np.histogram(l4_plot, bins=600, range=(200, 230))
    axs[1].bar((bins_l4[:-1]+bins_l4[1:])/2, counts_l4, width=0.05, color='gold', alpha=0.9)
axs[1].axvline(x=TBP_THRESHOLD, color='red', linestyle='dashed', linewidth=2.5)
axs[1].set_title('Channel 2: 4-Lepton System Total Transverse Momentum (Bin: 0.05 GeV)', fontweight='bold', color='white')
axs[1].set_xlabel('Total $p_T$ (GeV)', color='white')
axs[1].set_ylabel('Events / 0.05 GeV', color='white')
axs[1].set_xlim(210, 220)
axs[1].set_facecolor('black')

# Plot 3: MET
if results["met_distribution"]:
    met_all = np.concatenate(results["met_distribution"])
    met_plot = met_all[(met_all > 150) & (met_all < 250)]
    counts_met, bins_met = np.histogram(met_plot, bins=100, range=(150, 250))
    axs[2].bar((bins_met[:-1]+bins_met[1:])/2, counts_met, width=1.0, color='orange', edgecolor='black', alpha=0.9)
axs[2].axvline(x=TBP_THRESHOLD, color='red', linestyle='dashed', linewidth=3)
axs[2].set_title(r'Channel 3: Missing Transverse Energy $E_T^{miss}$ (Bin: 1.0 GeV)', fontweight='bold', color='white')
axs[2].set_xlabel('Missing Transverse Energy (GeV)', color='white')
axs[2].set_ylabel('Events / 1.0 GeV', color='white')
axs[2].set_xlim(180, 250)
axs[2].set_facecolor('black')

# Global Formatting
for ax in axs:
    ax.tick_params(colors='white')
    ax.spines['bottom'].set_color('white')
    ax.spines['left'].set_color('white')

fig.patch.set_facecolor('#121212')

# EXPORT TO FILE (CRITICAL FOR HEADLESS EXECUTION)
output_filename = "TBP_MultiChannel_Result.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
print(f"[INFO] Rendering complete. Graph saved locally as '{output_filename}'.")
print("[INFO] Execution terminated successfully.")
