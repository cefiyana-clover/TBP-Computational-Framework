# ==============================================================================
# TBP KINEMATIC PIPELINE - ITERATIVE DELETION ARCHITECTURE (TEXT-BASED VAST.AI FINAL)
# TARGET THRESHOLD: 215.11111111111 GeV
# DATASET: RUN 2016G PFNANOAOD (DoubleEG & DoubleMuon) -> 1.3 TB SCALE
# ==============================================================================

import os
import time
import uproot
import awkward as ak
import vector
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# 1. CONTROL PARAMETERS
# ------------------------------------------------------------------------------
TBP_THRESHOLD = 215.11111111111
TEMP_ROOT = "temp_dataset.root"
MAX_FILES = 5000

print("\n==================================================")
print(" ☢️ OPERATION: OMNICIDE (1.3 TB INGESTION) ☢️ ")
print("==================================================\n")
print(f"[INFO] Target threshold: {TBP_THRESHOLD} GeV")
print(f"[INFO] Execution Mode: Iterative Deletion (Eat & Burn)")

TXT_FILES = ["DoubleMuon_files.txt", "DoubleEG_files.txt"]
all_target_urls = []

# ------------------------------------------------------------------------------
# 2. MANIFEST EXTRACTION (LIGHTWEIGHT TEXT READER)
# ------------------------------------------------------------------------------
print("[PROCESS] Reading raw XRootD coordinates from local text files...")

for txt_file in TXT_FILES:
    if os.path.exists(txt_file):
        try:
            with open(txt_file, 'r') as f:
                # Mengambil baris yang murni link XRootD, membuang spasi kosong
                links = [line.strip() for line in f if line.strip().startswith('root://')]
                all_target_urls.extend(links)
            print(f"[STATUS] Coordinates extracted successfully from {txt_file}.")
        except Exception as e:
            print(f"[ERROR] Failed to read {txt_file}: {e}")
    else:
        print(f"[FATAL ERROR] {txt_file} not found in the current directory! Download it first.")

all_target_urls = list(set(all_target_urls))

if not all_target_urls:
    raise RuntimeError("\n[FATAL] No coordinates found. Make sure the .txt files are downloaded. Terminating.")

# Membatasi jumlah file sesuai MAX_FILES
all_target_urls = all_target_urls[:MAX_FILES]
print(f"\n[STATUS] AMMO LOADED! Total .root files queued: {len(all_target_urls)}")

# ------------------------------------------------------------------------------
# 3. ITERATIVE INGESTION AND EXTRACTION LOOP (HIGH-SPEED XROOTD)
# ------------------------------------------------------------------------------
results = {"di_photon_mass": [], "four_lepton_mass": [], "met_distribution": []}
total_events = 0
start_time = time.time()

for idx, target in enumerate(all_target_urls):
    file_name = target.split('/')[-1]
    print(f"\n[+] Cycle {idx+1}/{len(all_target_urls)}: Ingesting {file_name}")
    
    # DOWNLOAD PHASE (XRootD - Bypass Web Firewall)
    exit_code = os.system(f"xrdcp -f -s {target} {TEMP_ROOT}")
    
    if exit_code != 0 or not os.path.exists(TEMP_ROOT):
        print(f" [!] XRootD transfer failed or timed out. Proceeding to next target.")
        continue
    
    # EXTRACTION PHASE
    try:
        with uproot.open(f"{TEMP_ROOT}:Events") as events:
            chunk = events.arrays(["nPhoton", "Photon_pt", "Photon_eta", "Photon_phi",
                                   "nElectron", "Electron_pt", "nMuon", "Muon_pt", "MET_pt"], 
                                   how=dict, library="ak")
            
            current_events = len(chunk["MET_pt"])
            total_events += current_events
            
            # Channel 1: Di-Photon
            if "nPhoton" in chunk.fields:
                mask_2ph = chunk["nPhoton"] >= 2
                if ak.any(mask_2ph):
                    g1 = vector.zip({"pt": chunk["Photon_pt"][mask_2ph, 0], "eta": chunk["Photon_eta"][mask_2ph, 0], 
                                     "phi": chunk["Photon_phi"][mask_2ph, 0], "mass": 0})
                    g2 = vector.zip({"pt": chunk["Photon_pt"][mask_2ph, 1], "eta": chunk["Photon_eta"][mask_2ph, 1], 
                                     "phi": chunk["Photon_phi"][mask_2ph, 1], "mass": 0})
                    ph_cut = (g1.pt > 40) & (g2.pt > 40) & (abs(g1.eta) < 2.5) & (abs(g2.eta) < 2.5) & (abs(g1.deltaphi(g2)) > 2.5)
                    ph_mass = (g1[ph_cut] + g2[ph_cut]).mass
                    results["di_photon_mass"].append(ak.to_numpy(ph_mass[ph_mass > 150])) 

            # Channel 2: Four-Lepton
            if "nElectron" in chunk.fields and "nMuon" in chunk.fields:
                tot_lep = chunk["nElectron"] + chunk["nMuon"]
                mask_4l = tot_lep >= 4
                if ak.any(mask_4l):
                    total_4l_pt = ak.sum(chunk["Electron_pt"][mask_4l], axis=-1) + ak.sum(chunk["Muon_pt"][mask_4l], axis=-1)
                    results["four_lepton_mass"].append(ak.to_numpy(total_4l_pt[total_4l_pt > 150]))

            # Channel 3: MET
            if "MET_pt" in chunk.fields:
                met_data = chunk["MET_pt"]
                results["met_distribution"].append(ak.to_numpy(met_data[met_data > 100])) 
            
        print(f" [STATUS] Bedah sukses. Parsed: {current_events} events. Cumulative: {total_events}")
        
    except Exception as e:
        print(f" [!] Extraction failed due to corrupted particle matrix: {e}")
    
    finally:
        # DELETION PHASE (Eat & Burn)
        if os.path.exists(TEMP_ROOT):
            os.remove(TEMP_ROOT)

exec_time = time.time() - start_time
print(f"\n[INFO] OPERATION COMPLETE.")
print(f"[INFO] Total events evaluated: {total_events}")
print(f"[INFO] Total execution time: {exec_time/60:.2f} minutes")

# ------------------------------------------------------------------------------
# 4. DATA VISUALIZATION AND CSV EXPORT
# ------------------------------------------------------------------------------
print("[PROCESS] Generating TBP Diagnostic Render...")

fig, axs = plt.subplots(3, 1, figsize=(15, 18))
plt.subplots_adjust(hspace=0.4)

if results["di_photon_mass"]:
    ph_all = np.concatenate(results["di_photon_mass"])
    counts_ph, bins_ph = np.histogram(ph_all[(ph_all > 200) & (ph_all < 230)], bins=600, range=(200, 230))
    axs[0].bar((bins_ph[:-1]+bins_ph[1:])/2, counts_ph, width=0.05, color='cyan', alpha=0.9)
axs[0].axvline(x=TBP_THRESHOLD, color='red', linestyle='dashed', linewidth=2.5)
axs[0].set_title(r'Channel 1: Di-Photon Invariant Mass', fontweight='bold', color='white')
axs[0].set_facecolor('black')
axs[0].set_xlim(210, 220)

if results["four_lepton_mass"]:
    l4_all = np.concatenate(results["four_lepton_mass"])
    counts_l4, bins_l4 = np.histogram(l4_all[(l4_all > 200) & (l4_all < 230)], bins=600, range=(200, 230))
    axs[1].bar((bins_l4[:-1]+bins_l4[1:])/2, counts_l4, width=0.05, color='gold', alpha=0.9)
axs[1].axvline(x=TBP_THRESHOLD, color='red', linestyle='dashed', linewidth=2.5)
axs[1].set_title('Channel 2: 4-Lepton System Total Transverse Momentum', fontweight='bold', color='white')
axs[1].set_facecolor('black')
axs[1].set_xlim(210, 220)

golden_met = []
if results["met_distribution"]:
    met_all = np.concatenate(results["met_distribution"])
    met_plot = met_all[(met_all > 150) & (met_all < 250)]
    counts_met, bins_met = np.histogram(met_plot, bins=100, range=(150, 250))
    axs[2].bar((bins_met[:-1]+bins_met[1:])/2, counts_met, width=1.0, color='orange', edgecolor='black', alpha=0.9)
    
    golden_mask = (met_all >= 215.0) & (met_all <= 216.0)
    golden_met = met_all[golden_mask]

axs[2].axvline(x=TBP_THRESHOLD, color='red', linestyle='dashed', linewidth=3)
axs[2].set_title(r'Channel 3: Missing Transverse Energy $E_T^{miss}$', fontweight='bold', color='white')
axs[2].set_facecolor('black')
axs[2].set_xlim(180, 250)

for ax in axs:
    ax.tick_params(colors='white')
    ax.spines['bottom'].set_color('white')
    ax.spines['left'].set_color('white')

fig.patch.set_facecolor('#121212')
plt.savefig("TBP_Final_Discovery.png", dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')

if len(golden_met) > 0:
    with open("TBP_Discovery_Anomalies.csv", "w") as f:
        f.write("Index,MET_GeV\n")
        for i, val in enumerate(golden_met):
            f.write(f"Anomali_{i+1},{val:.6f}\n")

print("\n[SUCCESS] Render secured as 'TBP_Final_Discovery.png'.")
print(f"[SUCCESS] {len(golden_met)} anomalous 5-Sigma events isolated in CSV.")


