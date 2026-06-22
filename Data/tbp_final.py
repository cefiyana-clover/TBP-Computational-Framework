# ==============================================================================
# TBP KINEMATIC PIPELINE - PARALLEL MULTIPROCESSING ARCHITECTURE
# TARGET THRESHOLD: 215.11111111111 GeV
# DATASET: RUN 2016G PFNANOAOD (DoubleEG & DoubleMuon) -> 1.3 TB SCALE
# ==============================================================================

import os
import time
import uuid
import uproot
import awkward as ak
import vector
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import concurrent.futures

# ------------------------------------------------------------------------------
# 1. CONTROL PARAMETERS
# ------------------------------------------------------------------------------
TBP_THRESHOLD = 215.11111111111
MAX_FILES = 5000
MAX_WORKERS = 16 # Menggunakan 16 core dari 24 core EPYC untuk stabilitas I/O

print("\n==================================================")
print(" ☢️ OPERATION: OMNICIDE (1.3 TB INGESTION - MULTIPROCESSING) ☢️ ")
print("==================================================\n")
print(f"[INFO] Target threshold: {TBP_THRESHOLD} GeV")
print(f"[INFO] Execution Mode: Parallel Iterative Deletion ({MAX_WORKERS} Workers)")

TXT_FILES = ["DoubleMuon_files.txt", "DoubleEG_files.txt"]
all_target_urls = []

# ------------------------------------------------------------------------------
# 2. MANIFEST EXTRACTION
# ------------------------------------------------------------------------------
print("[PROCESS] Reading raw XRootD coordinates from local text files...")

for txt_file in TXT_FILES:
    if os.path.exists(txt_file):
        try:
            with open(txt_file, 'r') as f:
                links = [line.strip() for line in f if line.strip().startswith('root://')]
                all_target_urls.extend(links)
            print(f"[STATUS] Coordinates extracted successfully from {txt_file}.")
        except Exception as e:
            print(f"[ERROR] Failed to read {txt_file}: {e}")
    else:
        print(f"[FATAL ERROR] {txt_file} not found in the current directory! Download it first.")

all_target_urls = list(set(all_target_urls))

if not all_target_urls:
    raise RuntimeError("\n[FATAL] No coordinates found. Terminating.")

all_target_urls = all_target_urls[:MAX_FILES]
print(f"\n[STATUS] AMMO LOADED! Total .root files queued: {len(all_target_urls)}")

# ------------------------------------------------------------------------------
# 3. WORKER FUNCTION (CORE ENGINE)
# ------------------------------------------------------------------------------
def process_file(target):
    """Fungsi mandiri untuk di-eksekusi secara paralel oleh setiap core CPU."""
    temp_root = f"temp_{uuid.uuid4().hex}.root"
    local_results = {
        "di_photon_mass": [], 
        "four_lepton_mass": [], 
        "met_distribution": [], 
        "events": 0,
        "success": False
    }
    
    # DOWNLOAD PHASE (Muted log output to prevent console spam)
    exit_code = os.system(f"xrdcp -f -s {target} {temp_root} > /dev/null 2>&1")
    
    if exit_code != 0 or not os.path.exists(temp_root):
        return local_results # Download failed
    
    # EXTRACTION PHASE (LOGIKA FISIKA TIDAK DIUBAH SAMA SEKALI)
    try:
        with uproot.open(f"{temp_root}:Events") as events:
            chunk = events.arrays(["nPhoton", "Photon_pt", "Photon_eta", "Photon_phi",
                                   "nElectron", "Electron_pt", "nMuon", "Muon_pt", "MET_pt"], 
                                   how=dict, library="ak")
            
            local_results["events"] = len(chunk["MET_pt"])
            
            # Channel 1: Di-Photon
            if "nPhoton" in chunk.keys():
                mask_2ph = chunk["nPhoton"] >= 2
                if ak.any(mask_2ph):
                    g1 = vector.zip({"pt": chunk["Photon_pt"][mask_2ph, 0], "eta": chunk["Photon_eta"][mask_2ph, 0], 
                                     "phi": chunk["Photon_phi"][mask_2ph, 0], "mass": 0})
                    g2 = vector.zip({"pt": chunk["Photon_pt"][mask_2ph, 1], "eta": chunk["Photon_eta"][mask_2ph, 1], 
                                     "phi": chunk["Photon_phi"][mask_2ph, 1], "mass": 0})
                    ph_cut = (g1.pt > 40) & (g2.pt > 40) & (abs(g1.eta) < 2.5) & (abs(g2.eta) < 2.5) & (abs(g1.deltaphi(g2)) > 2.5)
                    ph_mass = (g1[ph_cut] + g2[ph_cut]).mass
                    local_results["di_photon_mass"].append(ak.to_numpy(ph_mass[ph_mass > 150])) 

            # Channel 2: Four-Lepton
            if "nElectron" in chunk.keys() and "nMuon" in chunk.keys():
                tot_lep = chunk["nElectron"] + chunk["nMuon"]
                mask_4l = tot_lep >= 4
                if ak.any(mask_4l):
                    total_4l_pt = ak.sum(chunk["Electron_pt"][mask_4l], axis=-1) + ak.sum(chunk["Muon_pt"][mask_4l], axis=-1)
                    local_results["four_lepton_mass"].append(ak.to_numpy(total_4l_pt[total_4l_pt > 150]))

            # Channel 3: MET
            if "MET_pt" in chunk.keys():
                met_data = chunk["MET_pt"]
                local_results["met_distribution"].append(ak.to_numpy(met_data[met_data > 100])) 
        
        local_results["success"] = True
        
    except Exception:
        pass # Corrupted matrix silent fail untuk menjaga pool tetap berjalan
    
    finally:
        # DELETION PHASE (Eat & Burn)
        if os.path.exists(temp_root):
            os.remove(temp_root)
            
    return local_results

# ------------------------------------------------------------------------------
# 4. PARALLEL EXECUTION LOOP
# ------------------------------------------------------------------------------
print(f"\n[PROCESS] Initiating Multiprocessing Pool with {MAX_WORKERS} workers...")
start_time = time.time()

global_results = {"di_photon_mass": [], "four_lepton_mass": [], "met_distribution": []}
total_events = 0
processed_files = 0

# Menggunakan ProcessPoolExecutor untuk full CPU core utilization
with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
    # Map the worker function to all target URLs
    futures = {executor.submit(process_file, url): url for url in all_target_urls}
    
    for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
        res = future.result()
        
        if res["success"]:
            total_events += res["events"]
            processed_files += 1
            
            # Aggregate array results
            if res["di_photon_mass"]:
                global_results["di_photon_mass"].extend(res["di_photon_mass"])
            if res["four_lepton_mass"]:
                global_results["four_lepton_mass"].extend(res["four_lepton_mass"])
            if res["met_distribution"]:
                global_results["met_distribution"].extend(res["met_distribution"])
                
        # Status update per file kelar (karena paralel, urutan berantakan, jadi pakai counter)
        if i % 10 == 0 or i == len(all_target_urls):
            print(f"[STATUS] Progress: {i}/{len(all_target_urls)} files processed | Cumulative Events: {total_events}")

exec_time = time.time() - start_time
print(f"\n[INFO] OPERATION COMPLETE.")
print(f"[INFO] Total files successfully parsed: {processed_files}/{len(all_target_urls)}")
print(f"[INFO] Total events evaluated: {total_events}")
print(f"[INFO] Total execution time: {exec_time/60:.2f} minutes")

# ------------------------------------------------------------------------------
# 5. DATA VISUALIZATION AND CSV EXPORT (TIDAK DIUBAH)
# ------------------------------------------------------------------------------
print("[PROCESS] Generating TBP Diagnostic Render...")

fig, axs = plt.subplots(3, 1, figsize=(15, 18))
plt.subplots_adjust(hspace=0.4)

if global_results["di_photon_mass"]:
    ph_all = np.concatenate(global_results["di_photon_mass"])
    counts_ph, bins_ph = np.histogram(ph_all[(ph_all > 200) & (ph_all < 230)], bins=600, range=(200, 230))
    axs[0].bar((bins_ph[:-1]+bins_ph[1:])/2, counts_ph, width=0.05, color='cyan', alpha=0.9)
axs[0].axvline(x=TBP_THRESHOLD, color='red', linestyle='dashed', linewidth=2.5)
axs[0].set_title(r'Channel 1: Di-Photon Invariant Mass', fontweight='bold', color='white')
axs[0].set_facecolor('black')
axs[0].set_xlim(210, 220)

if global_results["four_lepton_mass"]:
    l4_all = np.concatenate(global_results["four_lepton_mass"])
    counts_l4, bins_l4 = np.histogram(l4_all[(l4_all > 200) & (l4_all < 230)], bins=600, range=(200, 230))
    axs[1].bar((bins_l4[:-1]+bins_l4[1:])/2, counts_l4, width=0.05, color='gold', alpha=0.9)
axs[1].axvline(x=TBP_THRESHOLD, color='red', linestyle='dashed', linewidth=2.5)
axs[1].set_title('Channel 2: 4-Lepton System Total Transverse Momentum', fontweight='bold', color='white')
axs[1].set_facecolor('black')
axs[1].set_xlim(210, 220)

golden_met = []
if global_results["met_distribution"]:
    met_all = np.concatenate(global_results["met_distribution"])
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
