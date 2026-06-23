# ==============================================================================
# AUTOMATED KINEMATIC PIPELINE - OMNICIDE V6 (FULL SCALE)
# PROTOCOL: CAUSALITY EXPANSION & DIMENSIONAL LEAKAGE (PHASE 4)
# TARGET: 215.11111111111 GeV | FULL KINEMATIC TRANSVERSE MASS (M_T)
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
import concurrent.futures
import subprocess

# ------------------------------------------------------------------------------
# PARAMETER KONTROL SKALA PENUH
# ------------------------------------------------------------------------------
TBP_THRESHOLD = 215.11111111111
MAX_FILES = 100000   
MAX_WORKERS = 32     

print("\n==================================================")
print(" ☢️ OPERATION: OMNICIDE V6 (FULL KINEMATIC M_T) ☢️ ")
print("==================================================\n")

# ------------------------------------------------------------------------------
# 1. AKUISISI MANIFES OTONOM (API INTEGRATION)
# ------------------------------------------------------------------------------
RECORD_IDS = [31304, 31305, 31306, 31307, 31309, 31313, 31312, 390, 31303, 31308, 31311, 31314, 31310]
all_target_urls = []

print("[PROCESS] Mengakuisisi manifes XRootD langsung dari server CERN Open Data...")
for recid in RECORD_IDS:
    print(f" -> Menarik koordinat untuk Record ID: {recid}...")
    try:
        result = subprocess.run(
            ["cernopendata-client", "get-file-locations", "--recid", str(recid), "--protocol", "xrootd"],
            capture_output=True, text=True, check=True
        )
        urls = result.stdout.strip().split('\n')
        valid_urls = [url.strip() for url in urls if url.strip().startswith('root://')]
        all_target_urls.extend(valid_urls)
        print(f"    [STATUS] Terkunci {len(valid_urls)} file dari Record {recid}.")
    except Exception as e:
        print(f"    [ERROR] Kegagalan akses Record {recid}. Detail: {e}")

all_target_urls = list(set(all_target_urls))[:MAX_FILES]

if not all_target_urls:
    print("[FATAL] Tidak ada URL XRootD yang valid terdeteksi. Eksekusi dihentikan.")
    exit()

print(f"\n[STATUS] TOTAL AMUNISI TERKUNCI: {len(all_target_urls)} file XRootD\n", flush=True)

# ------------------------------------------------------------------------------
# 2. MESIN EKSTRAKSI ASINKRON (FULL KINEMATIC M_T)
# ------------------------------------------------------------------------------
def stream_and_analyze(url):
    local_results = {
        "di_photon_mass": [], "four_lepton_mass": [], "mt_distribution": [], 
        "events": 0, "success": False
    }
    
    for attempt in range(3):
        try:
            with uproot.open(f"{url}:Events", timeout=180) as events:
                chunk = events.arrays(["nPhoton", "Photon_pt", "Photon_eta", "Photon_phi",
                                       "nElectron", "Electron_pt", "Electron_eta", "Electron_phi", 
                                       "nMuon", "Muon_pt", "Muon_eta", "Muon_phi", 
                                       "MET_pt", "MET_phi"], 
                                       how=dict, library="ak")
                
                local_results["events"] = len(chunk["MET_pt"])
                
                # --- CHANNEL 1: DI-PHOTON ---
                if "nPhoton" in chunk.keys():
                    mask_2ph = chunk["nPhoton"] >= 2
                    if ak.any(mask_2ph):
                        g1 = vector.zip({"pt": chunk["Photon_pt"][mask_2ph, 0], "eta": chunk["Photon_eta"][mask_2ph, 0], "phi": chunk["Photon_phi"][mask_2ph, 0], "mass": 0})
                        g2 = vector.zip({"pt": chunk["Photon_pt"][mask_2ph, 1], "eta": chunk["Photon_eta"][mask_2ph, 1], "phi": chunk["Photon_phi"][mask_2ph, 1], "mass": 0})
                        ph_cut = (g1.pt > 40) & (g2.pt > 40) & (abs(g1.eta) < 2.5) & (abs(g2.eta) < 2.5) & (abs(g1.deltaphi(g2)) > 2.5)
                        ph_mass = (g1[ph_cut] + g2[ph_cut]).mass
                        local_results["di_photon_mass"].append(ak.to_numpy(ph_mass[ph_mass > 150])) 

                # --- CHANNEL 2: FOUR-LEPTON ---
                if "nElectron" in chunk.keys() and "nMuon" in chunk.keys():
                    tot_lep = chunk["nElectron"] + chunk["nMuon"]
                    mask_4l = tot_lep >= 4
                    if ak.any(mask_4l):
                        total_4l_pt = ak.sum(chunk["Electron_pt"][mask_4l], axis=-1) + ak.sum(chunk["Muon_pt"][mask_4l], axis=-1)
                        local_results["four_lepton_mass"].append(ak.to_numpy(total_4l_pt[total_4l_pt > 150]))

                # --- CHANNEL 3: FULL KINEMATIC TRANSVERSE MASS (M_T) ---
                if "MET_pt" in chunk.keys() and "MET_phi" in chunk.keys() and "nElectron" in chunk.keys() and "nMuon" in chunk.keys():
                    
                    # 1. Tameng Energi & Geometri: p_T > 30 GeV DAN Abs(eta) < 2.4
                    e_mask = (chunk["Electron_pt"] > 30.0) & (abs(chunk["Electron_eta"]) < 2.4)
                    mu_mask = (chunk["Muon_pt"] > 30.0) & (abs(chunk["Muon_eta"]) < 2.4)
                    
                    n_valid_e = ak.sum(e_mask, axis=-1)
                    n_valid_mu = ak.sum(mu_mask, axis=-1)
                    
                    # 2. Syarat Topologi: Tepat 2 lepton valid
                    mask_2l_strict = (n_valid_e + n_valid_mu) == 2
                    
                    if ak.any(mask_2l_strict):
                        chunk_pass = {k: v[mask_2l_strict] for k, v in chunk.items()}
                        e_pass = e_mask[mask_2l_strict]
                        mu_pass = mu_mask[mask_2l_strict]
                        
                        # Hitung Vektor Kartesian Transversal (px, py)
                        px_e = chunk_pass["Electron_pt"][e_pass] * np.cos(chunk_pass["Electron_phi"][e_pass])
                        py_e = chunk_pass["Electron_pt"][e_pass] * np.sin(chunk_pass["Electron_phi"][e_pass])
                        px_mu = chunk_pass["Muon_pt"][mu_pass] * np.cos(chunk_pass["Muon_phi"][mu_pass])
                        py_mu = chunk_pass["Muon_pt"][mu_pass] * np.sin(chunk_pass["Muon_phi"][mu_pass])
                        
                        # Gabungkan vektor untuk menjadi Sistem Dilepton
                        sum_px_ll = ak.sum(px_e, axis=-1) + ak.sum(px_mu, axis=-1)
                        sum_py_ll = ak.sum(py_e, axis=-1) + ak.sum(py_mu, axis=-1)
                        pt_ll = np.sqrt(sum_px_ll**2 + sum_py_ll**2)
                        
                        # Hitung Energi (E) dan Momentum Longitudinal (pz) asumsikan lepton massless
                        e_E = chunk_pass["Electron_pt"][e_pass] * np.cosh(chunk_pass["Electron_eta"][e_pass])
                        mu_E = chunk_pass["Muon_pt"][mu_pass] * np.cosh(chunk_pass["Muon_eta"][mu_pass])
                        pz_e = chunk_pass["Electron_pt"][e_pass] * np.sinh(chunk_pass["Electron_eta"][e_pass])
                        pz_mu = chunk_pass["Muon_pt"][mu_pass] * np.sinh(chunk_pass["Muon_eta"][mu_pass])
                        
                        sum_E_ll = ak.sum(e_E, axis=-1) + ak.sum(mu_E, axis=-1)
                        sum_pz_ll = ak.sum(pz_e, axis=-1) + ak.sum(pz_mu, axis=-1)
                        
                        # Massa Invarian Sistem Dilepton Kuadrat
                        m_ll_sq = sum_E_ll**2 - (sum_px_ll**2 + sum_py_ll**2 + sum_pz_ll**2)
                        m_ll_sq = ak.where(m_ll_sq < 0, 0, m_ll_sq) # Mencegah error presisi float
                        
                        # Energi Transversal Sistem Dilepton
                        ET_ll = np.sqrt(pt_ll**2 + m_ll_sq)
                        
                        # Vektor MET
                        met_pt = chunk_pass["MET_pt"]
                        met_phi = chunk_pass["MET_phi"]
                        met_px = met_pt * np.cos(met_phi)
                        met_py = met_pt * np.sin(met_phi)
                        
                        # KALKULASI TRANSVERSE MASS (M_T) ABSOLUT
                        mt_sq = (ET_ll + met_pt)**2 - (sum_px_ll + met_px)**2 - (sum_py_ll + met_py)**2
                        mt_sq = ak.where(mt_sq < 0, 0, mt_sq) # Mencegah NaN
                        mt_data = np.sqrt(mt_sq)
                        
                        # Simpan hasil perhitungan M_T penuh
                        local_results["mt_distribution"].append(ak.to_numpy(mt_data[mt_data > 100])) 
            
            local_results["success"] = True
            break 
            
        except Exception as e:
            if attempt == 2:
                print(f"[ERROR] Gagal membaca {url.split('/')[-1]} | Tipe: {type(e).__name__}", flush=True)
            time.sleep(3)
            
    return local_results

start_time = time.time()
global_results = {"di_photon_mass": [], "four_lepton_mass": [], "mt_distribution": []}
total_events = 0
processed_files = 0

print(f"[PROCESS] Memulai eksekusi threading ({MAX_WORKERS} Workers)...", flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    futures = {executor.submit(stream_and_analyze, url): url for url in all_target_urls}
    
    for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
        res = future.result()
        
        if res["success"]:
            total_events += res["events"]
            processed_files += 1
            if res["di_photon_mass"]: global_results["di_photon_mass"].extend(res["di_photon_mass"])
            if res["four_lepton_mass"]: global_results["four_lepton_mass"].extend(res["four_lepton_mass"])
            if res["mt_distribution"]: global_results["mt_distribution"].extend(res["mt_distribution"])
                
        if i % 50 == 0 or i == len(all_target_urls):
            print(f"[STATUS] Diproses: {i}/{len(all_target_urls)} | Sukses: {processed_files} | Akumulasi Events: {total_events:,}", flush=True)

exec_time = time.time() - start_time
print(f"\n[INFO] OPERASI SELESAI. Ekstraksi Berhasil: {processed_files}/{len(all_target_urls)}", flush=True)
print(f"[INFO] Waktu Eksekusi: {exec_time/3600:.2f} Jam", flush=True)

# ------------------------------------------------------------------------------
# 3. KONSOLIDASI DATA DAN VISUALISASI
# ------------------------------------------------------------------------------
print("[PROCESS] Mengkompilasi matriks diagnostik...", flush=True)
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

golden_mt = []
if global_results["mt_distribution"]:
    mt_all = np.concatenate(global_results["mt_distribution"])
    mt_plot = mt_all[(mt_all > 150) & (mt_all < 250)]
    counts_mt, bins_mt = np.histogram(mt_plot, bins=100, range=(150, 250))
    axs[2].bar((bins_mt[:-1]+bins_mt[1:])/2, counts_mt, width=1.0, color='orange', edgecolor='black', alpha=0.9)
    
    golden_mask = (mt_all >= 215.0) & (mt_all <= 216.0)
    golden_mt = mt_all[golden_mask]
    
    # Pengamanan absolut array Full M_T
    np.save("mt_all_distribution.npy", mt_all)

axs[2].axvline(x=TBP_THRESHOLD, color='red', linestyle='dashed', linewidth=3)
axs[2].set_title(r'Channel 3: Full Dilepton+MET Transverse Mass ($M_T$)', fontweight='bold', color='white')
axs[2].set_facecolor('black')
axs[2].set_xlim(180, 250)

for ax in axs: ax.tick_params(colors='white'); ax.spines['bottom'].set_color('white'); ax.spines['left'].set_color('white')

fig.patch.set_facecolor('#121212')
plt.savefig("TBP_Final_Discovery_V6.png", dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')

if len(golden_mt) > 0:
    with open("TBP_Discovery_Anomalies_V6.csv", "w") as f:
        f.write("Index,Full_MT_GeV\n")
        for i, val in enumerate(golden_mt): f.write(f"Anomali_{i+1},{val:.6f}\n")

print(f"[SUCCESS] {len(golden_mt)} anomali disimpan dan matriks mt_all_distribution.npy diamankan.", flush=True)
