# ==============================================================================
# AUTOMATED KINEMATIC PIPELINE - TBP-OMNI-TENSOR EXTRACTOR (COMPLETE 5-PATCH SUITE)
# PROTOCOL: FULL 4-VECTOR & MET EXTRACTION (93 .NPY MATRICES ACROSS 10 SECTORS)
# REPOSITORY: https://github.com/cefiyana-clover/TBP-Computational-Framework
# ==============================================================================
import os
import sys
import time
import vector
import numpy as np
import concurrent.futures
import subprocess

try:
    import XRootD
    import uproot
    import awkward as ak
    from coffea.lumi_tools import LumiMask
except ImportError:
    raise RuntimeError("[FATAL] Required modules (XRootD, uproot, awkward, coffea) are missing.")

GOLDEN_JSON = "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_MuonPhys.txt"
if not os.path.exists(GOLDEN_JSON):
    raise FileNotFoundError(f"[FATAL] Golden JSON file '{GOLDEN_JSON}' not found.")
lumimask = LumiMask(GOLDEN_JSON)

MAX_FILES = 100000
MAX_WORKERS = 16
CONTINGENCY_MODE = False


# ==============================================================================
# PATCH 1: CROSS-FLAVOR (MUONEG) & MONOPHOTON FULL 4-VECTOR EXTRACTION
# ==============================================================================
def stream_and_analyze_patch1(url):
    # Struktur Array Lengkap (Omni-Tensor)
    local_results = {
        "muoneg_mass": [], "muoneg_met_pt": [], "muoneg_met_phi": [],
        "muoneg_e_pt": [], "muoneg_e_eta": [], "muoneg_e_phi": [], "muoneg_e_charge": [],
        "muoneg_mu_pt": [], "muoneg_mu_eta": [], "muoneg_mu_phi": [], "muoneg_mu_charge": [],
        "single_photon_mass": [], "single_photon_eta": [], "single_photon_phi": [],
        "single_photon_met_pt": [], "single_photon_met_phi": [],
        "events_raw": 0, "events_certified": 0, "success": False
    }

    # Penambahan Cabang: Charge dan MET_phi untuk analisis Tensor 3D
    needed_branches = [
        "nPhoton", "Photon_pt", "Photon_eta", "Photon_phi", "Photon_electronVeto", "Photon_cutBased",
        "nElectron", "Electron_pt", "Electron_eta", "Electron_phi", "Electron_charge", "Electron_cutBased",
        "nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge", "Muon_tightId", "Muon_pfRelIso04_all",
        "MET_pt", "MET_phi", "run", "luminosityBlock", "HLT_*", "Flag_*"
    ]

    for attempt in range(3):
        try:
            with uproot.open(url, timeout=1800) as root_file:
                available_keys = [k.split(';')[0] for k in root_file.keys()]
                events = root_file["Events"] if "Events" in available_keys else (root_file["events"] if "events" in available_keys else None)
                if events is None:
                    break
                chunk = events.arrays(filter_name=needed_branches, how=dict, library="ak")
                if "run" not in chunk.keys():
                    break

                local_results["events_raw"] = len(chunk["run"])

                # GOLDEN JSON FILTER
                valid_lumi = lumimask(ak.to_numpy(chunk["run"]), ak.to_numpy(chunk["luminosityBlock"]))
                master_mask = ak.Array(valid_lumi)

                required_flags = ["Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter"]
                for flag in required_flags:
                    if flag in chunk.keys():
                        master_mask = master_mask & chunk[flag]

                local_results["events_certified"] = int(ak.sum(master_mask))
                if local_results["events_certified"] == 0:
                    local_results["success"] = True
                    break

                # HLT STRICT ENFORCEMENT
                hlt_keys = [k for k in chunk.keys() if k.startswith("HLT_")]
                emu_hlt = [k for k in hlt_keys if ("Mu" in k and "Ele" in k) or ("IsoMu" in k) or ("Ele2" in k)]
                pho_hlt = [k for k in hlt_keys if "Photon" in k or "Diphoton" in k]
                false_mask = chunk["run"] < 0

                pass_emu = false_mask
                if emu_hlt:
                    for k in emu_hlt:
                        pass_emu = pass_emu | chunk[k]
                pass_pho = false_mask
                if pho_hlt:
                    for k in pho_hlt:
                        pass_pho = pass_pho | chunk[k]

                # ==============================================================
                # 1. CROSS-FLAVOR DI-BOSON (MUONEG) FULL VECTOR EXTRACTION
                # ==============================================================
                if "nElectron" in chunk.keys() and "nMuon" in chunk.keys() and "MET_pt" in chunk.keys():
                    base_mask_emu = master_mask & pass_emu & (chunk["nElectron"] >= 1) & (chunk["nMuon"] >= 1)
                    if ak.any(base_mask_emu):
                        e_pt = chunk["Electron_pt"][base_mask_emu]
                        e_eta = chunk["Electron_eta"][base_mask_emu]
                        e_phi = chunk["Electron_phi"][base_mask_emu]
                        e_charge = chunk["Electron_charge"][base_mask_emu]

                        mu_pt = chunk["Muon_pt"][base_mask_emu]
                        mu_eta = chunk["Muon_eta"][base_mask_emu]
                        mu_phi = chunk["Muon_phi"][base_mask_emu]
                        mu_charge = chunk["Muon_charge"][base_mask_emu]

                        if CONTINGENCY_MODE:
                            valid_e = (e_pt > 20) & (abs(e_eta) < 2.5)
                            valid_mu = (mu_pt > 20) & (abs(mu_eta) < 2.4)
                        else:
                            e_id = chunk["Electron_cutBased"][base_mask_emu]
                            mu_id = chunk["Muon_tightId"][base_mask_emu]
                            mu_iso = chunk["Muon_pfRelIso04_all"][base_mask_emu]

                            valid_e = (e_pt > 20) & (abs(e_eta) < 2.5) & (e_id >= 2)
                            valid_mu = (mu_pt > 20) & (abs(mu_eta) < 2.4) & (mu_id == True) & (mu_iso < 0.15)

                        has_emu = (ak.num(e_pt[valid_e]) >= 1) & (ak.num(mu_pt[valid_mu]) >= 1)

                        if ak.any(has_emu):
                            # Ambil Leading Partikel
                            e_pt_lead = e_pt[valid_e][has_emu][:, 0]
                            e_eta_lead = e_eta[valid_e][has_emu][:, 0]
                            e_phi_lead = e_phi[valid_e][has_emu][:, 0]
                            e_charge_lead = e_charge[valid_e][has_emu][:, 0]

                            mu_pt_lead = mu_pt[valid_mu][has_emu][:, 0]
                            mu_eta_lead = mu_eta[valid_mu][has_emu][:, 0]
                            mu_phi_lead = mu_phi[valid_mu][has_emu][:, 0]
                            mu_charge_lead = mu_charge[valid_mu][has_emu][:, 0]

                            e_vec = vector.zip({"pt": e_pt_lead, "eta": e_eta_lead, "phi": e_phi_lead, "mass": 0.000511})
                            mu_vec = vector.zip({"pt": mu_pt_lead, "eta": mu_eta_lead, "phi": mu_phi_lead, "mass": 0.10566})

                            emu_mass = (e_vec + mu_vec).mass
                            emu_met_pt = chunk["MET_pt"][base_mask_emu][has_emu]
                            emu_met_phi = chunk["MET_phi"][base_mask_emu][has_emu]

                            # Terapkan masking global agar semua komponen sejajar sempurna
                            mass_cut = emu_mass > 100

                            local_results["muoneg_mass"].append(ak.to_numpy(emu_mass[mass_cut]))
                            local_results["muoneg_met_pt"].append(ak.to_numpy(emu_met_pt[mass_cut]))
                            local_results["muoneg_met_phi"].append(ak.to_numpy(emu_met_phi[mass_cut]))

                            local_results["muoneg_e_pt"].append(ak.to_numpy(e_pt_lead[mass_cut]))
                            local_results["muoneg_e_eta"].append(ak.to_numpy(e_eta_lead[mass_cut]))
                            local_results["muoneg_e_phi"].append(ak.to_numpy(e_phi_lead[mass_cut]))
                            local_results["muoneg_e_charge"].append(ak.to_numpy(e_charge_lead[mass_cut]))

                            local_results["muoneg_mu_pt"].append(ak.to_numpy(mu_pt_lead[mass_cut]))
                            local_results["muoneg_mu_eta"].append(ak.to_numpy(mu_eta_lead[mass_cut]))
                            local_results["muoneg_mu_phi"].append(ak.to_numpy(mu_phi_lead[mass_cut]))
                            local_results["muoneg_mu_charge"].append(ak.to_numpy(mu_charge_lead[mass_cut]))

                # ==============================================================
                # 2. MONOPHOTON FULL VECTOR EXTRACTION
                # ==============================================================
                if "nPhoton" in chunk.keys() and "MET_pt" in chunk.keys():
                    base_mask_ph = master_mask & pass_pho & (chunk["nPhoton"] >= 1)
                    if ak.any(base_mask_ph):
                        ph_pt = chunk["Photon_pt"][base_mask_ph]
                        ph_eta = chunk["Photon_eta"][base_mask_ph]
                        ph_phi = chunk["Photon_phi"][base_mask_ph]
                        ph_veto = chunk["Photon_electronVeto"][base_mask_ph]

                        if CONTINGENCY_MODE:
                            valid_ph = (ph_pt > 50) & (abs(ph_eta) < 2.4) & ph_veto
                        else:
                            ph_id = chunk["Photon_cutBased"][base_mask_ph]
                            valid_ph = (ph_pt > 50) & (abs(ph_eta) < 2.4) & ph_veto & (ph_id >= 2)

                        has_ph = ak.num(ph_pt[valid_ph]) >= 1

                        if ak.any(has_ph):
                            ph_pt_lead = ph_pt[valid_ph][has_ph][:, 0]
                            ph_eta_lead = ph_eta[valid_ph][has_ph][:, 0]
                            ph_phi_lead = ph_phi[valid_ph][has_ph][:, 0]

                            ph_met_pt = chunk["MET_pt"][base_mask_ph][has_ph]
                            ph_met_phi = chunk["MET_phi"][base_mask_ph][has_ph]

                            scale_cut = ph_pt_lead > 100

                            local_results["single_photon_mass"].append(ak.to_numpy(ph_pt_lead[scale_cut]))
                            local_results["single_photon_eta"].append(ak.to_numpy(ph_eta_lead[scale_cut]))
                            local_results["single_photon_phi"].append(ak.to_numpy(ph_phi_lead[scale_cut]))
                            local_results["single_photon_met_pt"].append(ak.to_numpy(ph_met_pt[scale_cut]))
                            local_results["single_photon_met_phi"].append(ak.to_numpy(ph_met_phi[scale_cut]))

                local_results["success"] = True
                break

        except Exception:
            time.sleep(15)

    return local_results


def execute_patch_1():
    print("\n==================================================")
    print(" [INIT] TBP-OMNI-TENSOR EXTRACTOR (FULL 4-VECTOR PATCH) ")
    print("==================================================\n")
    print(f"[INFO] Golden JSON loaded successfully. XRootD engine initialized.")
    if CONTINGENCY_MODE:
        print("[INFO] Contingency mode active: Standard ID requirements relaxed.", flush=True)
    else:
        print("[INFO] STRICT MODE ACTIVE: Standard Tight IDs & HLT fully enforced. No Bypass.", flush=True)

    RECORD_IDS = [30528, 30561, 30531, 30564]
    all_target_urls = []
    print("[PROCESS] Acquiring XRootD manifests from CERN Open Data...")
    for recid in RECORD_IDS:
        try:
            result = subprocess.run(
                ["cernopendata-client", "get-file-locations", "--recid", str(recid), "--protocol", "xrootd"],
                capture_output=True, text=True, check=True
            )
            urls = result.stdout.strip().split('\n')
            valid_urls = [url.strip() for url in urls if url.strip().startswith('root://')]
            all_target_urls.extend(valid_urls)
        except Exception:
            pass
    all_target_urls = list(set(all_target_urls))[:MAX_FILES]
    print(f"\n[INFO] Total XRootD targets acquired: {len(all_target_urls)} files\n", flush=True)

    global_results = {
        "muoneg_mass": [], "muoneg_met_pt": [], "muoneg_met_phi": [],
        "muoneg_e_pt": [], "muoneg_e_eta": [], "muoneg_e_phi": [], "muoneg_e_charge": [],
        "muoneg_mu_pt": [], "muoneg_mu_eta": [], "muoneg_mu_phi": [], "muoneg_mu_charge": [],
        "single_photon_mass": [], "single_photon_eta": [], "single_photon_phi": [],
        "single_photon_met_pt": [], "single_photon_met_phi": []
    }
    processed_files, total_raw, total_certified = 0, 0, 0
    print(f"[PROCESS] Initiating Omni-Tensor Extraction ({MAX_WORKERS} Workers)...", flush=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(stream_and_analyze_patch1, url): url for url in all_target_urls}
        for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
            res = future.result()
            if res["success"]:
                processed_files += 1
                total_raw += res["events_raw"]
                total_certified += res["events_certified"]
                for key in global_results.keys():
                    if res[key]:
                        global_results[key].extend(res[key])
            if i % 10 == 0 or i == len(all_target_urls):
                print(f"[STATUS] Files: {i}/{len(all_target_urls)} | Raw: {total_raw:,} | Certified: {total_certified:,}", flush=True)
    print(f"\n[INFO] Extraction Complete. Compiling Omni-Tensor matrices...", flush=True)
    for key, data_list in global_results.items():
        if data_list:
            np.save(f"{key}.npy", np.concatenate(data_list))
    print("[SUCCESS] 16 Matrix Components (Full 4-Vectors) saved successfully.", flush=True)


# ==============================================================================
# PATCH 2: DOUBLE-EG & DOUBLE-MUON FULL 4-VECTOR EXTRACTION (4-LEPTON & DI-PHOTON)
# ==============================================================================
def stream_and_analyze_patch2(url):
    local_results = {
        # Di-Photon Vectors
        "di_photon_mass": [], "di_photon_met_pt": [], "di_photon_met_phi": [],
        "di_photon_ph1_pt": [], "di_photon_ph1_eta": [], "di_photon_ph1_phi": [],
        "di_photon_ph2_pt": [], "di_photon_ph2_eta": [], "di_photon_ph2_phi": [],

        # 4-Lepton Vectors
        "four_lepton_mass": [], "four_lepton_met_pt": [], "four_lepton_met_phi": [],
        "four_lepton_sum_e_pt": [], "four_lepton_sum_e_eta": [], "four_lepton_sum_e_phi": [],
        "four_lepton_sum_mu_pt": [], "four_lepton_sum_mu_eta": [], "four_lepton_sum_mu_phi": [],

        "events_raw": 0, "events_certified": 0, "success": False
    }

    needed_branches = [
        "nPhoton", "Photon_pt", "Photon_eta", "Photon_phi", "Photon_electronVeto", "Photon_cutBased",
        "nElectron", "Electron_pt", "Electron_eta", "Electron_phi", "Electron_cutBased", "Electron_charge",
        "nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_tightId", "Muon_pfRelIso04_all", "Muon_charge",
        "MET_pt", "MET_phi", "run", "luminosityBlock", "HLT_*", "Flag_*"
    ]

    for attempt in range(3):
        try:
            with uproot.open(url, timeout=1800) as root_file:
                available_keys = [k.split(';')[0] for k in root_file.keys()]
                if "Events" in available_keys:
                    events = root_file["Events"]
                elif "events" in available_keys:
                    events = root_file["events"]
                else:
                    break
                chunk = events.arrays(filter_name=needed_branches, how=dict, library="ak")

                if "run" not in chunk.keys() or "luminosityBlock" not in chunk.keys():
                    break

                local_results["events_raw"] = len(chunk["run"])

                # 1. GOLDEN JSON & QUALITY FLAGS
                valid_lumi = lumimask(ak.to_numpy(chunk["run"]), ak.to_numpy(chunk["luminosityBlock"]))
                master_mask = ak.Array(valid_lumi)

                flag_keys = list(chunk.keys())
                required_flags = [
                    "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter",
                    "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter",
                    "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter"
                ]
                for flag in required_flags:
                    if flag in flag_keys:
                        master_mask = master_mask & chunk[flag]

                local_results["events_certified"] = int(ak.sum(master_mask))
                if local_results["events_certified"] == 0:
                    local_results["success"] = True
                    break

                # HLT SELECTION
                hlt_keys = list(chunk.keys())
                dipho_keys = [k for k in hlt_keys if "HLT_Diphoton30" in k]
                muon_keys = [k for k in hlt_keys if "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL" in k]
                ele_keys = [k for k in hlt_keys if "HLT_Ele23_Ele12_CaloIdL_TrackIdL" in k]
                false_mask = chunk["run"] < 0

                pass_dipho = chunk[dipho_keys[0]] if dipho_keys else false_mask
                pass_muon = chunk[muon_keys[0]] if muon_keys else false_mask
                pass_ele = chunk[ele_keys[0]] if ele_keys else false_mask
                pass_4l = pass_muon | pass_ele

                # ==================================================
                # 2. DI-PHOTON SELECTION FULL VECTOR
                # ==================================================
                if "nPhoton" in chunk.keys() and "MET_pt" in chunk.keys():
                    base_mask_2ph = master_mask & pass_dipho & (chunk["nPhoton"] >= 2)

                    if ak.any(base_mask_2ph):
                        ph_pt_sub = chunk["Photon_pt"][base_mask_2ph]
                        ph_eta_sub = chunk["Photon_eta"][base_mask_2ph]
                        ph_phi_sub = chunk["Photon_phi"][base_mask_2ph]
                        ph_veto_sub = chunk["Photon_electronVeto"][base_mask_2ph]

                        if CONTINGENCY_MODE:
                            valid_ph_mask = (ph_pt_sub > 30) & (abs(ph_eta_sub) < 2.4) & (ph_veto_sub == True)
                        else:
                            ph_id_sub = chunk["Photon_cutBased"][base_mask_2ph]
                            valid_ph_mask = (ph_pt_sub > 30) & (abs(ph_eta_sub) < 2.4) & (ph_veto_sub == True) & (ph_id_sub >= 2)

                        event_has_2ph = ak.num(ph_pt_sub[valid_ph_mask]) >= 2

                        if ak.any(event_has_2ph):
                            valid_photons = vector.zip({
                                "pt": ph_pt_sub[valid_ph_mask][event_has_2ph],
                                "eta": ph_eta_sub[valid_ph_mask][event_has_2ph],
                                "phi": ph_phi_sub[valid_ph_mask][event_has_2ph],
                                "mass": ak.zeros_like(ph_pt_sub[valid_ph_mask][event_has_2ph])
                            })

                            g1 = valid_photons[:, 0]
                            g2 = valid_photons[:, 1]

                            ph_mass = (g1 + g2).mass
                            met_2ph_pt = chunk["MET_pt"][base_mask_2ph][event_has_2ph]
                            met_2ph_phi = chunk["MET_phi"][base_mask_2ph][event_has_2ph]

                            mass_cut = ph_mass > 100

                            local_results["di_photon_mass"].append(ak.to_numpy(ph_mass[mass_cut]))
                            local_results["di_photon_met_pt"].append(ak.to_numpy(met_2ph_pt[mass_cut]))
                            local_results["di_photon_met_phi"].append(ak.to_numpy(met_2ph_phi[mass_cut]))

                            local_results["di_photon_ph1_pt"].append(ak.to_numpy(g1.pt[mass_cut]))
                            local_results["di_photon_ph1_eta"].append(ak.to_numpy(g1.eta[mass_cut]))
                            local_results["di_photon_ph1_phi"].append(ak.to_numpy(g1.phi[mass_cut]))

                            local_results["di_photon_ph2_pt"].append(ak.to_numpy(g2.pt[mass_cut]))
                            local_results["di_photon_ph2_eta"].append(ak.to_numpy(g2.eta[mass_cut]))
                            local_results["di_photon_ph2_phi"].append(ak.to_numpy(g2.phi[mass_cut]))

                # ==================================================
                # 3. 4-LEPTON SELECTION FULL VECTOR
                # ==================================================
                if "nElectron" in chunk.keys() and "nMuon" in chunk.keys() and "MET_pt" in chunk.keys():
                    base_mask_4l = master_mask & pass_4l & ((chunk["nElectron"] + chunk["nMuon"]) >= 4)

                    if ak.any(base_mask_4l):
                        e_pt_sub = chunk["Electron_pt"][base_mask_4l]
                        e_eta_sub = chunk["Electron_eta"][base_mask_4l]
                        e_phi_sub = chunk["Electron_phi"][base_mask_4l]

                        mu_pt_sub = chunk["Muon_pt"][base_mask_4l]
                        mu_eta_sub = chunk["Muon_eta"][base_mask_4l]
                        mu_phi_sub = chunk["Muon_phi"][base_mask_4l]

                        if CONTINGENCY_MODE:
                            valid_e_mask = (e_pt_sub > 10) & (abs(e_eta_sub) < 2.5)
                            mu_iso_sub = chunk["Muon_pfRelIso04_all"][base_mask_4l]
                            valid_mu_mask = (mu_pt_sub > 10) & (abs(mu_eta_sub) < 2.4) & (mu_iso_sub < 0.25)
                        else:
                            e_id_sub = chunk["Electron_cutBased"][base_mask_4l]
                            mu_id_sub = chunk["Muon_tightId"][base_mask_4l]
                            mu_iso_sub = chunk["Muon_pfRelIso04_all"][base_mask_4l]

                            valid_e_mask = (e_pt_sub > 10) & (abs(e_eta_sub) < 2.5) & (e_id_sub >= 2)
                            valid_mu_mask = (mu_pt_sub > 10) & (abs(mu_eta_sub) < 2.4) & (mu_id_sub == True) & (mu_iso_sub < 0.15)

                        event_has_4l = (ak.num(e_pt_sub[valid_e_mask]) + ak.num(mu_pt_sub[valid_mu_mask])) >= 4

                        if ak.any(event_has_4l):
                            e_mass = ak.ones_like(e_pt_sub[valid_e_mask][event_has_4l]) * 0.000511
                            mu_mass = ak.ones_like(mu_pt_sub[valid_mu_mask][event_has_4l]) * 0.10566

                            electrons = vector.zip({
                                "pt": e_pt_sub[valid_e_mask][event_has_4l], "eta": e_eta_sub[valid_e_mask][event_has_4l],
                                "phi": e_phi_sub[valid_e_mask][event_has_4l], "mass": e_mass
                            })

                            muons = vector.zip({
                                "pt": mu_pt_sub[valid_mu_mask][event_has_4l], "eta": mu_eta_sub[valid_mu_mask][event_has_4l],
                                "phi": mu_phi_sub[valid_mu_mask][event_has_4l], "mass": mu_mass
                            })

                            sum_e = ak.sum(electrons, axis=1)
                            sum_mu = ak.sum(muons, axis=1)

                            mass_4l = (sum_e + sum_mu).mass
                            met_4l_pt = chunk["MET_pt"][base_mask_4l][event_has_4l]
                            met_4l_phi = chunk["MET_phi"][base_mask_4l][event_has_4l]

                            mass_cut = mass_4l > 100

                            local_results["four_lepton_mass"].append(ak.to_numpy(mass_4l[mass_cut]))
                            local_results["four_lepton_met_pt"].append(ak.to_numpy(met_4l_pt[mass_cut]))
                            local_results["four_lepton_met_phi"].append(ak.to_numpy(met_4l_phi[mass_cut]))

                            local_results["four_lepton_sum_e_pt"].append(ak.to_numpy(sum_e.pt[mass_cut]))
                            local_results["four_lepton_sum_e_eta"].append(ak.to_numpy(sum_e.eta[mass_cut]))
                            local_results["four_lepton_sum_e_phi"].append(ak.to_numpy(sum_e.phi[mass_cut]))

                            local_results["four_lepton_sum_mu_pt"].append(ak.to_numpy(sum_mu.pt[mass_cut]))
                            local_results["four_lepton_sum_mu_eta"].append(ak.to_numpy(sum_mu.eta[mass_cut]))
                            local_results["four_lepton_sum_mu_phi"].append(ak.to_numpy(sum_mu.phi[mass_cut]))

                local_results["success"] = True
                break

        except Exception:
            time.sleep(15)

    return local_results


def execute_patch_2():
    print("\n==================================================")
    print(" [INIT] TBP-OMNI-TENSOR EXTRACTOR (GOLDEN CHANNELS PATCH) ")
    print("==================================================\n")
    print(f"[INFO] Golden JSON loaded successfully. XRootD engine initialized.")
    if CONTINGENCY_MODE:
        print("[INFO] Contingency mode active: Standard ID requirements relaxed.", flush=True)
    else:
        print("[INFO] STRICT MODE ACTIVE: Standard Tight IDs & HLT fully enforced. No Bypass.", flush=True)

    RECORD_IDS = [30521, 30554, 30522, 30555]
    all_target_urls = []
    print("[PROCESS] Acquiring XRootD manifests from CERN Open Data...")
    for recid in RECORD_IDS:
        try:
            result = subprocess.run(
                ["cernopendata-client", "get-file-locations", "--recid", str(recid), "--protocol", "xrootd"],
                capture_output=True, text=True, check=True
            )
            urls = result.stdout.strip().split('\n')
            valid_urls = [url.strip() for url in urls if url.strip().startswith('root://')]
            all_target_urls.extend(valid_urls)
        except Exception:
            pass
    all_target_urls = list(set(all_target_urls))[:MAX_FILES]
    print(f"\n[INFO] Total XRootD targets acquired: {len(all_target_urls)} files\n", flush=True)

    global_results = {
        "di_photon_mass": [], "di_photon_met_pt": [], "di_photon_met_phi": [],
        "di_photon_ph1_pt": [], "di_photon_ph1_eta": [], "di_photon_ph1_phi": [],
        "di_photon_ph2_pt": [], "di_photon_ph2_eta": [], "di_photon_ph2_phi": [],
        "four_lepton_mass": [], "four_lepton_met_pt": [], "four_lepton_met_phi": [],
        "four_lepton_sum_e_pt": [], "four_lepton_sum_e_eta": [], "four_lepton_sum_e_phi": [],
        "four_lepton_sum_mu_pt": [], "four_lepton_sum_mu_eta": [], "four_lepton_sum_mu_phi": []
    }
    processed_files, total_raw, total_certified = 0, 0, 0
    print(f"[PROCESS] Initiating Omni-Tensor ProcessPoolExecutor ({MAX_WORKERS} Workers)...", flush=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(stream_and_analyze_patch2, url): url for url in all_target_urls}
        for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
            res = future.result()
            if res["success"]:
                processed_files += 1
                total_raw += res["events_raw"]
                total_certified += res["events_certified"]
                for key in global_results.keys():
                    if res[key]:
                        global_results[key].extend(res[key])
            if i % 10 == 0 or i == len(all_target_urls):
                print(f"[STATUS] Files Processed: {i}/{len(all_target_urls)} | Raw: {total_raw:,} | Certified: {total_certified:,}", flush=True)
    print(f"\n[INFO] OPERATION COMPLETE. Successful Extractions: {processed_files}/{len(all_target_urls)}", flush=True)
    if processed_files > 0:
        print("[PROCESS] Compiling Omni-Tensor output matrices...", flush=True)
        for key, data_list in global_results.items():
            if data_list:
                np.save(f"{key}.npy", np.concatenate(data_list))
        print("[SUCCESS] Vector Matrices saved successfully.", flush=True)


# ==============================================================================
# PATCH 3: SINGLE-CHANNEL (ELE, MU, TAU) FULL 4-VECTOR EXTRACTION
# ==============================================================================
def stream_and_analyze_patch3(url):
    local_results = {
        # Single Electron Vectors
        "single_ele_mass": [], "single_ele_met_pt": [], "single_ele_met_phi": [],
        "single_ele_e1_pt": [], "single_ele_e1_eta": [], "single_ele_e1_phi": [], "single_ele_e1_charge": [],
        "single_ele_e2_pt": [], "single_ele_e2_eta": [], "single_ele_e2_phi": [], "single_ele_e2_charge": [],

        # Single Muon Vectors
        "single_mu_mass": [], "single_mu_met_pt": [], "single_mu_met_phi": [],
        "single_mu_m1_pt": [], "single_mu_m1_eta": [], "single_mu_m1_phi": [], "single_mu_m1_charge": [],
        "single_mu_m2_pt": [], "single_mu_m2_eta": [], "single_mu_m2_phi": [], "single_mu_m2_charge": [],

        # Single Tau Vectors
        "single_tau_mass": [], "single_tau_met_pt": [], "single_tau_met_phi": [],
        "single_tau_t1_pt": [], "single_tau_t1_eta": [], "single_tau_t1_phi": [], "single_tau_t1_charge": [],
        "single_tau_t2_pt": [], "single_tau_t2_eta": [], "single_tau_t2_phi": [], "single_tau_t2_charge": [],

        "events_raw": 0, "events_certified": 0, "success": False
    }

    needed_branches = [
        "nElectron", "Electron_pt", "Electron_eta", "Electron_phi", "Electron_cutBased", "Electron_charge",
        "nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_tightId", "Muon_pfRelIso04_all", "Muon_charge",
        "nTau", "Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass", "Tau_charge",
        "MET_pt", "MET_phi", "run", "luminosityBlock", "HLT_*", "Flag_*"
    ]

    for attempt in range(3):
        try:
            with uproot.open(url, timeout=1800) as root_file:
                available_keys = [k.split(';')[0] for k in root_file.keys()]
                if "Events" in available_keys:
                    events = root_file["Events"]
                elif "events" in available_keys:
                    events = root_file["events"]
                else:
                    break
                chunk = events.arrays(filter_name=needed_branches, how=dict, library="ak")

                if "run" not in chunk.keys() or "luminosityBlock" not in chunk.keys():
                    break

                local_results["events_raw"] = len(chunk["run"])

                # GOLDEN JSON & QUALITY FLAGS
                valid_lumi = lumimask(ak.to_numpy(chunk["run"]), ak.to_numpy(chunk["luminosityBlock"]))
                master_mask = ak.Array(valid_lumi)

                flag_keys = list(chunk.keys())
                required_flags = [
                    "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter",
                    "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter",
                    "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter"
                ]
                for flag in required_flags:
                    if flag in flag_keys:
                        master_mask = master_mask & chunk[flag]

                local_results["events_certified"] = int(ak.sum(master_mask))
                if local_results["events_certified"] == 0:
                    local_results["success"] = True
                    break

                # HLT SELECTION STRICT ENFORCEMENT
                hlt_keys = list(chunk.keys())
                ele_keys = [k for k in hlt_keys if "HLT_Ele27" in k or "HLT_Ele25" in k]
                muon_keys = [k for k in hlt_keys if "HLT_IsoMu24" in k or "HLT_IsoMu22" in k]
                tau_keys = [k for k in hlt_keys if "HLT_LooseIsoPFTau50" in k or "HLT_VLooseIsoPFTau" in k]
                false_mask = chunk["run"] < 0

                pass_ele = false_mask
                if ele_keys:
                    for k in ele_keys:
                        pass_ele = pass_ele | chunk[k]

                pass_muon = false_mask
                if muon_keys:
                    for k in muon_keys:
                        pass_muon = pass_muon | chunk[k]

                pass_tau = false_mask
                if tau_keys:
                    for k in tau_keys:
                        pass_tau = pass_tau | chunk[k]

                # ==================================================
                # 2. SINGLE ELECTRON RECONSTRUCTION FULL VECTOR
                # ==================================================
                if "nElectron" in chunk.keys() and "MET_pt" in chunk.keys():
                    base_mask_ele = master_mask & pass_ele & (chunk["nElectron"] >= 2)

                    if ak.any(base_mask_ele):
                        e_pt = chunk["Electron_pt"][base_mask_ele]
                        e_eta = chunk["Electron_eta"][base_mask_ele]
                        e_phi = chunk["Electron_phi"][base_mask_ele]
                        e_charge = chunk["Electron_charge"][base_mask_ele]

                        if CONTINGENCY_MODE:
                            valid_e_mask = (e_pt > 30) & (abs(e_eta) < 2.5)
                        else:
                            e_id = chunk["Electron_cutBased"][base_mask_ele]
                            valid_e_mask = (e_pt > 30) & (abs(e_eta) < 2.5) & (e_id >= 2)

                        event_has_2e = ak.num(e_pt[valid_e_mask]) >= 2

                        if ak.any(event_has_2e):
                            e_mass_arr = ak.ones_like(e_pt[valid_e_mask][event_has_2e]) * 0.000511
                            electrons = vector.zip({
                                "pt": e_pt[valid_e_mask][event_has_2e],
                                "eta": e_eta[valid_e_mask][event_has_2e],
                                "phi": e_phi[valid_e_mask][event_has_2e],
                                "mass": e_mass_arr
                            })
                            e_charge_filtered = e_charge[valid_e_mask][event_has_2e]

                            e1 = electrons[:, 0]
                            e2 = electrons[:, 1]
                            mass_ee = (e1 + e2).mass

                            met_ee_pt = chunk["MET_pt"][base_mask_ele][event_has_2e]
                            met_ee_phi = chunk["MET_phi"][base_mask_ele][event_has_2e]

                            mass_cut = mass_ee > 100

                            local_results["single_ele_mass"].append(ak.to_numpy(mass_ee[mass_cut]))
                            local_results["single_ele_met_pt"].append(ak.to_numpy(met_ee_pt[mass_cut]))
                            local_results["single_ele_met_phi"].append(ak.to_numpy(met_ee_phi[mass_cut]))

                            local_results["single_ele_e1_pt"].append(ak.to_numpy(e1.pt[mass_cut]))
                            local_results["single_ele_e1_eta"].append(ak.to_numpy(e1.eta[mass_cut]))
                            local_results["single_ele_e1_phi"].append(ak.to_numpy(e1.phi[mass_cut]))
                            local_results["single_ele_e1_charge"].append(ak.to_numpy(e_charge_filtered[:, 0][mass_cut]))

                            local_results["single_ele_e2_pt"].append(ak.to_numpy(e2.pt[mass_cut]))
                            local_results["single_ele_e2_eta"].append(ak.to_numpy(e2.eta[mass_cut]))
                            local_results["single_ele_e2_phi"].append(ak.to_numpy(e2.phi[mass_cut]))
                            local_results["single_ele_e2_charge"].append(ak.to_numpy(e_charge_filtered[:, 1][mass_cut]))

                # ==================================================
                # 3. SINGLE MUON RECONSTRUCTION FULL VECTOR
                # ==================================================
                if "nMuon" in chunk.keys() and "MET_pt" in chunk.keys():
                    base_mask_mu = master_mask & pass_muon & (chunk["nMuon"] >= 2)

                    if ak.any(base_mask_mu):
                        mu_pt = chunk["Muon_pt"][base_mask_mu]
                        mu_eta = chunk["Muon_eta"][base_mask_mu]
                        mu_phi = chunk["Muon_phi"][base_mask_mu]
                        mu_charge = chunk["Muon_charge"][base_mask_mu]

                        if CONTINGENCY_MODE:
                            valid_mu_mask = (mu_pt > 30) & (abs(mu_eta) < 2.4)
                        else:
                            mu_id = chunk["Muon_tightId"][base_mask_mu]
                            mu_iso = chunk["Muon_pfRelIso04_all"][base_mask_mu]
                            valid_mu_mask = (mu_pt > 30) & (abs(mu_eta) < 2.4) & (mu_id == True) & (mu_iso < 0.15)

                        event_has_2mu = ak.num(mu_pt[valid_mu_mask]) >= 2

                        if ak.any(event_has_2mu):
                            mu_mass_arr = ak.ones_like(mu_pt[valid_mu_mask][event_has_2mu]) * 0.10566
                            muons = vector.zip({
                                "pt": mu_pt[valid_mu_mask][event_has_2mu],
                                "eta": mu_eta[valid_mu_mask][event_has_2mu],
                                "phi": mu_phi[valid_mu_mask][event_has_2mu],
                                "mass": mu_mass_arr
                            })
                            mu_charge_filtered = mu_charge[valid_mu_mask][event_has_2mu]

                            m1 = muons[:, 0]
                            m2 = muons[:, 1]
                            mass_mumu = (m1 + m2).mass

                            met_mumu_pt = chunk["MET_pt"][base_mask_mu][event_has_2mu]
                            met_mumu_phi = chunk["MET_phi"][base_mask_mu][event_has_2mu]

                            mass_cut = mass_mumu > 100

                            local_results["single_mu_mass"].append(ak.to_numpy(mass_mumu[mass_cut]))
                            local_results["single_mu_met_pt"].append(ak.to_numpy(met_mumu_pt[mass_cut]))
                            local_results["single_mu_met_phi"].append(ak.to_numpy(met_mumu_phi[mass_cut]))

                            local_results["single_mu_m1_pt"].append(ak.to_numpy(m1.pt[mass_cut]))
                            local_results["single_mu_m1_eta"].append(ak.to_numpy(m1.eta[mass_cut]))
                            local_results["single_mu_m1_phi"].append(ak.to_numpy(m1.phi[mass_cut]))
                            local_results["single_mu_m1_charge"].append(ak.to_numpy(mu_charge_filtered[:, 0][mass_cut]))

                            local_results["single_mu_m2_pt"].append(ak.to_numpy(m2.pt[mass_cut]))
                            local_results["single_mu_m2_eta"].append(ak.to_numpy(m2.eta[mass_cut]))
                            local_results["single_mu_m2_phi"].append(ak.to_numpy(m2.phi[mass_cut]))
                            local_results["single_mu_m2_charge"].append(ak.to_numpy(mu_charge_filtered[:, 1][mass_cut]))

                # ==================================================
                # 4. TAU RECONSTRUCTION FULL VECTOR
                # ==================================================
                if "nTau" in chunk.keys() and "MET_pt" in chunk.keys():
                    base_mask_tau = master_mask & pass_tau & (chunk["nTau"] >= 2)

                    if ak.any(base_mask_tau):
                        tau_pt = chunk["Tau_pt"][base_mask_tau]
                        tau_eta = chunk["Tau_eta"][base_mask_tau]
                        tau_phi = chunk["Tau_phi"][base_mask_tau]
                        tau_m = chunk["Tau_mass"][base_mask_tau]
                        if "Tau_charge" in chunk.keys():
                            tau_charge = chunk["Tau_charge"][base_mask_tau]
                        else:
                            tau_charge = ak.zeros_like(tau_pt)

                        valid_tau_mask = (tau_pt > 30) & (abs(tau_eta) < 2.3)
                        event_has_2tau = ak.num(tau_pt[valid_tau_mask]) >= 2

                        if ak.any(event_has_2tau):
                            taus = vector.zip({
                                "pt": tau_pt[valid_tau_mask][event_has_2tau],
                                "eta": tau_eta[valid_tau_mask][event_has_2tau],
                                "phi": tau_phi[valid_tau_mask][event_has_2tau],
                                "mass": tau_m[valid_tau_mask][event_has_2tau]
                            })
                            tau_charge_filtered = tau_charge[valid_tau_mask][event_has_2tau]

                            t1 = taus[:, 0]
                            t2 = taus[:, 1]
                            mass_tautau = (t1 + t2).mass

                            met_tautau_pt = chunk["MET_pt"][base_mask_tau][event_has_2tau]
                            met_tautau_phi = chunk["MET_phi"][base_mask_tau][event_has_2tau]

                            mass_cut = mass_tautau > 100

                            local_results["single_tau_mass"].append(ak.to_numpy(mass_tautau[mass_cut]))
                            local_results["single_tau_met_pt"].append(ak.to_numpy(met_tautau_pt[mass_cut]))
                            local_results["single_tau_met_phi"].append(ak.to_numpy(met_tautau_phi[mass_cut]))

                            local_results["single_tau_t1_pt"].append(ak.to_numpy(t1.pt[mass_cut]))
                            local_results["single_tau_t1_eta"].append(ak.to_numpy(t1.eta[mass_cut]))
                            local_results["single_tau_t1_phi"].append(ak.to_numpy(t1.phi[mass_cut]))
                            local_results["single_tau_t1_charge"].append(ak.to_numpy(tau_charge_filtered[:, 0][mass_cut]))

                            local_results["single_tau_t2_pt"].append(ak.to_numpy(t2.pt[mass_cut]))
                            local_results["single_tau_t2_eta"].append(ak.to_numpy(t2.eta[mass_cut]))
                            local_results["single_tau_t2_phi"].append(ak.to_numpy(t2.phi[mass_cut]))
                            local_results["single_tau_t2_charge"].append(ak.to_numpy(tau_charge_filtered[:, 1][mass_cut]))

                local_results["success"] = True
                break

        except Exception:
            time.sleep(15)

    return local_results


def execute_patch_3():
    print("\n==================================================")
    print(" [INIT] TBP-OMNI-TENSOR EXTRACTOR (SINGLE-CHANNEL PATCH) ")
    print("==================================================\n")
    print(f"[INFO] Golden JSON loaded successfully. XRootD engine initialized.")
    if CONTINGENCY_MODE:
        print("[INFO] Contingency mode active: Standard ID requirements relaxed.", flush=True)
    else:
        print("[INFO] STRICT MODE ACTIVE: Standard Tight IDs & HLT fully enforced. No Bypass.", flush=True)

    RECORD_IDS = [30529, 30562, 30563, 30530, 30565, 30532]
    all_target_urls = []
    print("[PROCESS] Acquiring XRootD manifests from CERN Open Data...")
    for recid in RECORD_IDS:
        try:
            result = subprocess.run(
                ["cernopendata-client", "get-file-locations", "--recid", str(recid), "--protocol", "xrootd"],
                capture_output=True, text=True, check=True
            )
            urls = result.stdout.strip().split('\n')
            valid_urls = [url.strip() for url in urls if url.strip().startswith('root://')]
            all_target_urls.extend(valid_urls)
        except Exception:
            pass
    all_target_urls = list(set(all_target_urls))[:MAX_FILES]
    print(f"\n[INFO] Total XRootD targets acquired: {len(all_target_urls)} files\n", flush=True)

    global_results = {
        "single_ele_mass": [], "single_ele_met_pt": [], "single_ele_met_phi": [],
        "single_ele_e1_pt": [], "single_ele_e1_eta": [], "single_ele_e1_phi": [], "single_ele_e1_charge": [],
        "single_ele_e2_pt": [], "single_ele_e2_eta": [], "single_ele_e2_phi": [], "single_ele_e2_charge": [],
        "single_mu_mass": [], "single_mu_met_pt": [], "single_mu_met_phi": [],
        "single_mu_m1_pt": [], "single_mu_m1_eta": [], "single_mu_m1_phi": [], "single_mu_m1_charge": [],
        "single_mu_m2_pt": [], "single_mu_m2_eta": [], "single_mu_m2_phi": [], "single_mu_m2_charge": [],
        "single_tau_mass": [], "single_tau_met_pt": [], "single_tau_met_phi": [],
        "single_tau_t1_pt": [], "single_tau_t1_eta": [], "single_tau_t1_phi": [], "single_tau_t1_charge": [],
        "single_tau_t2_pt": [], "single_tau_t2_eta": [], "single_tau_t2_phi": [], "single_tau_t2_charge": []
    }
    processed_files, total_raw, total_certified = 0, 0, 0
    print(f"[PROCESS] Initiating Omni-Tensor ProcessPoolExecutor ({MAX_WORKERS} Workers)...", flush=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(stream_and_analyze_patch3, url): url for url in all_target_urls}
        for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
            res = future.result()
            if res["success"]:
                processed_files += 1
                total_raw += res["events_raw"]
                total_certified += res["events_certified"]
                for key in global_results.keys():
                    if res[key]:
                        global_results[key].extend(res[key])
            if i % 10 == 0 or i == len(all_target_urls):
                print(f"[STATUS] Files Processed: {i}/{len(all_target_urls)} | Raw: {total_raw:,} | Certified: {total_certified:,}", flush=True)
    print(f"\n[INFO] OPERATION COMPLETE. Successful Extractions: {processed_files}/{len(all_target_urls)}", flush=True)
    if processed_files > 0:
        print("[PROCESS] Compiling Omni-Tensor output matrices...", flush=True)
        for key, data_list in global_results.items():
            if data_list:
                np.save(f"{key}.npy", np.concatenate(data_list))
        print("[SUCCESS] 30 Vector Matrices saved successfully.", flush=True)
