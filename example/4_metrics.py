import pandas as pd
import json
import os
import numpy as np
from pathlib import Path
from multiprocessing import Pool, cpu_count
from functools import partial
from tqdm import tqdm
from collections import defaultdict
import traceback 
def get_interface_res_from_cif(cif_file, dist_cutoff=10):
    """Get interface residues from CIF file using robust parsing"""
    chain_coords = defaultdict(lambda: defaultdict(dict))
    
    try:
        with open(cif_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    atom_name = line[12:16].strip()
                    if atom_name not in ['"CA"', 'CA']:
                        continue
                    
                    chain_id = line[21:22].strip()
                    # residue_id = int(line[24:28].strip())
                    residue_id_str = line[24:28].strip().split()[0]
                    residue_id = int(residue_id_str)
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    
                    chain_coords[chain_id][residue_id] = np.array([x, y, z])
        
        if "A" not in chain_coords or len(chain_coords) < 2:
            return [], []

        chain_A_res = sorted(chain_coords["A"].keys())
        chain_B_res = sorted(list(set(res for chain, res_dict in chain_coords.items() 
                                    if chain != "A" for res in res_dict.keys())))

        if not chain_A_res or not chain_B_res:
            return [], []

        chain_A_coords = np.array([chain_coords["A"][res] for res in chain_A_res])
        
        # Combine all non-A chains into chain B
        chain_B_coords_list = []
        for chain, res_dict in chain_coords.items():
            if chain != "A":
                for res_id in sorted(res_dict.keys()):
                    chain_B_coords_list.append(res_dict[res_id])
        chain_B_coords = np.array(chain_B_coords_list)

        dist = np.sqrt(np.sum((chain_A_coords[:, None, :] - chain_B_coords[None, :, :])**2, axis=2))
        interface_residues_indices = np.where(dist < dist_cutoff)
        
        interface_A = sorted(set(chain_A_res[i] for i in interface_residues_indices[0]))
        interface_B = sorted(set(chain_B_res[i] for i in interface_residues_indices[1]))
        return interface_A, interface_B
    except Exception:
        # If parsing fails for any reason, return empty lists
        return [], []

def get_metrics_from_json(base_path, target, seed, sample):
    """为特定的 target, seed, 和 sample 获取所有指标"""
    # 动态构建文件夹名称
    seed_folder = f"seed-{seed}_sample-{sample}"
    base_json_path = Path(base_path) / target / seed_folder
    try:
        with open(base_json_path / "summary_confidences.json") as f:
            summary_data = json.load(f)
            complex_iptm = summary_data.get("iptm", np.nan)
            complex_ptm = summary_data.get("ptm", np.nan)
            monomer_ptm = summary_data.get("chain_ptm", [np.nan])[0] 
            ranking_score = summary_data.get("ranking_score", np.nan)
        chain_ca = defaultdict(dict)
        with open(base_json_path / "model.cif", 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    parts = line.split()
                    if len(parts) < 15: continue
                    residue_id = int(parts[8])
                    atom_name = parts[3].strip()
                    plddt = float(parts[14])
                    chain_id = parts[6]
                    if atom_name == "CA":
                        chain_ca[chain_id][residue_id] = plddt
        
        a_chain_ca_plddts = list(chain_ca.get("A", {}).values())
        monomer_ca_plddt = np.mean(a_chain_ca_plddts) if a_chain_ca_plddts else np.nan
        
        all_ca_plddts = [plddt for chain_plddts in chain_ca.values() for plddt in chain_plddts.values()]
        complex_ca_plddt = np.mean(all_ca_plddts) if all_ca_plddts else np.nan
        
        interface_A, interface_B = get_interface_res_from_cif(base_json_path / "model.cif")
        
        with open(base_json_path / "confidences.json") as f:
            conf_data = json.load(f)
            pae_matrix = np.array(conf_data["pae"])
            token_chain_ids = conf_data["token_chain_ids"]
            token_res_ids = conf_data["token_res_ids"]
            
            interface_A_indices = [i for i, (r, c) in enumerate(zip(token_res_ids, token_chain_ids)) if c == "A" and r in interface_A]
            interface_B_indices = [i for i, (r, c) in enumerate(zip(token_res_ids, token_chain_ids)) if c != "A" and r in interface_B]
            chain_A_indices = [i for i, c in enumerate(token_chain_ids) if c == "A"]
            chain_B_indices = [i for i, c in enumerate(token_chain_ids) if c != "A"]

            ipae = np.mean(pae_matrix[np.ix_(interface_A_indices, interface_B_indices)]) if interface_A_indices and interface_B_indices else np.nan
            monomer_pae = np.mean(pae_matrix[np.ix_(chain_A_indices, chain_A_indices)]) if chain_A_indices else np.nan
            pae_interaction = np.mean(pae_matrix[np.ix_(chain_A_indices, chain_B_indices)]) if chain_A_indices and chain_B_indices else np.nan
            complex_pae = np.mean(pae_matrix) if pae_matrix.size > 0 else np.nan
            
        return {
            "HighScore_monomer_ca_plddt": monomer_ca_plddt, "HighScore_monomer_pae": monomer_pae,
            "HighScore_monomer_ptm": monomer_ptm, "HighScore_complex_ca_plddt": complex_ca_plddt,
            "HighScore_complex_pae": complex_pae, "HighScore_complex_ptm": complex_ptm,
            "HighScore_complex_iptm": complex_iptm, "HighScore_pae_interaction": pae_interaction,
            "HighScore_ipae": ipae, "HighScore_ranking_score": ranking_score
        }, None
    except Exception as e:
        return None, f"Error processing {target}/{seed_folder}: {e}"


def process_row(base_path, row):
    target = str(row['target']).strip()

    seed = int(row['seed'])
    sample = int(row['sample'])
    metrics, error = get_metrics_from_json(base_path, target, seed, sample)
    return {'idx': row.name, 'metrics': metrics, 'error': error}


def update_sc_file():
    
    df = pd.read_csv("/data/dataset/output/complex_chain_sequences.csv", sep=',', low_memory=False)
    
    print("Columns in file:", df.columns.tolist())

    if 'target' not in df.columns or 'description' not in df.columns:
        raise ValueError(f"Missing required 'target' or 'description' columns in CSV file\nAvailable columns: {df.columns.tolist()}")
    
    df['target'] = df['target'].astype(str)
    df['description'] = df['description'].astype(str)

    metrics_list = [
        'HighScore_monomer_ca_plddt', 'HighScore_monomer_pae', 'HighScore_monomer_ptm',
        'HighScore_complex_ca_plddt', 'HighScore_complex_pae', 'HighScore_complex_ptm', 
        'HighScore_complex_iptm', 'HighScore_pae_interaction', 'HighScore_ipae', 'HighScore_ranking_score'
    ]
    for metric in metrics_list:
        df[metric] = np.nan
    
    print(f"\nTotal {len(df)} records to process")
    
    base_path = "/data/dataset/output/"

    num_processes = max(1, int(cpu_count() * 0.8))
    print(f"Using {num_processes} processes")
    
    with Pool(num_processes) as pool:
        process_func = partial(process_row, base_path)
        results = list(tqdm(
            pool.imap(process_func, [row for _, row in df.iterrows()]),
            total=len(df),
            desc="Processing"
        ))
    
    failed_records = []
    success_count = {metric: 0 for metric in metrics_list}
    
    for result in results:
        if result['metrics'] is not None:
            for metric, value in result['metrics'].items():
                df.at[result['idx'], metric] = value
                if not pd.isna(value):
                    success_count[metric] += 1
        else:
            failed_records.append(f"line {result['idx']+2}: {result['error']}")
    
    with open('failed_records.txt', 'w') as f:
        f.write(f"Total processed rows: {len(df)}\n")
        f.write(f"Failed rows: {len(failed_records)}\n\n")
        f.write("Detailed failure records:\n")
        for record in failed_records:
            f.write(f"{record}\n")
    
    df.to_csv("/data/dataset/output/complex_chain_sequences.csv", sep=",", index=False)

    print("\nProcessing completion statistics:")
    print(f"Total entries: {len(df)}")
    for metric in metrics_list:
        print(f"Successfully updated {metric} count: {success_count[metric]}")
    print(f"Failed entries: {len(failed_records)}")
    print(f"Failed records written to: failed_records.txt")

if __name__ == "__main__":
    update_sc_file()