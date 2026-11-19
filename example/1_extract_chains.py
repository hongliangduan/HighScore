from config import (
    PDB_DIR,
    CIF_DIR,
    CSV_PATH,
    STANDARDIZE_MODIFIED_RESIDUES_FOR_SCORING,
)

from Bio import PDB
from Bio.PDB import Structure, Model, Chain
from Bio.PDB import PDBParser, MMCIFIO
import os
import json
import pandas as pd
from tqdm import tqdm
import multiprocessing as mp
from pathlib import Path
from typing import Dict, Iterable, Tuple
from ccd_residues import resolve_polymer_residue

def get_sequence_and_modifications(chain):
    """Return parent sequence plus AF3 modifications for one chain."""
    sequence = []
    modifications = []
    for residue in chain:
        resname = residue.get_resname().upper()
        resolved = resolve_polymer_residue(resname)
        if resolved is not None:
            one_letter, is_modified = resolved
            sequence.append(one_letter)
            if is_modified and not STANDARDIZE_MODIFIED_RESIDUES_FOR_SCORING:
                
                modifications.append({
                    'ptmType': resname,
                    'ptmPosition': len(sequence),
                })
    return ''.join(sequence), modifications

def assign_json_chain_ids(chain_sequences: Dict[str, str]) -> Dict[str, str]:
    """Map original PDB chain IDs to JSON chain IDs.

    For protein-peptide complexes, the longest chain is treated as protein A and
    the next chain as peptide B. This also maps common peptide chain X to B.
    """
    ordered = sorted(
        chain_sequences.items(),
        key=lambda item: (-len(item[1]), item[0]),
    )
    return {
        original_chain_id: chr(ord('A') + idx)
        for idx, (original_chain_id, _) in enumerate(ordered)
    }

def rewrite_atom_site_chain_ids(cif_path: str, json_chain_id: str) -> None:
    """Normalize chain IDs and polymer residue fields in a generated CIF.

    Biopython writes non-standard polymer residues such as HYP as HETATM with
    an empty label_seq_id. AlphaFold 3 expects polymer residues to have a
    numeric, continuous label_seq_id, so normalize recognized polymer residues
    while leaving unrelated ligands untouched.
    """
    with open(cif_path, 'r') as f:
        lines = f.readlines()

    rewritten = []
    residue_indices = {}
    next_residue_index = 1
    for line in lines:
        parts = line.split()
        if parts and parts[0] in {'ATOM', 'HETATM'} and len(parts) > 16:
            residue_name = parts[5].upper()
            if resolve_polymer_residue(residue_name) is not None:
                residue_key = (parts[6], parts[15], parts[5])
                if residue_key not in residue_indices:
                    residue_indices[residue_key] = next_residue_index
                    next_residue_index += 1
                parts[0] = 'ATOM'
                parts[8] = str(residue_indices[residue_key])
            parts[6] = json_chain_id
            parts[16] = json_chain_id
            line = ' '.join(parts) + '\n'
        rewritten.append(line)

    with open(cif_path, 'w') as f:
        f.writelines(rewritten)

def process_single_pdb(args):
    input_pdb, output_dir_cif = args
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("structure", input_pdb)
        base_name = os.path.splitext(os.path.basename(input_pdb))[0]
        
        original_chain_sequences = {}
        original_chain_modifications = {}
        
        for chain in structure[0]:
            chain_id = chain.id
            sequence, modifications = get_sequence_and_modifications(chain)
            original_chain_sequences[chain_id] = sequence
            original_chain_modifications[chain_id] = modifications

        chain_id_map = assign_json_chain_ids(original_chain_sequences)
        chain_sequences = {}
        chain_modifications = {}

        for chain in structure[0]:
            original_chain_id = chain.id
            json_chain_id = chain_id_map[original_chain_id]
            sequence = original_chain_sequences[original_chain_id]
            chain_sequences[json_chain_id] = sequence
            chain_modifications[json_chain_id] = original_chain_modifications[original_chain_id]
            new_structure = Structure.Structure("new_structure")
            new_model = Model.Model(0)
            new_structure.add(new_model)
            chain_copy = chain.copy()
            chain_copy.id = json_chain_id
            new_model.add(chain_copy)
            
            cif_io = MMCIFIO()
            cif_io.set_structure(new_structure)
            cif_output = os.path.join(output_dir_cif, f"{base_name}_chain_{json_chain_id}.cif")
            cif_io.save(cif_output)
            rewrite_atom_site_chain_ids(cif_output, json_chain_id)

            if original_chain_id != json_chain_id:
                print(
                    f"{base_name}: renamed chain {original_chain_id} -> "
                    f"{json_chain_id} (length={len(sequence)})"
                )
        
        return base_name, chain_sequences, chain_modifications
    except Exception as e:
        print(f"\nError processing {input_pdb}: {str(e)}")
        return None, None, None

def main():

    input_dir = str(PDB_DIR)
    output_dir_cif = str(CIF_DIR)

    # Create output directory
    os.makedirs(output_dir_cif, exist_ok=True)
    
    # Get all PDB files
    pdb_files = list(Path(input_dir).glob("*.pdb"))
    
    # Prepare parameters for process pool
    args = [(str(f), output_dir_cif) for f in pdb_files]
    
    # Process files using process pool
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = list(tqdm(
            pool.imap(process_single_pdb, args),
            total=len(pdb_files),
            desc="Processing PDB files"
        ))
    
    # Collect results
    sequences_dict = {}
    modifications_dict = {}
    for base_name, chain_sequences, chain_modifications in results:
        if base_name is not None:
            sequences_dict[base_name] = chain_sequences
            modifications_dict[base_name] = chain_modifications
    
    # Find all possible chain IDs and sort them by custom order
    all_chain_ids = set()
    for complex_data in sequences_dict.values():
        all_chain_ids.update(complex_data.keys())
    
    all_chain_ids = sorted(list(all_chain_ids))
    
    # Create DataFrame
    rows = []
    for complex_name, chain_data in sequences_dict.items():
        row = {'complex': complex_name}
        for chain_id in all_chain_ids:
            row[f'chain_{chain_id}_seq'] = chain_data.get(chain_id, '')
            modifications = modifications_dict[complex_name].get(chain_id, [])
            row[f'chain_{chain_id}_modifications'] = json.dumps(modifications)
        rows.append(row)
    
    df = pd.DataFrame(rows)
    cols = ['complex'] + [col for col in df.columns if col != 'complex']
    df = df[cols]

    df.to_csv(CSV_PATH, index=False)

    print(f"\nSequence information has been saved to {CSV_PATH}")



if __name__ == "__main__":
    main()
