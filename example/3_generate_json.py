from config import CSV_PATH, JSON_DIR, CIF_DIR


import json
import pandas as pd
import os

ATOM_SITE_LABEL_ASYM_ID = 6
ATOM_SITE_AUTH_ASYM_ID = 16

def format_msa_sequence(sequence):
    """Format MSA sequence"""
    return f">query\n{sequence}\n"

def get_chain_sequences(row):
    """Get all non-empty chain sequences from row data"""
    chain_sequences = []
    chain_columns = [col for col in row.index if col.startswith('chain_') and col.endswith('_seq')]
    chain_columns = sorted(chain_columns, key=lambda col: col.split('_')[1])
    for col in chain_columns:
        if pd.notna(row[col]) and row[col] != '':
            chain_id = col.split('_')[1]
            chain_sequences.append((chain_id, row[col]))
    return chain_sequences

def rewrite_atom_site_chain_ids(cif_path, json_chain_id):
    """Force _atom_site label_asym_id/auth_asym_id to match JSON chain id."""
    with open(cif_path, 'r') as f:
        lines = f.readlines()

    rewritten = []
    for line in lines:
        parts = line.split()
        if parts and parts[0] in {'ATOM', 'HETATM'} and len(parts) > ATOM_SITE_AUTH_ASYM_ID:
            parts[ATOM_SITE_LABEL_ASYM_ID] = json_chain_id
            parts[ATOM_SITE_AUTH_ASYM_ID] = json_chain_id
            line = ' '.join(parts) + '\n'
        rewritten.append(line)

    with open(cif_path, 'w') as f:
        f.writelines(rewritten)

def first_atom_chain_ids(mmcif_text):
    for line in mmcif_text.splitlines():
        parts = line.split()
        if parts and parts[0] in {'ATOM', 'HETATM'} and len(parts) > ATOM_SITE_AUTH_ASYM_ID:
            return line, parts[ATOM_SITE_LABEL_ASYM_ID], parts[ATOM_SITE_AUTH_ASYM_ID]
    raise ValueError('template mmCIF contains no ATOM/HETATM rows')

def _sequence_id_and_data(sequence_entry):
    chain_type, chain_data = next(iter(sequence_entry.items()))
    return chain_type, chain_data['id'], chain_data

def validate_json_template_chain_ids(json_path):
    """Validate JSON sequence IDs, template chain IDs, and template mappings."""
    with open(json_path, 'r') as f:
        data = json.load(f)

    for sequence_entry in data.get('sequences', []):
        _, sequence_id, chain_data = _sequence_id_and_data(sequence_entry)
        sequence = chain_data['sequence']
        expected_indices = list(range(len(sequence)))
        templates = chain_data.get('templates', [])
        if len(sequence_id) != 1:
            raise ValueError(f"{json_path}: chain id {sequence_id!r} is not a single chain id")
        if not templates:
            raise ValueError(f"{json_path}: chain {sequence_id} has no templates")

        template = templates[0]
        query_indices = template.get('queryIndices', [])
        template_indices = template.get('templateIndices', [])
        if query_indices != expected_indices:
            raise ValueError(
                f"{json_path}: chain {sequence_id} queryIndices must be 0..{len(sequence)-1}, "
                f"got length {len(query_indices)}"
            )
        if template_indices != expected_indices:
            raise ValueError(
                f"{json_path}: chain {sequence_id} templateIndices must be 0..{len(sequence)-1}, "
                f"got length {len(template_indices)}"
            )

        if 'mmcif' in template:
            mmcif_text = template['mmcif']
        else:
            with open(template['mmcifPath'], 'r') as f:
                mmcif_text = f.read()

        first_atom, label_asym_id, auth_asym_id = first_atom_chain_ids(mmcif_text)
        if label_asym_id != sequence_id:
            raise ValueError(
                f"{json_path}: chain {sequence_id} template label_asym_id={label_asym_id}, "
                f"expected {sequence_id}. first ATOM: {first_atom}"
            )
        if auth_asym_id != sequence_id:
            raise ValueError(
                f"{json_path}: chain {sequence_id} template auth_asym_id={auth_asym_id}, "
                f"expected {sequence_id}. first ATOM: {first_atom}"
            )

        print(
            f"[validate_json_template_chain_ids] chain={sequence_id} "
            f"length={len(sequence)} first_atom='{first_atom}' "
            f"label_asym_id={label_asym_id} auth_asym_id={auth_asym_id} "
            f"queryIndices_len={len(query_indices)} "
            f"templateIndices_len={len(template_indices)}"
        )

    return True

def generate_json_files(csv_path, output_dir, cif_dir):
    """Generate JSON files from CSV file"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Read CSV file
    df = pd.read_csv(csv_path)
    json_count = 0
    
    # Process each row
    for _, row in df.iterrows():
        complex_name = row['complex']  # Get name from complex column
        
        # Get sequences for all chains
        chain_sequences = get_chain_sequences(row)
        chain_sequences = [(chain_id, seq) for chain_id, seq in chain_sequences if chain_id in {'A', 'B'}]
        
        if not chain_sequences:  # Skip if no valid chain sequences
            print(f"Warning: {complex_name} has no valid chain sequences")
            continue
        
        # Create a list of all chain sequences
        sequences = []
        for chain_id, sequence in chain_sequences:
            # Build cif file path
            cif_filename = f"{complex_name}_chain_{chain_id}.cif"
            cif_path = os.path.join(cif_dir, cif_filename)
            
            # Check if cif file exists
            if not os.path.exists(cif_path):
                print(f"Warning: {cif_filename} does not exist")
                continue

            rewrite_atom_site_chain_ids(cif_path, chain_id)
                
            sequences.append({
                "protein": {
                    "id": chain_id,
                    "sequence": sequence,
                    "modifications": [],
                    "unpairedMsa": format_msa_sequence(sequence),
                    "pairedMsa": format_msa_sequence(sequence),
                    "templates": [{
                        "mmcifPath": cif_path,
                        "queryIndices": list(range(len(sequence))),
                        "templateIndices": list(range(len(sequence)))
                    }]
                }
            })
        
        if not sequences:  # Skip if no valid sequence data
            print(f"Warning: {complex_name} has no valid sequence data")
            continue
        
        # Create complete JSON data
        json_data = {
            "dialect": "alphafold3",
            "version": 1,
            "name": complex_name,
            "sequences": sequences,
            "modelSeeds": [40],
            "bondedAtomPairs": None,
            "userCCD": None
        }
        
        # Generate output file name - REMOVED _data suffix
        output_filename = f"{complex_name}.json"  # Removed _data suffix to match H5 filename
        output_path = os.path.join(output_dir, output_filename)
        
        # Write JSON file
        with open(output_path, 'w') as f:
            json.dump(json_data, f, indent=2)

        validate_json_template_chain_ids(output_path)
        
        chain_ids = [
            next(iter(sequence_entry.values()))['id']
            for sequence_entry in sequences
        ]
        print(f"Generated JSON file: {output_filename} (chains: {', '.join(chain_ids)})")
        json_count += 1
    
    print(f"\nComplete, generated {json_count} JSON files")


if __name__ == "__main__":
    generate_json_files(
        csv_path=str(CSV_PATH),
        output_dir=str(JSON_DIR),
        cif_dir=str(CIF_DIR)
    )
