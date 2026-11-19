from config import CSV_PATH, JSON_DIR, CIF_DIR, DISULFIDE_TOPOLOGY_CSV


import json
import pandas as pd
import os
import argparse
import ast

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


def _is_blank(value):
    """Return True for empty CSV cells, including pandas NaN values."""
    return pd.isna(value) or str(value).strip() == ''


def get_row_modifications(row):
    """Read non-standard residues from the metadata CSV.

    Positions follow the AlphaFold 3 convention: 1-based positions in the
    parent (standard-letter) sequence.  Empty columns deliberately produce an
    empty list, preserving the old natural-amino-acid behavior.
    """
    residue_column = 'unstandard_residue'
    position_column = 'modified_position'
    if residue_column not in row.index or position_column not in row.index:
        return []

    residues_blank = _is_blank(row[residue_column])
    positions_blank = _is_blank(row[position_column])
    if residues_blank and positions_blank:
        return []
    if residues_blank != positions_blank:
        raise ValueError(
            'unstandard_residue and modified_position must both be filled '
            'or both be empty'
        )

    residue_names = [
        item.strip().upper()
        for item in str(row[residue_column]).split(',')
        if item.strip()
    ]
    try:
        positions = [
            int(item.strip())
            for item in str(row[position_column]).split(',')
            if item.strip()
        ]
    except ValueError as exc:
        raise ValueError(
            f'Invalid modified_position: {row[position_column]!r}'
        ) from exc

    if len(residue_names) != len(positions):
        raise ValueError(
            'unstandard_residue and modified_position have different lengths: '
            f'{residue_names!r} vs {positions!r}'
        )

    return [
        {'ptmType': residue_name, 'ptmPosition': position}
        for residue_name, position in zip(residue_names, positions)
    ]


def get_chain_modifications(row, chain_id):
    """Read modifications emitted by 1_extract_chains.py, if present."""
    column = f'chain_{chain_id}_modifications'
    if column not in row.index or _is_blank(row[column]):
        return []
    try:
        modifications = json.loads(row[column])
    except json.JSONDecodeError as exc:
        raise ValueError(f'Invalid {column}: {row[column]!r}') from exc
    if not isinstance(modifications, list):
        raise ValueError(f'{column} must contain a JSON list')
    return modifications


def get_disulfide_pairs(row):
    """Parse the optional s2s_cycle column as residue-number pairs.

    Values are written as ``7-28,13-33`` (or ``7:28,13:33``).  Missing
    values preserve the previous JSON output with no explicit bonds.
    """
    if 's2s_cycle' not in row.index or _is_blank(row['s2s_cycle']):
        return []
    pairs = []
    for item in str(row['s2s_cycle']).replace(';', ',').split(','):
        item = item.strip()
        if not item:
            continue
        fields = item.replace(':', '-').split('-')
        if len(fields) != 2:
            raise ValueError(f'Invalid s2s_cycle value: {row["s2s_cycle"]!r}')
        try:
            first, second = (int(field.strip()) for field in fields)
        except ValueError as exc:
            raise ValueError(
                f'Invalid s2s_cycle value: {row["s2s_cycle"]!r}'
            ) from exc
        if first < 1 or second < 1:
            raise ValueError(f'Invalid s2s_cycle residue numbers: {item!r}')
        pairs.append([first, second])
    return pairs


def load_disulfide_topologies(topology_csv):
    """Load per-structure topology encoded as [[chain, res1, res2, ...]]."""
    if not topology_csv or not os.path.exists(topology_csv):
        return {}
    topology_df = pd.read_csv(topology_csv)
    required = {'structure_id', 'disulfide_chain_res'}
    if not required.issubset(topology_df.columns):
        raise ValueError(f'{topology_csv} must contain columns {sorted(required)}')
    topologies = {}
    for _, topology_row in topology_df.iterrows():
        structure_id = str(topology_row['structure_id']).strip()
        raw_value = topology_row['disulfide_chain_res']
        if not structure_id or _is_blank(raw_value):
            continue
        try:
            topology = ast.literal_eval(str(raw_value))
        except (SyntaxError, ValueError) as exc:
            raise ValueError(
                f'{structure_id}: invalid disulfide_chain_res {raw_value!r}'
            ) from exc
        if not isinstance(topology, list):
            raise ValueError(f'{structure_id}: topology must be a list')
        topologies[structure_id] = topology
    return topologies

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

def generate_json_files(
    csv_path, output_dir, cif_dir, head_to_tail=False,
    disulfide_topology_csv=None,
):
    """Generate JSON files from CSV file"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Read CSV file
    df = pd.read_csv(csv_path)
    disulfide_topologies = load_disulfide_topologies(disulfide_topology_csv)
    json_count = 0
    
    # Process each row
    for _, row in df.iterrows():
        # New ncAA metadata uses pdb_id; retain support for the old CSV format.
        complex_name = row['pdb_id'] if 'pdb_id' in row.index else row['complex']
        
        # Get sequences for all chains
        chain_sequences = get_chain_sequences(row)
        if not chain_sequences and 'sequence' in row.index and not _is_blank(row['sequence']):
            # The supplied sequence column is the parent sequence with ncAAs
            # replaced by their origin_residue.  Its template is chain A.
            chain_sequences = [('A', str(row['sequence']).strip().upper())]
        else:
            chain_sequences = [
                (chain_id, seq) for chain_id, seq in chain_sequences
                if chain_id in {'A', 'B'}
            ]

        modifications = get_row_modifications(row)
        if modifications and len(chain_sequences) != 1:
            raise ValueError(
                f'{complex_name}: row-level modifications require one protein '
                'sequence (the sequence column)'
            )
        
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

            chain_modifications = (
                modifications if modifications else get_chain_modifications(row, chain_id)
            )
            for modification in chain_modifications:
                position = modification['ptmPosition']
                if not isinstance(position, int) or not 0 < position <= len(sequence):
                    raise ValueError(
                        f'{complex_name} chain {chain_id}: invalid modification '
                        f'position {position!r} for sequence length {len(sequence)}'
                    )
                
            sequences.append({
                "protein": {
                    "id": chain_id,
                    "sequence": sequence,
                    "modifications": chain_modifications,
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

        # A cyclic bond is added only when explicitly requested.  Keeping the
        # default false preserves all previous linear evaluations.
        bonded_atom_pairs = []
        if head_to_tail:
            for sequence_entry in sequences:
                protein = sequence_entry['protein']
                if protein['id'] == 'A':
                    bonded_atom_pairs.append([
                        ['A', 1, 'N'],
                        ['A', len(protein['sequence']), 'C'],
                    ])
                    break
        # Disulfide bridges are explicit SG-SG covalent bonds.  They are
        # independent of head-to-tail cyclisation and therefore must not be
        # encoded as an N(1)-C(n) bond.
        topology = disulfide_topologies.get(complex_name)
        if topology is None:
            topology = [[1, *sum(get_disulfide_pairs(row), [])]] if get_disulfide_pairs(row) else []
        chain_ids = [entry['protein']['id'] for entry in sequences]
        for chain_topology in topology:
            if len(chain_topology) < 3 or (len(chain_topology) - 1) % 2:
                raise ValueError(
                    f'{complex_name}: topology must be [chain, res1, res2, ...]'
                )
            chain_number = int(chain_topology[0])
            if not 1 <= chain_number <= len(chain_ids):
                raise ValueError(f'{complex_name}: invalid chain {chain_number}')
            chain_id = chain_ids[chain_number - 1]
            for first, second in zip(chain_topology[1::2], chain_topology[2::2]):
                bonded_atom_pairs.append([
                    [chain_id, int(first), 'SG'],
                    [chain_id, int(second), 'SG'],
                ])
        
        # Create complete JSON data
        json_data = {
            "dialect": "alphafold3",
            "version": 1,
            "name": complex_name,
            "sequences": sequences,
            "modelSeeds": [40],
            "bondedAtomPairs": bonded_atom_pairs or None,
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
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--head_to_tail',
        action='store_true',
        help='Add the explicit A:1:N to A:last:C cyclization bond.',
    )
    parser.add_argument(
        '--disulfide_topology_csv',
        default=str(DISULFIDE_TOPOLOGY_CSV),
        help='CSV mapping structure_id to disulfide_chain_res.',
    )
    args = parser.parse_args()
    generate_json_files(
        csv_path=str(CSV_PATH),
        output_dir=str(JSON_DIR),
        cif_dir=str(CIF_DIR),
        head_to_tail=args.head_to_tail,
        disulfide_topology_csv=args.disulfide_topology_csv,
    )
