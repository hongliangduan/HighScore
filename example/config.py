# config.py
from pathlib import Path

# ===== Base directory =====
BASE_DIR = Path(__file__).resolve().parent

# ===== Input =====
PDB_DIR = BASE_DIR / "/data/pdb"

# ===== Outputs =====
CIF_DIR = BASE_DIR / "/data/pdb_chain_cifs"
H5_DIR = BASE_DIR / "/data/h5"
JSON_DIR = BASE_DIR / "/data/json"

# ===== Intermediate files =====
CSV_PATH = BASE_DIR / "/data/chainid.csv"
# ===== H5 parameters =====
NUM_COPIES = 1
MAX_LENGTH = 3072

# ===== Cycle setting =====
HEAD_TO_TAIL = False
STANDARDIZE_MODIFIED_RESIDUES_FOR_SCORING = False
