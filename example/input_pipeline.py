# pipeline.py
# pipeline.py
from config import (
   PDB_DIR, CIF_DIR, H5_DIR, JSON_DIR, CSV_PATH,
   NUM_COPIES, MAX_LENGTH
)

import subprocess
import sys

STEPS = [
   ("PDB → CIF + CSV", "python 1_extract_chains.py"),
   ("PDB → H5",        "python 2_pdb2jax.py"),
   ("CIF + CSV → JSON","python 3_generate_json.py"),
]

def run_step(name, command):
   print(f"\n{'='*60}")
   print(f"{name}")
   print(f"{'='*60}")
   
   result = subprocess.run(command, shell=True)
   
   if result.returncode != 0:
       print(f"\n Step failed: {name}")
       sys.exit(1)
       
   print(f" Finished: {name}")

def main():
   print("\n Starting full pipeline\n")
   
   for name, cmd in STEPS:
       run_step(name, cmd)

   print("\n All pipeline steps completed successfully!")

if __name__ == "__main__":
   main()



# # pipeline.py
# from config import PDB_ROOT_DIR
# import subprocess
# import sys
# import os
# from pathlib import Path

# STEPS = [
#     ("PDB → CIF + CSV", "python 1_extract_chains.py"),
#     ("PDB → H5",        "python 2_pdb2jax.py"),
#     ("CIF + CSV → JSON","python 3_generate_json.py"),
# ]

# def run_step(name, command, env):
#     print(f"\n{'='*60}")
#     print(f"▶ {name}")
#     print(f"{'='*60}")

#     result = subprocess.run(command, shell=True, env=env)

#     if result.returncode != 0:
#         print(f"\n Step failed: {name}")
#         sys.exit(1)

#     print(f" Finished: {name}")


# def has_pdb_files(folder: Path):
#     return any(folder.glob("*.pdb"))


# def main():
#     print("\n Starting full pipeline\n")

#     cluster_dirs = [
#         d for d in PDB_ROOT_DIR.iterdir()
#         if d.is_dir() and has_pdb_files(d)
#     ]

#     if not cluster_dirs:
#         print(f"No cluster folders with PDB files found in: {PDB_ROOT_DIR}")
#         sys.exit(1)

#     cluster_dirs = sorted(cluster_dirs)

#     print(f"Found {len(cluster_dirs)} cluster folders:")
#     for d in cluster_dirs:
#         print(f" - {d.name}")

#     for cluster_dir in cluster_dirs:
#         print(f"\n\n{'#'*70}")
#         print(f" Processing cluster: {cluster_dir.name}")
#         print(f" Input PDB folder: {cluster_dir}")
#         print(f"{'#'*70}")

#         env = os.environ.copy()
#         env["CURRENT_PDB_DIR"] = str(cluster_dir)

#         for name, cmd in STEPS:
#             run_step(name, cmd, env)

#     print("\n All cluster pipeline steps completed successfully!")


# if __name__ == "__main__":
#     main()