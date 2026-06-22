# HighScore

HighScore is an AlphaFold 3 based scoring pipeline for cyclic peptide structures. 

## Repository Layout

```text
.
├── run_highscore.py
├── fetch_databases.sh
├── example/
│   ├── 1_extract_chains.py
│   ├── 2_pdb2jax.py
│   ├── 3_generate_json.py
│   ├── 4_metrics.py
│   ├── config.py
│   └── input_pipeline.py
├── docs/
├── src/alphafold3/
├── pyproject.toml
├── requirements.txt
└── environment.yml
```

## Environment

This codebase targets:

- Linux
- Python 3.11
- NVIDIA GPU for inference
- AlphaFold 3 model weights
- AlphaFold 3 sequence/template databases

The packaged dependencies are defined in `pyproject.toml`, `requirements.txt`,
and `environment.yml`.

### Recommended setup

```bash
conda create -n highscore python=3.11
conda activate highscore
conda install -c conda-forge \
  cxx-compiler c-compiler cmake ninja pkg-config make hmmer \
  biopython h5py pandas tqdm pytest
pip install -r dev-requirements.txt
pip install --no-deps -e .
build_data
```

### Optional convenience path

If you prefer, you can still create the same base environment from
`environment.yml`:

```bash
conda env create -f environment.yml
conda activate highscore
pip install -r dev-requirements.txt
pip install --no-deps -e .
build_data
```

`environment.yml` is not the primary installation path. It should be treated as
an optional, minimal base environment only.

If you do not use `environment.yml`, install at least:

- a C/C++ build toolchain
- HMMER binaries (`jackhmmer`, `nhmmer`, `hmmalign`, `hmmsearch`, `hmmbuild`)
- Python dependencies from `pyproject.toml`
- extra helper-script dependencies such as `biopython`, `h5py`, and `pandas`

### Databases

Download AlphaFold 3 databases with:

```bash
bash fetch_databases.sh <DB_DIR>
```

`run_highscore.py` expects AlphaFold 3 style database layout under `--db_dir`.

### Model parameters

You must separately obtain AlphaFold 3 model parameters and pass their path via
`--model_dir`.

## Main Inference Entry Point

The main scoring entry point is:

```bash
python run_highscore.py --output_dir=/path/to/output ...
```

This script supports two modes.

### 1. Single-input mode

Example with one JSON file plus one initial structure H5:

```bash
python run_highscore.py \
  --db_dir=/path/to/alphafold_databases \
  --model_dir=/path/to/alphafold3_model_parameters \
  --json_path=/path/to/input.json \
  --path=/path/to/input.h5 \
  --output_dir=/path/to/output \
  --run_data_pipeline=false \
  --run_inference=true \
  --init_guess=true
```

### 2. Batch mode

Use both `--batch_json_dir` and `--batch_h5_dir`.

Each JSON file must have a matching H5 file with the same basename.

```bash
python run_highscore.py \
  --db_dir=/path/to/alphafold_databases \
  --model_dir=/path/to/alphafold3_model_parameters \
  --batch_json_dir=/path/to/json_dir \
  --batch_h5_dir=/path/to/h5_dir \
  --output_dir=/path/to/output \
  --run_data_pipeline=false \
  --run_inference=true \
  --init_guess=true
```

### Important flags implemented in `run_highscore.py`

- `--json_path`: single JSON input
- `--input_dir`: directory of JSON inputs
- `--batch_json_dir`: batch JSON directory
- `--batch_h5_dir`: batch H5 directory
- `--path`: initial structure H5 path used in single-input mode
- `--output_dir`: required output directory
- `--db_dir`: AlphaFold 3 database root; may be passed multiple times
- `--model_dir`: model parameter directory
- `--run_data_pipeline`: whether to run MSA/template pipeline
- `--run_inference`: whether to run model inference
- `--init_guess`: whether to use the initial structure guess
- `--num_samples`: number of samples per seed
- `--head_to_tail`: cyclic peptide option

## Helper Input-Preparation Scripts

The helper scripts live in `example/`:

- `example/1_extract_chains.py`
- `example/2_pdb2jax.py`
- `example/3_generate_json.py`
- `example/input_pipeline.py`

Their intended workflow is:

1. extract per-chain CIF files and sequence CSV from PDBs
2. convert PDB structures into H5 coordinate tensors
3. generate AlphaFold 3 JSON files using the extracted sequences and CIFs

### Actual commands

Run them from the `example/` directory:

```bash
cd example
python 1_extract_chains.py
python 2_pdb2jax.py
python 3_generate_json.py
```

Or run the combined wrapper:

```bash
cd example
python input_pipeline.py
```

### Current configuration caveat

`example/config.py` contains hard-coded paths written as:

```python
PDB_DIR = BASE_DIR / "/data/pdb"
CIF_DIR = BASE_DIR / "/data/pdb_chain_cifs"
H5_DIR = BASE_DIR / "/data/h5"
JSON_DIR = BASE_DIR / "/data/json"
CSV_PATH = BASE_DIR / "/data/chainid.csv"
```

Because the right-hand side starts with `/`, these resolve to absolute paths
under `/data/...`, not to paths inside this repository. Before using the helper
scripts, update `example/config.py` to match your real input/output locations.

## Metrics Extraction

The metrics helper script in this repository is:

```bash
python example/4_metrics.py
```

It extracts:

- `HighScore_monomer_ca_plddt`
- `HighScore_monomer_pae`
- `HighScore_monomer_ptm`
- `HighScore_complex_ca_plddt`
- `HighScore_complex_pae`
- `HighScore_complex_ptm`
- `HighScore_complex_iptm`
- `HighScore_pae_interaction`
- `HighScore_ipae`
- `HighScore_ranking_score`

### Current caveat

`example/4_metrics.py` is not a general CLI tool yet. It contains hard-coded
input/output paths under `/data/dataset/output/` and expects specific CSV
columns such as `target`, `seed`, `sample`, and `description`. It should be
treated as a project-specific post-processing script unless you adapt it.

## Output Structure

The scoring outputs follow AlphaFold 3 style directories under `--output_dir`,
typically one subdirectory per target name, containing:

- `seed-<seed>_sample-<sample>/`
- `model.cif`
- `confidences.json`
- `summary_confidences.json`

See `docs/output.md` for the detailed AlphaFold 3 output description.

## Additional Documentation

- `docs/installation.md`
- `docs/input.md`
- `docs/output.md`
- `docs/performance.md`
- `docs/known_issues.md`

These documents mostly describe the bundled AlphaFold 3 codebase and remain the
best reference for database layout, JSON schema, and output semantics.


## Reference

For more information about AlphaFold 3, please visit their GitHub repository:
<https://github.com/google-deepmind/alphafold3>


