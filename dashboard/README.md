# Blood–MSG Transcriptomics Explorer

This Streamlit application is a read-only narrative layer over results already produced by the R pipeline. It does not execute R, train modules, cluster patients, or recompute correlations, P values, or FDR. Patient identifiers are retained only in memory for validated joins and are never displayed.

## Required outputs

The Overview requires `00_pairing/paired_manifest.tsv`. Endotypes requires the Stage 02 endotype assignments and standardized MSG eigengenes. Blood Projection requires MSG eigengenes, Blood projected scores, and the official concordance table. Optional unmatched-sample, PCA, silhouette, and Stage 04 image files degrade to explicit availability messages. Stage `04_endotype_de/` is not required to start or use the three core views.

Only one `OUTPUT_DIR` is read at a time. Until the pipeline publishes immutable run manifests, the UI labels it **Unversioned output directory**; creation of validated, read-only snapshots with an atomically updated `latest.json` is a future pipeline task.

## Biowulf installation and launch

Run on a compute node, **not** the login node. First obtain a node and tunnel:

```bash
sinteractive --tunnel --mem=8g --cpus-per-task=2
hostname
echo "$PORT1"
```

On that compute node, from the repository root:

```bash
module avail python
# module load <institutionally available Python module>
python3 -m venv .venv_dashboard
source .venv_dashboard/bin/activate
python -m pip install -r dashboard/requirements.txt
python --version
python -m pip freeze

export OUTPUT_DIR=/path/to/the/completed/R/output
export STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
python -m streamlit run dashboard/app.py \
  --server.address 127.0.0.1 \
  --server.port "$PORT1" \
  --server.headless true
```

On the local workstation, use the exact command/port emitted by `sinteractive --tunnel`, then open `http://localhost:$PORT1`. Keep both the Slurm job and SSH tunnel alive. Never bind this institutional-data dashboard to `0.0.0.0`.

## Tests

Tests contain clearly synthetic identifiers and no NIH data:

```bash
python -m pytest dashboard/tests
```

For a headless smoke test on a compute node:

```bash
timeout 20s python -m streamlit run dashboard/app.py \
  --server.address 127.0.0.1 --server.port 8501 --server.headless true
```

The app reports missing/corrupt required outputs in the UI. It never falls back to synthetic or older files. `INPUT_DIR` is neither imported nor accessed.
