# Suite Integration (Easylab Suite)

This folder contains the integration contract for bundling this repo into the `easylab-suite` desktop launcher.

## What the suite expects
- Streamlit entrypoint at `app.py` (bundled into the suite under `apps/qpcr-analysis/app.py`).
- Supporting files: `qpcr_core.py`, `requirements.txt`, `sample-data/`, `.streamlit/` (optional).
- The suite starts the server on `http://127.0.0.1:8501` by default.

## Module metadata
See `suite/module.json`.
