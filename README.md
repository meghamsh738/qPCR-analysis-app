# EasyLab Next · qPCR Analysis

Separate development version for the approved EasyLab improvement plan. The existing Suite, standalone apps and original checkouts are preserved. This branch has not been released or installed over an existing app.

## Run locally

From this Next module checkout:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run_next.py
```

Open `http://127.0.0.1:18501`. The launcher refuses an occupied port. In the Next worktree, data, settings and exports live under `easylab-next/runtime`, separately from the existing applications. `EASYLAB_DATA_PATH`, `EASYLAB_SETTINGS_PATH` and `EASYLAB_EXPORTS_PATH` must remain inside the separate `EASYLAB_RUNTIME_ROOT`. No HOME override or copied account/profile is needed. The Suite's own launcher uses reserved Next ports and isolated profiles.

Use the **Open practice experiment** button for fictional data, or import CSV/TSV/TXT/XLSX Results files. Explicitly confirm columns, identities, retained wells and methods. Imported rows and original files remain available in the experiment package.

**Independent reference sets.** In Methods, select the same target assays and biological calibrators, then name each reference set. Put each reference gene in its own set to compare them separately; select multiple genes in a single set to combine them. Choose a set in Results and Graphs & statistics. The assay-by-assay workbook uses the chosen set; detailed exports include the set ID and name for each result. Phase-1 packages open under their original reference set and retain their saved results.

## Workflow

1. **Import:** preview original rows and confirm columns. Native Results headers and separate sample-information tables are supported.
2. **Map samples:** confirm biological SampleId, Group and optional PairId/Tissue/Batch. Default identities include plate, label and supplied context. Explicit IDs can link the same sample across target/reference plates; conflicting metadata is held for review.
3. **Review wells:** retain undetermined/excluded observations, give well-specific exclusion reasons, and inspect technical replicate summaries. Warnings do not automatically remove observations.
4. **Methods:** review assay/plate-specific standards, select references and calibrator samples, and explicitly review efficiency and range assumptions. The named personal preset starts at 16384 AU with fourfold dilutions. Non-contiguous standard selections remain supported.
5. **Results:** compare plate-specific standard-curve ratios, standard-curve calibrator-relative fold changes and ΔΔCq fold changes. Technical Cqs are averaged before biological summaries; calibrator samples have equal weight. Missing or withheld inputs remain visible with reasons.
6. **Graphs & statistics:** biological dots and means, optional uncertainty, editable labels/colours/order/dimensions and linear/log axes. Welch, explicitly paired tests and one-way ANOVA with pooled-variance comparisons use selected families and Holm adjustment. Axes do not transform statistical values. Multigroup repeated-measures and two-factor designs are outside this feature set.
7. **Save & export:** versioned `.easylab.json` packages, local saves, CSV/TSV, Excel, PNG/SVG/PDF figures and selected-results PDF reports. Save a package to preserve sources, decisions, results and figure settings. New experiments retain display/recipe preferences but require new scientific decisions.

The **Assay-by-assay manual review** Excel option preserves compact lower-case assay sheets, standards/curves, retained technical Cqs and a `normalised` sheet. It contains live formulas with initial cached values. Excel edits are separate from the app experiment. To return edits to the app, export a reviewed flat well table and confirm mappings; directly importing this review workbook is deliberately rejected so its original raw block cannot overwrite reviewed changes. Formula editing requires a compatible spreadsheet engine; native Excel recalculation and installed desktop behavior are separate verification steps.

## Development checks

```sh
.venv/bin/python -m pip install pytest==8.4.2
.venv/bin/python -m pytest -q
npm ci
npm run test:e2e
```

Browser tests start their own isolated server on 18501 and refuse an existing listener. Stop only a server you own before running them. Tests use fictional fixtures and temporary data. `qpcr_core.py` remains for legacy regression checks; the Next app uses `qpcr_workflow.py` and `qpcr_statistics.py`.

The standalone Electron shell has a distinct Next name/profile and uses this checkout's environment in development. Packaging/runtime validation on macOS and Windows remains separate from the source/browser checks. No replacement of an installed app is implied.

License: All Rights Reserved.
