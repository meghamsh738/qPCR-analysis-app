# qPCR Analysis App

Streamlit dashboard for quickly reviewing qPCR plates end-to-end – from pasted wells tables through replicate QC, standard curve fitting, absolute quantity estimation, ΔΔCt normalisation, and Excel export. The app is designed to work directly from the tables produced by most qPCR instruments, so you can validate a run without pre-processing in spreadsheets.

Status: this Streamlit app is the canonical version. A prior React/FastAPI experiment that lacked chart export has been removed to keep the repo lean.

## Highlights
- Paste CSV/TSV/Excel-style well tables or vendor exports; delimiter detection and case-insensitive column matching keep the input flexible.
- Automatically reshapes vendor-specific columns, derives replicate numbers, and flags negative or control wells.
- Per-gene tabs surface replicate statistics (mean Cq, SD, ΔCq), standard curve plots with slope/intercept/R²/efficiency, and concentration back-calculations.
- Calibrator-aware normalisation plus relative expression via 2^-ΔΔCt with selectable reference genes and calibrator pools.
- One-click Excel export (`qpcr_results.xlsx`) that bundles cleaned wells, replicate summaries, standard curve points & plots, quantities, and relative-expression tables.
- Includes `mock_wells.csv` so you can try the workflow without real instrument output.

## Screenshot

Example run using the bundled `sample-data/qpcr_example.csv` (select **Example** in the sidebar) and default options (refreshed Dec 12, 2025):

![qPCR analysis dashboard screenshot](screenshots/example_run.png)

Additional views from the full workflow:

![Overview](screenshots/overview.png)
![Replicates](screenshots/replicates.png)
![Standards](screenshots/standards.png)
![Curves](screenshots/curves.png)
![Quant + Normalize](screenshots/quant_normalize.png)
![Export](screenshots/export.png)

## Quick Start
1. **Prerequisites:** Python 3.9+ and `pip`. Creating a virtual environment is recommended.
2. **Install dependencies**
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate  # Windows: .venv\Scripts\activate
   pip install -r requirements.txt
   ```
3. **Run the dashboard**
   ```bash
   streamlit run app.py
   ```
4. Streamlit prints a local URL (typically http://localhost:8501). Open it in your browser.

## Using the Dashboard
The workflow is split across tabs: **Wells → Replicates → Standards → Curves → Quantify → Normalize → Export**.

1. **Wells:** Paste or upload wells. Minimum columns are `Gene`, `Label`, and `Cq`; optional metadata (Plate, Well, Type, Replicate) is auto-detected. Use **Example** in the sidebar to load `sample-data/qpcr_example.csv`.
2. **Replicates:** Review mean/SD/CV and verify outliers.
3. **Standards:** Fill the standards map (Label → Concentration) and optionally auto-fill a serial dilution series.
4. **Curves:** Standards are fitted per gene (or per gene×plate); review slope/intercept/R²/efficiency and inspect the plot.
5. **Quantify / Normalize:** Convert sample Cq values into absolute quantities using the fitted curve, then normalize to the reference gene.
6. **Export:** Download an Excel workbook that bundles all derived tables.

Use the download button at the bottom of the page to collect all derived tables (plus rendered standard-curve plots) into a single Excel workbook for record keeping.

## Data format & quick AI helper
- Best-fit columns (CSV/TSV/Excel): `Plate, Well, Gene, Type, Label, Replicate, Cq`. Optional extras are ignored safely. See `sample-data/qpcr_example.csv`.
- To reshape quickly, use: [ChatGPT](https://chat.openai.com/), [Gemini](https://gemini.google.com/app), [Grok](https://grok.com/).
- Prompt: "Convert my table to CSV with headers: Plate, Well, Gene, Type, Label, Replicate, Cq. Keep my data, no invented rows, output UTF-8 CSV text only."
- Save as `qpcr_example.csv` (any name) and upload/paste. Visual helper: `screenshots/data-format-helper.svg`.

## Tips
- Paste data exactly as exported; the app trims whitespace, harmonises case, and interprets common "NA"/"Undetermined" tokens automatically.
- Keep at least two standard levels with known concentrations per gene to enable curve fitting.
- When experimenting, use the **Example** input source to load `sample-data/qpcr_example.csv` (includes standards and samples).
- If you only want replicate QC (no standards), `mock_wells.csv` is a minimal wells table you can paste/upload.

## How calculations are done
1) **Clean & outliers** – ΔCq per gene/label is median‐based; wells with ΔCq above the sidebar threshold are flagged. If a well is unchecked, its Cq becomes NaN.  
2) **Replicate collapse** – Wells sharing Plate/Well/Gene/Type/Label are averaged (mean Cq) before any downstream math. Extra columns keep the first value.  
3) **Standards** – For each gene (or Gene×Plate if selected):  
   - log10(conc) from the standards map vs mean Cq → linear fit (slope, intercept).  
   - R² and efficiency = (10^(−1/slope) − 1)·100.  
4) **Quantities (absolute)** – For each kept sample well:  
   - pred_log10Q = (Cq − intercept) / slope; Quantity = 10^pred_log10Q.  
5) **Normalisation to ref gene** – Per label, RefQty is the mean Quantity of the chosen ref gene; Norm_Qty = Quantity / RefQty.  
6) **Exports** – `PerWell_Normalized` carries one row per Plate/Well/Gene/Label with slope/intercept/pred_log10Q/Quantity/RefQty/Norm_Qty plus your metadata.  

### ΔΔCt option?
The app currently performs absolute quantification with ref-gene normalisation (step 5). A 2^-ΔΔCt readout is not yet exposed in the UI; we can add a toggle to pick a calibrator label and emit ΔΔCt/fold change if you want.

## Development Notes
- The Streamlit script lives in `app.py`; update it and re-run `streamlit run app.py` to see changes.
- Dependencies are listed in `requirements.txt`. Pin additional libraries there if you extend the workflow.
- No dedicated tests ship with the project yet; consider adding regression tests (e.g., with `pytest`) if you automate calculations outside the UI.
