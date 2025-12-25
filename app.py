# qPCR Analysis App (Streamlit) — Paste support + auto Cq detection + replicate averages
# Run: streamlit run app.py

import streamlit as st
import pandas as pd
import numpy as np
import re
from pathlib import Path
from datetime import datetime

from qpcr_core import (
    VALID_COLS,
    DEFAULT_TOP_CONC,
    DEFAULT_DILUTION_FACTOR,
    _guess_label_col,
    _guess_rep_col,
    _best_numeric_column_as_cq,
    coerce_columns,
    fit_standard_curve,
    quantify_samples,
    mark_outliers,
    replicate_stats,
    normalize_to_ref,
    download_excel,
    read_from_paste,
    read_from_file,
    _serial_dilution,
    make_std_curve_figure,
)

APP_DIR = Path(__file__).resolve().parent
EXAMPLE_WELLS_PATH = APP_DIR / "sample-data" / "qpcr_example.csv"

st.set_page_config(page_title="qPCR Analysis", page_icon="🧬", layout="centered")

st.markdown(
    """
    <style>
        @import url("https://fonts.googleapis.com/css2?family=Chakra+Petch:wght@500;600;700&family=IBM+Plex+Mono:wght@400;500;600&family=Space+Grotesk:wght@400;500;600;700&display=swap");

        :root{
          --bg: #F7F7F8;
          --surface: #FFFFFF;
          --surface-2: #FBFBFC;
          --text: #111113;
          --text-2: #5E5E66;
          --text-3: #8A8A94;
          --border: #E7E7EA;
          --border-2: #D7D7DD;

          --accent: #4F7CF7;
          --accent-weak: rgba(79,124,247,0.14);
          --focus: rgba(79,124,247,0.25);

          --radius-sm: 10px;
          --radius-md: 12px;

          --s-1: 4px;
          --s-2: 8px;
          --s-3: 12px;
          --s-4: 16px;
          --s-5: 24px;
          --s-6: 32px;

          --shadow-1: 0 1px 2px rgba(0,0,0,0.06);

          --font-display: "Chakra Petch", "Space Grotesk", system-ui, sans-serif;
          --font-body: "Space Grotesk", system-ui, sans-serif;
          --font-sans: var(--font-body);
          --font-mono: "IBM Plex Mono", ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", monospace;

          --tracking-display: 0.06em;
          --tracking-label: 0.12em;
        }

        @media (prefers-color-scheme: dark){
          :root:not([data-theme="light"]){
            --bg: #0B0B0C;
            --surface: #141417;
            --surface-2: #101012;
            --text: #F1F1F3;
            --text-2: #B7B7BF;
            --text-3: #8F8F99;
            --border: #232326;
            --border-2: #2F2F35;

            --accent-weak: rgba(79,124,247,0.18);
            --focus: rgba(79,124,247,0.30);
            --shadow-1: 0 1px 2px rgba(0,0,0,0.35);

            color-scheme: dark;
          }
        }

        [data-theme="dark"]{
          --bg: #0B0B0C;
          --surface: #141417;
          --surface-2: #101012;
          --text: #F1F1F3;
          --text-2: #B7B7BF;
          --text-3: #8F8F99;
          --border: #232326;
          --border-2: #2F2F35;

          --accent-weak: rgba(79,124,247,0.18);
          --focus: rgba(79,124,247,0.30);
          --shadow-1: 0 1px 2px rgba(0,0,0,0.35);

          color-scheme: dark;
        }

        html, body, [data-testid="stAppViewContainer"]{
          background: var(--bg);
          color: var(--text);
          font-family: var(--font-sans);
          font-weight: 500;
          letter-spacing: 0.01em;
        }

        h1, h2, h3, h4, h5{
          font-family: var(--font-display);
          font-weight: 700;
          text-transform: uppercase;
          letter-spacing: var(--tracking-display);
        }

        .block-container{
          max-width: 720px;
          padding-top: 32px;
          padding-bottom: 64px;
        }

        [data-testid="stSidebar"]{
          background: var(--surface);
          border-right: 1px solid var(--border);
          color: var(--text);
        }

        [data-testid="stSidebar"] h1,
        [data-testid="stSidebar"] h2,
        [data-testid="stSidebar"] h3,
        [data-testid="stSidebar"] h4{
          color: var(--text);
        }

        [data-testid="stSidebar"] [data-testid="stWidgetLabel"] > label{
          color: var(--text);
          font-family: var(--font-mono);
          font-size: 12px;
          font-weight: 600;
          letter-spacing: var(--tracking-label);
          text-transform: uppercase;
        }

        [data-testid="stSidebar"] .stCaption,
        [data-testid="stSidebar"] .stMarkdown,
        [data-testid="stSidebar"] .stMarkdown p,
        [data-testid="stSidebar"] .stMarkdown span,
        [data-testid="stSidebar"] .stMarkdown li{
          color: var(--text-2);
        }

        [data-testid="stSidebar"] .stSlider > div > div > span{
          color: var(--text);
        }

        [data-testid="stSidebar"] .stRadio div[role="radiogroup"]{
          display: grid;
          gap: 6px;
        }

        [data-testid="stSidebar"] .stRadio div[role="radiogroup"] label{
          border: 1px solid var(--border);
          background: var(--surface-2);
          border-radius: 10px;
          padding: 6px 10px;
        }

        [data-testid="stSidebar"] .stRadio div[role="radiogroup"] label span{
          color: var(--text);
          font-weight: 600;
          opacity: 1;
        }

        [data-testid="stSidebar"] .stRadio div[role="radiogroup"] label:has(input:checked){
          border-color: var(--accent);
          background: var(--accent-weak);
          box-shadow: 0 0 0 2px var(--focus);
        }

        [data-testid="stSidebar"] input::placeholder,
        [data-testid="stSidebar"] textarea::placeholder{
          color: var(--text-3);
        }

        .streamlit-expanderHeader{
          color: var(--text);
          font-weight: 600;
        }

        [data-testid="stMetricLabel"]{
          color: var(--text-2);
          font-family: var(--font-mono);
          letter-spacing: var(--tracking-label);
          text-transform: uppercase;
          font-weight: 600;
        }

        [data-testid="stMetricValue"]{
          color: var(--text);
        }

        a{
          color: var(--accent);
          text-decoration: underline;
          text-underline-offset: 0.18em;
          text-decoration-color: var(--border-2);
        }

        a:hover{
          text-decoration-color: var(--accent);
        }

        code, pre{
          font-family: var(--font-mono);
        }

        .stButton > button{
          height: 40px;
          border-radius: var(--radius-sm);
          border: 1px solid var(--border);
          background: var(--surface);
          color: var(--text);
          font-weight: 500;
          transition: background-color 140ms ease-out, border-color 140ms ease-out, box-shadow 140ms ease-out, transform 140ms ease-out, filter 140ms ease-out;
        }

        .stButton > button:hover{
          background: var(--surface-2);
        }

        .stButton > button:active{
          transform: translateY(1px);
        }

        .stButton > button:focus-visible{
          outline: none;
          border-color: var(--border-2);
          box-shadow: 0 0 0 2px var(--focus);
        }

        input, textarea, select{
          border-radius: var(--radius-sm) !important;
          border: 1px solid var(--border) !important;
          background: var(--surface) !important;
          color: var(--text) !important;
        }

        input:focus, textarea:focus, select:focus{
          border-color: var(--border-2) !important;
          box-shadow: 0 0 0 2px var(--focus) !important;
        }

        .hero{
          padding: 16px 20px;
          border-radius: var(--radius-md);
          background: var(--surface);
          color: var(--text);
          border: 1px solid var(--border);
          margin-bottom: 12px;
        }

        .hero h4{
          margin: 6px 0 6px 0;
          font-weight: 700;
          text-transform: uppercase;
          letter-spacing: var(--tracking-display);
        }

        .hero p{
          margin: 0;
          color: var(--text-2);
        }

        .pill{
          display: inline-flex;
          align-items: center;
          gap: 6px;
          padding: 4px 10px;
          border-radius: 999px;
          background: var(--accent-weak);
          border: 1px solid var(--border);
          color: var(--text);
          font-weight: 600;
          font-size: 12px;
          letter-spacing: 0;
          font-family: var(--font-mono);
          letter-spacing: var(--tracking-label);
          text-transform: uppercase;
        }

        .stButton > button{
          font-family: var(--font-mono);
          letter-spacing: var(--tracking-label);
          text-transform: uppercase;
          font-weight: 600;
        }
    </style>
    """,
    unsafe_allow_html=True,
)

# ------------- UI -------------
st.sidebar.title("qPCR Analysis")
st.sidebar.caption("Load example wells, upload a file, or paste a table.")

input_mode = st.sidebar.radio("Input source", ["Example (sample-data/qpcr_example.csv)", "Upload file", "Paste table"], index=0)

uploaded_file = None
pasted_text = ""
if input_mode == "Upload file":
    uploaded_file = st.sidebar.file_uploader("CSV / TSV / Excel", type=["csv", "tsv", "txt", "xlsx"])
elif input_mode == "Paste table":
    st.sidebar.caption("Paste Excel/Sheets selection or CSV/TSV (first row must be headers).")
    pasted_text = st.sidebar.text_area(
        "Paste table",
        height=220,
        placeholder="Plate\tWell\tGene\tType\tLabel\tReplicate\tGroup\tCq\nPlate 1\tA1\thcar2\tSample\tA3\t1\tUNDETERMINED\t30.23937102\n..."
    )

st.sidebar.markdown("---")
outlier_thresh = st.sidebar.slider("Outlier ΔCq threshold", 0.1, 3.0, 0.75, 0.05)
plate_scope = st.sidebar.selectbox("Fit standard curve by", ["Gene (all plates)", "Gene × Plate"])
st.sidebar.markdown("---")
ref_gene = st.sidebar.text_input("Reference gene", value="gapdh")
quant_mode = st.sidebar.radio(
    "Quantification mode",
    ["Absolute (std curve)", "ΔΔCt (relative)"],
    horizontal=True,
    help="Absolute uses standard curves; ΔΔCt uses 2^-ΔΔCt relative expression."
)

st.title("qPCR Analysis")
st.markdown(
    """
    <div class="hero">
        <div class="pill">Standard-curve workflow</div>
        <h4>Clean → replicate QC → fit → quantify → export</h4>
        <p>Load example wells, paste instrument output, flag outliers, prefill serial dilutions, and export a tidy Excel report.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

with st.expander("Need to reformat your wells? AI helper prompt"):
    st.markdown(
        """
[ChatGPT](https://chat.openai.com/) · [Gemini](https://gemini.google.com/app) · [Grok](https://grok.com/)

**Prompt**
```text
Convert my table to CSV with headers: Plate, Well, Gene, Type, Label, Replicate, Cq. Normalize Well to A1-style, Replicate numeric, Cq numeric or NA, keep all rows, no invented data. Output CSV only.
```
Paste the returned CSV into the uploader or sidebar paste box.
        """
    )

# Load
try:
    if input_mode.startswith("Example"):
        if not EXAMPLE_WELLS_PATH.exists():
            st.error(f"Missing example file: {EXAMPLE_WELLS_PATH}")
            st.stop()
        df_raw = pd.read_csv(EXAMPLE_WELLS_PATH)
    elif input_mode == "Upload file":
        if uploaded_file is None:
            st.info("Upload a CSV/Excel with columns: Plate, Well, Gene, Type, Label, Replicate, (Group optional), Cq.")
            st.stop()
        df_raw = read_from_file(uploaded_file)
    else:
        if not pasted_text.strip():
            st.info("Paste your table in the sidebar. Your last column (even if unnamed) will be auto-detected as Cq.")
            st.stop()
        df_raw = read_from_paste(pasted_text)
except Exception as e:
    st.error(f"Could not read the data: {e}")
    st.stop()

# Column mapping UI (allows manual selection of label/Cq/replicate headers)
st.sidebar.markdown("---")
st.sidebar.subheader("Column mapping")

cols_available = list(df_raw.columns)
label_default = _guess_label_col(df_raw)
cq_default = "Cq" if "Cq" in df_raw.columns else _best_numeric_column_as_cq(df_raw) or cols_available[-1]
rep_guess = _guess_rep_col(df_raw)

label_col = st.sidebar.selectbox("Sample label column", options=cols_available, index=cols_available.index(label_default) if label_default in cols_available else 0)
cq_col = st.sidebar.selectbox("Cq column", options=cols_available, index=cols_available.index(cq_default) if cq_default in cols_available else len(cols_available)-1)
rep_options = ["<none>"] + cols_available
rep_col_choice = st.sidebar.selectbox("Replicate column (optional)", options=rep_options, index=rep_options.index(rep_guess) if rep_guess in cols_available else 0)

df_mapped = df_raw.copy()
df_mapped["Label"] = df_raw[label_col]
df_mapped["Cq"] = df_raw[cq_col]
if rep_col_choice != "<none>":
    df_mapped["Replicate"] = df_raw[rep_col_choice]
else:
    df_mapped["Replicate"] = 1

raw_df = coerce_columns(df_mapped)

# Keep column + outliers
raw_df["keep"] = True
raw_df = mark_outliers(raw_df, threshold=outlier_thresh)

k1, k2, k3, k4 = st.columns(4)
k1.metric("Wells loaded", raw_df.shape[0])
k2.metric("Genes", raw_df["Gene"].nunique())
k3.metric("Standards", int((raw_df["Type"].str.lower()=="standard").sum()))
k4.metric("Outliers flagged", int(raw_df["Outlier"].sum()))

st.subheader("1) Review & clean wells")
st.caption("Toggle **keep** to drop bad replicates. Outliers by ΔCq are flagged.")
edited = st.data_editor(
    raw_df,
    key="editor",
    num_rows="dynamic",
    column_config={
        "keep": st.column_config.CheckboxColumn("keep", default=True),
        "Outlier": st.column_config.CheckboxColumn("Outlier", disabled=True),
        "DeltaCq": st.column_config.NumberColumn("ΔCq", format="%.3f", disabled=True),
        "Cq": st.column_config.NumberColumn("Cq", format="%.3f")
    },
    hide_index=True,
    width="stretch",
)
clean_df = edited.copy()
clean_df.loc[~clean_df["keep"], "Cq"] = np.nan

# Collapse exact duplicate wells (same Plate/Well/Gene/Type/Label); average Cq across kept replicates.
extra_cols = [c for c in clean_df.columns if c not in VALID_COLS + ["keep","Outlier","DeltaCq"]]
well_keys = ["Plate","Well","Gene","Type","Label"]
agg_map = {"Cq": "mean", "keep": "any", "DeltaCq": "mean", "Outlier": "any"}
agg_map.update({c: "first" for c in extra_cols})
collapsed_df = clean_df.groupby(well_keys, as_index=False, dropna=False).agg(agg_map)

# ------------- replicate averages -------------
st.subheader("2) Replicate averages (auto-calculated)")
rep_stats = replicate_stats(collapsed_df)
st.dataframe(rep_stats, width="stretch", height=440)

# ------------- standards mapping -------------
if quant_mode == "Absolute (std curve)":
    st.subheader("3) Standards map (Label → Concentration)")
    std_labels = clean_df.loc[(clean_df["Type"].str.lower()=="standard"), "Label"].dropna().unique().tolist()
    std_labels = sorted(std_labels, key=lambda x: (re.sub(r"\D","",x)=="", re.sub(r"\D","",x), x))
    auto_fill = _serial_dilution(std_labels, DEFAULT_TOP_CONC, DEFAULT_DILUTION_FACTOR, highest_first=True)
    default_map = pd.DataFrame({"Label": std_labels, "Concentration": auto_fill})
    st.caption("Default mapping uses a 4-fold serial dilution from the top standard; edit if your plate differs.")

    with st.expander("Auto-fill serial dilution"):
        top = st.number_input("Top concentration", min_value=0.0, value=DEFAULT_TOP_CONC, step=1.0, format="%.6f")
        factor = st.number_input("Dilution factor", min_value=1.0, value=DEFAULT_DILUTION_FACTOR, step=0.5, format="%.3f")
        order = st.selectbox("Order of labels", ["Std1 highest → StdN lowest", "Std1 lowest → StdN highest"])
        if st.button("Fill mapping from settings"):
            default_map["Concentration"] = _serial_dilution(
                std_labels,
                top,
                factor,
                highest_first=order.startswith("Std1 highest")
            )

    map_df = st.data_editor(
        default_map,
        key="stdmap",
        column_config={"Concentration": st.column_config.NumberColumn("Concentration", format="%.6f")},
        hide_index=True,
        width="stretch",
    )
    if map_df["Concentration"].isna().any():
        st.warning("Some standard concentrations are missing. Fill them before fitting curves.")

if quant_mode == "Absolute (std curve)":
    # ------------- fit curves -------------
    st.subheader("4) Fit standard curves")
    if plate_scope == "Gene × Plate":
        std_input = collapsed_df[(collapsed_df["Type"].str.lower()=="standard") & (collapsed_df["keep"])].copy()
        std_input["Gene"] = std_input["Gene"].astype(str) + " | " + std_input["Plate"].astype(str)
    else:
        std_input = collapsed_df[(collapsed_df["Type"].str.lower()=="standard") & (collapsed_df["keep"])].copy()

    curves_df, std_points = fit_standard_curve(std_input, map_df)
    st.dataframe(curves_df, width="stretch", height=260)

    if curves_df.empty:
        st.info("No standard curves to fit yet. Add standards (Type=Standard) and fill the standards map.")
    else:
        bad = curves_df[(curves_df["R2"] < 0.98) | (~curves_df["Efficiency_%"].between(90, 110))]
        if bad.shape[0] > 0:
            st.warning("Some curves have R² < 0.98 or efficiency outside 90–110%. Consider revising outliers or concentrations.")

        show_plots = st.toggle("Show curve plots", value=True)
        if show_plots:
            st.caption("Each point: mean Cq at that standard level. Line: linear fit. Efficiency computed from slope.")
            for _, row in curves_df.iterrows():
                gene = row["Gene"]
                sub = std_points[std_points["Gene"] == gene]
                if sub.empty or pd.isna(row["slope"]):
                    continue
                x = sub["log10_conc"].values
                y = sub["meanCq"].values
                fig = make_std_curve_figure(
                    gene=str(gene),
                    x=x,
                    y=y,
                    slope=float(row["slope"]),
                    intercept=float(row["intercept"]),
                    r2=float(row["R2"]),
                    eff=float(row["Efficiency_%"]),
                )
                st.pyplot(fig, width="stretch")

    # ------------- quantify samples -------------
    st.subheader("5) Quantify samples")
    if plate_scope == "Gene × Plate":
        samp_input = collapsed_df[(collapsed_df["Type"].str.lower()=="sample") & (collapsed_df["keep"])].copy()
        samp_input["Gene"] = samp_input["Gene"].astype(str) + " | " + samp_input["Plate"].astype(str)
    else:
        samp_input = collapsed_df[(collapsed_df["Type"].str.lower()=="sample") & (collapsed_df["keep"])].copy()

    quant_df = quantify_samples(samp_input, curves_df)
    st.dataframe(quant_df.head(60), width="stretch", height=440)

    # ------------- normalize -------------
    st.subheader("6) Normalize to reference gene")
    norm_df = normalize_to_ref(quant_df, ref_gene=ref_gene)
    if norm_df["RefQty"].isna().all():
        st.info(f"Reference gene '{ref_gene}' missing in standards/samples; add a ref gene to see normalized values.")
    st.dataframe(norm_df, width="stretch", height=440)

    # Per-well normalized view (attach all metadata for export)
    per_well_norm = quant_df.merge(norm_df[["Label","RefQty"]], on="Label", how="left")
    per_well_norm["Norm_Qty"] = per_well_norm["Quantity"] / per_well_norm["RefQty"]
else:
    # ΔΔCt workflow
    st.subheader("4) ΔΔCt (relative)")
    ref_rows = collapsed_df[collapsed_df["Gene"].str.lower() == ref_gene.lower()].copy()
    ref_mean = ref_rows.groupby("Label", dropna=False)["Cq"].mean().reset_index().rename(columns={"Cq":"RefCq"})
    dd_df = collapsed_df.merge(ref_mean, on="Label", how="left")
    dd_df = dd_df[dd_df["Gene"].str.lower() != ref_gene.lower()].copy()
    dd_df["DeltaCt"] = dd_df["Cq"] - dd_df["RefCq"]
    all_labels = sorted(dd_df["Label"].dropna().unique().tolist())
    calibrators = st.multiselect("Calibrator labels", options=all_labels, default=all_labels[:1])
    exclude = st.multiselect("Exclude from calibrator", options=calibrators, default=[])
    calib_set = [l for l in calibrators if l not in exclude]
    calib_means = (
        dd_df[dd_df["Label"].isin(calib_set)]
        .groupby("Gene")["DeltaCt"]
        .mean()
        .reset_index()
        .rename(columns={"DeltaCt":"CalibDeltaCt"})
    )
    dd_df = dd_df.merge(calib_means, on="Gene", how="left")
    dd_df["DeltaDeltaCt"] = dd_df["DeltaCt"] - dd_df["CalibDeltaCt"]
    dd_df["FoldChange"] = 2 ** (-dd_df["DeltaDeltaCt"])
    st.dataframe(dd_df[["Plate","Well","Gene","Label","Cq","RefCq","DeltaCt","CalibDeltaCt","DeltaDeltaCt","FoldChange"]], width="stretch", height=440)
    dd_csv = dd_df[["Plate","Well","Gene","Label","Cq","RefCq","DeltaCt","CalibDeltaCt","DeltaDeltaCt","FoldChange"]].to_csv(index=False).encode()
    st.download_button("Download ΔΔCt table (CSV)", dd_csv, file_name="ddct_results.csv", mime="text/csv")

    # placeholders to keep downstream variables defined
    map_df = pd.DataFrame()
    curves_df = pd.DataFrame()
    std_points = pd.DataFrame()
    quant_df = pd.DataFrame()
    norm_df = dd_df.rename(columns={"FoldChange":"Norm_Qty"})
    per_well_norm = dd_df[["Plate","Well","Gene","Label","Cq","RefCq","DeltaCt","CalibDeltaCt","DeltaDeltaCt","FoldChange"]].copy()
    per_well_norm["slope"] = np.nan
    per_well_norm["intercept"] = np.nan
    per_well_norm["pred_log10Q"] = np.nan
    per_well_norm["Quantity"] = np.nan
    per_well_norm["RefQty"] = np.nan
    per_well_norm["Norm_Qty"] = per_well_norm["FoldChange"]

# Build a metadata-rich per-well export (one row per well)
meta_cols = ["Plate","Well","Gene","Type","Label","Cq","keep","DeltaCq","Outlier"] + extra_cols
meta_frame = collapsed_df[meta_cols].copy()

# Attach curve + quantity metadata so the Excel "PerWell_Normalized" sheet is fully informative.
merge_cols = ["Plate","Well","Gene","Label"]
per_well_export_cols = merge_cols + ["slope","intercept","pred_log10Q","Quantity","RefQty","Norm_Qty"]
export_norm_df = (
    meta_frame
    .merge(per_well_norm[per_well_export_cols], on=merge_cols, how="left")
    # Ensure only one row per well/gene/label (collapse any accidental duplicate copies)
    .drop_duplicates(subset=["Plate","Well","Gene","Label"])
    # Order by gene first so reference and targets are grouped, not interleaved.
    .sort_values(["Gene","Plate","Well","Label"])
)

# ------------- export -------------
st.subheader("7) Export")
# Pass per-sample normalized table separately from the per-well export table
excel_bytes = download_excel(clean_df, rep_stats, map_df, curves_df, std_points, quant_df, norm_df, export_norm_df)
st.download_button(
    label="Download Excel report",
    data=excel_bytes,
    file_name=f"qpcr_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx",
    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    type="primary",
)

st.success("Done. Paste/upload → clean → averages → fit → quantify → normalize → export.")
