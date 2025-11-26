# qPCR Analysis App (Streamlit) — Paste support + auto Cq detection + replicate averages
# Run: streamlit run app.py

import streamlit as st
import pandas as pd
import numpy as np
import io, re
from datetime import datetime
import matplotlib.pyplot as plt

st.set_page_config(page_title="qPCR Analysis", layout="wide")

VALID_COLS = ["Plate", "Well", "Gene", "Type", "Label", "Replicate", "Cq"]
DEFAULT_TOP_CONC = 1000.0
DEFAULT_DILUTION_FACTOR = 4.0

st.markdown(
    """
    <style>
        .hero {
            padding: 1rem 1.25rem;
            border-radius: 14px;
            background: linear-gradient(120deg, #0f172a, #1c2740);
            color: #e2e8f0;
            border: 1px solid #1f2a44;
            margin-bottom: 0.75rem;
        }
        .hero h4 {
            margin: 0.15rem 0 0.25rem 0;
            font-weight: 700;
        }
        .hero p {
            margin: 0;
        }
        .pill {
            display: inline-flex;
            align-items: center;
            gap: 0.35rem;
            padding: 0.25rem 0.7rem;
            border-radius: 999px;
            background: #38bdf8;
            color: #0b1324;
            font-weight: 700;
            font-size: 0.85rem;
            letter-spacing: 0.01em;
        }
    </style>
    """,
    unsafe_allow_html=True,
)

# ------------- utilities -------------
def _to_num(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip().upper()
    if s in ["UNDETERMINED", "NA", "NAN", "INF", "#DIV/0!"]:
        return np.nan
    try:
        return float(s)
    except:
        return np.nan

def _to_int(x):
    try:
        return int(float(x))
    except:
        return np.nan


def _guess_label_col(df: pd.DataFrame) -> str:
    """Pick a likely sample label column."""
    for cand in ["Label", "Sample", "Sample ID", "Sample Name", "SampleID", "SampleName"]:
        for col in df.columns:
            if str(col).strip().lower() == cand.strip().lower():
                return col
    return df.columns[0]


def _guess_rep_col(df: pd.DataFrame) -> str | None:
    for col in df.columns:
        cl = str(col).lower()
        if "rep" in cl:
            return col
    return None

def _best_numeric_column_as_cq(df):
    """Pick the most likely numeric column to be Cq when no 'Cq' header exists.
       Strategy: choose non-core columns with the highest count of numeric-parsable values,
       prefer columns whose variance > 0."""
    candidates = []
    core = set(["Plate","Well","Gene","Type","Label","Replicate","Group","group"])
    for c in df.columns:
        if c in VALID_COLS or c in core:
            continue
        vals = df[c].apply(_to_num)
        nnum = vals.notna().sum()
        var = float(vals.var(skipna=True)) if nnum > 1 else 0.0
        candidates.append((c, nnum, var))
    if not candidates:
        # try Unnamed columns
        for c in df.columns:
            if str(c).lower().startswith("unnamed"):
                vals = df[c].apply(_to_num)
                nnum = vals.notna().sum()
                var = float(vals.var(skipna=True)) if nnum > 1 else 0.0
                candidates.append((c, nnum, var))
    if not candidates:
        return None
    # choose by nnum then variance
    candidates.sort(key=lambda t: (t[1], t[2]))
    return candidates[-1][0] if candidates[-1][1] > 0 else None

def coerce_columns(df: pd.DataFrame) -> pd.DataFrame:
    # Flexible renaming
    rename = {}
    for c in df.columns:
        cl = str(c).strip().lower()
        if cl in ["plate id", "plate"]:
            rename[c] = "Plate"
        elif cl in ["well id", "well"]:
            rename[c] = "Well"
        elif cl in ["gene", "assay"]:
            rename[c] = "Gene"
        elif cl in ["type", "sample type"]:
            rename[c] = "Type"
        elif cl in ["sample", "label", "sample id", "target name", "sample name"]:
            rename[c] = "Label"
        elif cl in ["rep", "replicate"]:
            rename[c] = "Replicate"
        elif cl in ["cq", "ct", "c[t]", "cq average", "cq mean"]:
            rename[c] = "Cq"
    df = df.rename(columns=rename)

    # If Cq missing, try detect from unnamed / numeric column
    if "Cq" not in df.columns:
        guess = _best_numeric_column_as_cq(df)
        if guess is not None:
            df = df.rename(columns={guess: "Cq"})

    # Create required columns if absent
    for k in VALID_COLS:
        if k not in df.columns:
            df[k] = np.nan

    # Order main columns first, keep extras at end
    extras = [c for c in df.columns if c not in VALID_COLS]
    df = df[VALID_COLS + extras]

    # Clean data types
    df["Cq"] = df["Cq"].apply(_to_num)
    df["Replicate"] = df["Replicate"].apply(_to_int)
    for col in ["Plate", "Well", "Gene", "Type", "Label"]:
        df[col] = df[col].astype(str).str.strip()

    return df

def r2_score(y_true, y_pred):
    y_true = np.array(y_true, dtype=float)
    y_pred = np.array(y_pred, dtype=float)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    return 1 - ss_res / ss_tot if ss_tot != 0 else np.nan

def fit_standard_curve(std_df: pd.DataFrame, mapping: pd.DataFrame):
    m = mapping.dropna(subset=["Label", "Concentration"]).copy()
    m["Label"] = m["Label"].astype(str)
    m["Concentration"] = m["Concentration"].astype(float)
    m["log10_conc"] = np.log10(m["Concentration"])
    df = std_df.merge(m, on="Label", how="inner")
    g = df.groupby(["Gene", "Label", "log10_conc"], as_index=False)["Cq"].mean().rename(columns={"Cq":"meanCq"})
    curves = []
    for gene, sub in g.groupby("Gene"):
        if sub.shape[0] < 2:
            curves.append({"Gene": gene, "slope": np.nan, "intercept": np.nan, "R2": np.nan, "Efficiency_%": np.nan, "n_points": sub.shape[0]})
            continue
        x = sub["log10_conc"].values
        y = sub["meanCq"].values
        slope, intercept = np.polyfit(x, y, 1)
        pred = slope * x + intercept
        R2 = r2_score(y, pred)
        eff = (10 ** (-1.0 / slope) - 1.0) * 100.0
        curves.append({
            "Gene": gene,
            "slope": slope,
            "intercept": intercept,
            "R2": R2,
            "Efficiency_%": eff,
            "n_points": sub.shape[0]
        })
    return pd.DataFrame(curves), g

def quantify_samples(samples_df: pd.DataFrame, curves_df: pd.DataFrame):
    df = samples_df.merge(curves_df[["Gene","slope","intercept"]], on="Gene", how="left")
    df["pred_log10Q"] = (df["Cq"] - df["intercept"]) / df["slope"]
    df["Quantity"] = 10 ** df["pred_log10Q"]
    return df

def mark_outliers(df: pd.DataFrame, threshold=0.75):
    df = df.copy()
    grpkey = df["Gene"].astype(str) + "||" + df["Label"].astype(str)
    med = df.groupby(grpkey)["Cq"].transform("median")
    df["DeltaCq"] = (df["Cq"] - med).abs()
    df["Outlier"] = (df["DeltaCq"] > float(threshold)) & (~df["Cq"].isna())
    return df

def replicate_stats(clean_df: pd.DataFrame):
    """Average Cq by (Gene, Type, Label) for kept wells."""
    sub = clean_df[clean_df["keep"]].copy()
    agg = (
        sub.groupby(["Gene","Type","Label"], dropna=False)["Cq"]
        .agg(["mean","std","count"])
        .reset_index()
        .rename(columns={"mean":"AvgCq","std":"SdCq","count":"N"})
    )
    agg["CV%"] = (agg["SdCq"] / agg["AvgCq"]) * 100
    return agg

def normalize_to_ref(q_df: pd.DataFrame, ref_gene: str):
    """Return per-sample aggregates with normalization against ref gene."""
    if q_df.empty:
        cols = ["Gene","Label","Qty_mean","Qty_sd","n","RefQty","Norm_Qty"]
        return pd.DataFrame(columns=cols)
    agg = (
        q_df.groupby(["Gene", "Label"], dropna=False)["Quantity"]
        .agg(Qty_mean="mean", Qty_sd="std", n="count")
        .reset_index()
    )
    ref = (
        agg[agg["Gene"].str.lower() == str(ref_gene).lower()][["Label","Qty_mean"]]
        .rename(columns={"Qty_mean":"RefQty"})
    )
    merged = agg.merge(ref, on="Label", how="left")
    merged["Norm_Qty"] = merged["Qty_mean"] / merged["RefQty"]
    return merged

def download_excel(clean_df, rep_stats, std_map_df, curves_df, std_points_df, quant_df, norm_df, per_well_norm_df):
    out = io.BytesIO()
    with pd.ExcelWriter(out, engine="openpyxl") as xls:
        clean_df.to_excel(xls, sheet_name="Cleaned_Wells", index=False)
        rep_stats.to_excel(xls, sheet_name="Replicate_Cq_Stats", index=False)
        std_map_df.to_excel(xls, sheet_name="Standards_Map", index=False)
        std_points_df.to_excel(xls, sheet_name="StdCurve_Points", index=False)
        curves_df.to_excel(xls, sheet_name="StdCurve_Fits", index=False)
        quant_df.to_excel(xls, sheet_name="PerWell_Quant", index=False)
        norm_df.to_excel(xls, sheet_name="PerSample_Normalized", index=False)
        per_well_norm_df.to_excel(xls, sheet_name="PerWell_Normalized", index=False)

        try:
            from openpyxl.drawing.image import Image as XLImage
            plot_sheet = xls.book.create_sheet("StdCurve_Plots")
            row_cursor = 1
            for _, row in curves_df.iterrows():
                gene = row["Gene"]
                sub = std_points_df[std_points_df["Gene"] == gene]
                if sub.empty or pd.isna(row.get("slope")):
                    continue
                x = sub["log10_conc"].values
                y = sub["meanCq"].values
                xp = np.linspace(min(x), max(x), 100)
                yp = row["slope"] * xp + row["intercept"]
                fig, ax = plt.subplots()
                ax.scatter(x, y, label="Points")
                ax.plot(xp, yp, linestyle="--", label=f"Fit (R²={row['R2']:.3f}, Eff={row['Efficiency_%']:.1f}%)")
                ax.set_xlabel("log10(concentration)")
                ax.set_ylabel("Cq")
                ax.set_title(f"Std curve: {gene}")
                ax.legend()
                buf = io.BytesIO()
                fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
                plt.close(fig)
                buf.seek(0)
                img = XLImage(buf)
                plot_sheet.add_image(img, f"A{row_cursor}")
                row_cursor += 25
        except Exception:
            pass
    out.seek(0)
    return out

# ------------- reading -------------
def _detect_sep(text: str) -> str:
    head = "\n".join(text.strip().splitlines()[:5])
    cand = {'\t': head.count('\t'), ',': head.count(','), ';': head.count(';'), '|': head.count('|')}
    sep = max(cand, key=cand.get)
    return sep if cand[sep] > 0 else r"\s+"

def read_from_paste(pasted: str) -> pd.DataFrame:
    txt = pasted.strip()
    sep = _detect_sep(txt)
    if sep == r"\s+":
        return pd.read_csv(io.StringIO(txt), sep=sep, engine="python")
    return pd.read_csv(io.StringIO(txt), sep=sep)

def read_from_file(uploaded):
    if uploaded.name.lower().endswith(".xlsx"):
        return pd.read_excel(uploaded, engine="openpyxl")
    content = uploaded.read().decode("utf-8", errors="ignore")
    if "\t" in content and (content.count("\t") > content.count(",")):
        return pd.read_csv(io.StringIO(content), sep="\t")
    return pd.read_csv(io.StringIO(content))

def _serial_dilution(labels, top, factor, highest_first=True):
    """Return concentrations for labels based on serial dilution settings."""
    vals = []
    n = len(labels)
    if n == 0:
        return vals
    if highest_first:
        for i in range(n):
            vals.append(top / (factor ** i))
    else:
        for i in range(n):
            vals.append(top * (factor ** i))
    return vals

# ------------- UI -------------
st.sidebar.title("Settings")
input_mode = st.sidebar.radio("Input source", ["Upload file", "Paste table"], horizontal=True)

uploaded_file = None
pasted_text = ""
if input_mode == "Upload file":
    uploaded_file = st.sidebar.file_uploader("CSV or Excel", type=["csv","tsv","txt","xlsx"])
else:
    st.sidebar.caption("Paste Excel/Sheets selection or CSV/TSV (first row must be headers).")
    pasted_text = st.sidebar.text_area(
        "Paste table here",
        height=220,
        placeholder="Plate\tWell\tGene\tType\tLabel\tReplicate\tGroup\tCq\nPlate 1\tA1\thcar2\tSample\tA3\t1\tUNDETERMINED\t30.23937102\n..."
    )

st.sidebar.markdown("---")
outlier_thresh = st.sidebar.slider("Outlier ΔCq threshold", 0.1, 3.0, 0.75, 0.05)
plate_scope = st.sidebar.selectbox("Fit standard curve by", ["Gene (all plates)", "Gene × Plate"])
st.sidebar.markdown("---")
ref_gene = st.sidebar.text_input("Reference gene", value="gapdh")

st.title("qPCR Analysis (Standard-curve based)")
st.markdown(
    """
    <div class="hero">
        <div class="pill">4-fold standards by default</div>
        <h4>QC → standard curves → quant → export</h4>
        <p>Paste or upload wells, flag outliers, prefill serial dilutions, and export a tidy Excel report.</p>
    </div>
    """,
    unsafe_allow_html=True,
)

# Load
try:
    if input_mode == "Upload file":
        if uploaded_file is None:
            st.info("Upload a CSV/Excel with columns: Plate, Well, Gene, Type, Label, Replicate, (Group optional), Cq.")
            st.stop()
        df_raw = read_from_file(uploaded_file)
    else:
        if not pasted_text.strip():
            st.info("Paste your table on the left. Your last column (even if unnamed) will be auto-detected as Cq.")
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
)
clean_df = edited.copy()
clean_df.loc[~clean_df["keep"], "Cq"] = np.nan

# ------------- replicate averages -------------
st.subheader("2) Replicate averages (auto-calculated)")
rep_stats = replicate_stats(clean_df)
st.dataframe(rep_stats)

# ------------- standards mapping -------------
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
)
if map_df["Concentration"].isna().any():
    st.warning("⚠️ Some standard concentrations are missing. Fill them before fitting curves.")

# ------------- fit curves -------------
st.subheader("4) Fit standard curves")
if plate_scope == "Gene × Plate":
    std_input = clean_df[(clean_df["Type"].str.lower()=="standard") & (clean_df["keep"])].copy()
    std_input["Gene"] = std_input["Gene"].astype(str) + " | " + std_input["Plate"].astype(str)
else:
    std_input = clean_df[(clean_df["Type"].str.lower()=="standard") & (clean_df["keep"])].copy()

curves_df, std_points = fit_standard_curve(std_input, map_df)
st.dataframe(curves_df)

bad = curves_df[(curves_df["R2"] < 0.98) | (~curves_df["Efficiency_%"].between(90, 110))]
if bad.shape[0] > 0:
    st.warning("Some curves have R² < 0.98 or efficiency outside 90–110%. Consider revising outliers or concentrations.")

st.caption("Each point: mean Cq at that standard level. Line: linear fit. Efficiency computed from slope.")
for _, row in curves_df.iterrows():
    gene = row["Gene"]
    sub = std_points[std_points["Gene"] == gene]
    if sub.empty or pd.isna(row["slope"]):
        continue
    x = sub["log10_conc"].values
    y = sub["meanCq"].values
    xp = np.linspace(min(x), max(x), 100)
    yp = row["slope"]*xp + row["intercept"]
    fig, ax = plt.subplots()
    ax.scatter(x, y, label="Points")
    ax.plot(xp, yp, linestyle="--", label=f"Fit (R²={row['R2']:.3f}, Eff={row['Efficiency_%']:.1f}%)")
    ax.set_xlabel("log10(concentration)")
    ax.set_ylabel("Cq")
    ax.set_title(f"Standard curve: {gene}")
    ax.legend()
    st.pyplot(fig)

# ------------- quantify samples -------------
st.subheader("5) Quantify samples")
if plate_scope == "Gene × Plate":
    samp_input = clean_df[(clean_df["Type"].str.lower()=="sample") & (clean_df["keep"])].copy()
    samp_input["Gene"] = samp_input["Gene"].astype(str) + " | " + samp_input["Plate"].astype(str)
else:
    samp_input = clean_df[(clean_df["Type"].str.lower()=="sample") & (clean_df["keep"])].copy()

quant_df = quantify_samples(samp_input, curves_df)
st.dataframe(quant_df.head(30))

# ------------- normalize -------------
st.subheader("6) Normalize to reference gene")
norm_df = normalize_to_ref(quant_df, ref_gene=ref_gene)
if norm_df["RefQty"].isna().all():
    st.info(f"Reference gene '{ref_gene}' missing in standards/samples; add a ref gene to see normalized values.")
st.dataframe(norm_df)

# Per-well normalized view (attach all metadata for export)
per_well_norm = quant_df.merge(norm_df[["Label","RefQty"]], on="Label", how="left")
per_well_norm["Norm_Qty"] = per_well_norm["Quantity"] / per_well_norm["RefQty"]

# Build a metadata-rich per-well export (Plate/Well/Gene/Type/Label/Replicate/Cq + extras + keep/outlier)
extra_cols = [c for c in clean_df.columns if c not in VALID_COLS + ["keep","Outlier","DeltaCq"]]
meta_cols = ["Plate","Well","Gene","Type","Label","Replicate","Cq","keep","DeltaCq","Outlier"] + extra_cols
meta_frame = clean_df[meta_cols].copy()

# Attach curve + quantity metadata so the Excel "PerWell_Normalized" sheet is fully informative.
merge_cols = ["Plate","Well","Gene","Label","Replicate"]
per_well_export_cols = merge_cols + ["slope","intercept","pred_log10Q","Quantity","RefQty","Norm_Qty"]
export_norm_df = (
    meta_frame
    .merge(per_well_norm[per_well_export_cols], on=merge_cols, how="left")
    # Ensure only one row per well/gene/label (collapse any accidental duplicate copies)
    .drop(columns=["Replicate"])
    .drop_duplicates(subset=["Plate","Well","Gene","Label","Cq"])
)

# ------------- export -------------
st.subheader("7) Export")
excel_bytes = download_excel(clean_df, rep_stats, map_df, curves_df, std_points, quant_df, export_norm_df, per_well_norm)
st.download_button(
    label="Download Excel report",
    data=excel_bytes,
    file_name=f"qpcr_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx",
    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
)

st.success("Done. Paste/upload → clean → averages → fit → quantify → normalize → export.")
