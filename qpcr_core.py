from __future__ import annotations

import io
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

VALID_COLS = ["Plate", "Well", "Gene", "Type", "Label", "Replicate", "Cq"]
DEFAULT_TOP_CONC = 1000.0
DEFAULT_DILUTION_FACTOR = 4.0
PLOT_BG = "#0b1224"


def _to_num(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip().upper()
    if s in ["UNDETERMINED", "NA", "NAN", "INF", "#DIV/0!"]:
        return np.nan
    try:
        return float(s)
    except Exception:
        return np.nan


def _to_int(x):
    try:
        return int(float(x))
    except Exception:
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


def _best_numeric_column_as_cq(df: pd.DataFrame):
    """Pick the most likely numeric column to be Cq when no 'Cq' header exists.
       Strategy: choose non-core columns with the highest count of numeric-parsable values,
       prefer columns whose variance > 0."""
    candidates = []
    core = set(["Plate", "Well", "Gene", "Type", "Label", "Replicate", "Group", "group"])
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


def make_std_curve_figure(gene: str, x: np.ndarray, y: np.ndarray, slope: float, intercept: float, r2: float, eff: float):
    xp = np.linspace(float(np.min(x)), float(np.max(x)), 100)
    yp = slope * xp + intercept

    fig, ax = plt.subplots(figsize=(6.6, 4.0), dpi=140)
    fig.patch.set_facecolor(PLOT_BG)
    ax.set_facecolor(PLOT_BG)
    ax.grid(True, color=(1, 1, 1, 0.08), linewidth=0.8)

    ax.scatter(x, y, s=56, color="#38bdf8", edgecolor="#070b14", linewidth=0.6, label="Std mean Cq", zorder=3)
    ax.plot(xp, yp, linestyle="--", linewidth=2.2, color="#a78bfa", label=f"Fit (R²={r2:.3f}, Eff={eff:.1f}%)")

    ax.set_xlabel("log10(concentration)", color="#cbd5e1")
    ax.set_ylabel("Cq", color="#cbd5e1")
    ax.set_title(f"Standard curve — {gene}", color="#e8eefc", pad=10)

    for spine in ax.spines.values():
        spine.set_color((1, 1, 1, 0.18))
    ax.tick_params(colors="#cbd5e1")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    return fig


def fit_standard_curve(std_df: pd.DataFrame, mapping: pd.DataFrame):
    m = mapping.dropna(subset=["Label", "Concentration"]).copy()
    m["Label"] = m["Label"].astype(str)
    m["Concentration"] = m["Concentration"].astype(float)
    m["log10_conc"] = np.log10(m["Concentration"])
    df = std_df.merge(m, on="Label", how="inner")
    g = df.groupby(["Gene", "Label", "log10_conc"], as_index=False)["Cq"].mean().rename(columns={"Cq": "meanCq"})
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
    curves_df = pd.DataFrame(curves)
    if curves_df.empty:
        curves_df = pd.DataFrame(columns=["Gene", "slope", "intercept", "R2", "Efficiency_%", "n_points"])
    return curves_df, g


def quantify_samples(samples_df: pd.DataFrame, curves_df: pd.DataFrame):
    df = samples_df.merge(curves_df[["Gene", "slope", "intercept"]], on="Gene", how="left")
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
        sub.groupby(["Gene", "Type", "Label"], dropna=False)["Cq"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"mean": "AvgCq", "std": "SdCq", "count": "N"})
    )
    agg["CV%"] = (agg["SdCq"] / agg["AvgCq"]) * 100
    return agg


def normalize_to_ref(q_df: pd.DataFrame, ref_gene: str):
    """Return per-sample aggregates with normalization against ref gene."""
    if q_df.empty:
        cols = ["Gene", "Label", "Qty_mean", "Qty_sd", "n", "RefQty", "Norm_Qty"]
        return pd.DataFrame(columns=cols)
    agg = (
        q_df.groupby(["Gene", "Label"], dropna=False)["Quantity"]
        .agg(Qty_mean="mean", Qty_sd="std", n="count")
        .reset_index()
    )
    ref = (
        agg[agg["Gene"].str.lower() == str(ref_gene).lower()][["Label", "Qty_mean"]]
        .rename(columns={"Qty_mean": "RefQty"})
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


def _detect_sep(text: str) -> str:
    head = "\n".join(text.strip().splitlines()[:5])
    cand = {"\t": head.count("\t"), ",": head.count(","), ";": head.count(";"), "|": head.count("|")}
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
