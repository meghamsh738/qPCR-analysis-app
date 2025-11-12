import hashlib
import io
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Set

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st


NA_TOKENS = ["UNDETERMINED", "#DIV/0!", "NA", "NaN", "INF", "inf", "nan"]
CANDIDATE_SEPARATORS = [",", "\t", ";", "|"]
DEFAULT_REFERENCE_GENE = "gapdh"
NEGATIVE_LABEL_TOKENS = {"rt-", "rna-", "blank", "ntc", "water", "pcr-", "no template", "negative"}


@dataclass
class StandardCurve:
    """Container for per-gene standard curve attributes."""

    gene: str
    slope: float
    intercept: float
    r2: float
    efficiency: float
    n_points: int
    points: pd.DataFrame


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def detect_delimiter(text: str) -> Tuple[str, Dict[str, bool]]:
    """Detect a reasonable delimiter for pasted tabular data."""
    lines = [ln for ln in text.strip().splitlines() if ln.strip()]
    if not lines:
        return ",", {}
    best_sep, best_score = None, -1
    for sep in CANDIDATE_SEPARATORS:
        splits = [ln.count(sep) for ln in lines[:10]]
        score = np.mean(splits)
        if score > best_score:
            best_score, best_sep = score, sep
    if best_sep and best_score > 0:
        return best_sep, {}
    return r"\s+", {"engine": "python"}


def read_pasted_dataframe(name: str, raw_text: str, required: bool = True) -> Optional[pd.DataFrame]:
    """Parse pasted text into a DataFrame, warning if parsing fails."""
    if not raw_text.strip():
        if required:
            st.warning(f"Paste {name} data to continue.")
        return None
    sep, extra = detect_delimiter(raw_text)
    buffer = io.StringIO(raw_text.strip())
    try:
        if sep == r"\s+":
            df = pd.read_csv(buffer, sep=sep, engine="python", na_values=NA_TOKENS, keep_default_na=False)
        else:
            df = pd.read_csv(buffer, sep=sep, na_values=NA_TOKENS, keep_default_na=False, **extra)
    except Exception as exc:  # noqa: BLE001
        st.error(f"Could not parse {name} input: {exc}")
        return None
    return df.applymap(lambda x: x.strip() if isinstance(x, str) else x)


def rename_case_insensitive(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    """Align column names to a canonical case-insensitive set."""
    lower_map = {col.lower(): col for col in df.columns}
    rename_map: Dict[str, str] = {}
    for canonical in columns:
        lower = canonical.lower()
        actual = lower_map.get(lower)
        if actual and actual != canonical:
            rename_map[actual] = canonical
    if rename_map:
        df = df.rename(columns=rename_map)
    return df


def map_task_to_type(task_value: object) -> str:
    """Map vendor Task codes into simplified type labels."""
    if pd.isna(task_value):
        return "Sample"
    value = str(task_value).strip().lower()
    mapping = {
        "unknown": "Sample",
        "sample": "Sample",
        "unk": "Sample",
        "standard": "Standard",
        "std": "Standard",
        "calibrator": "Standard",
        "calibration": "Standard",
        "negative": "Negative",
        "ntc": "Negative",
        "blank": "Negative",
        "control": "Control",
        "reference": "Sample",
    }
    return mapping.get(value, "Sample")


def derive_type_from_label(label: object) -> str:
    """Infer assay type from the sample label when no vendor metadata is available."""
    if pd.isna(label):
        return "Sample"
    value = str(label).strip().lower().replace("−", "-")
    if value.startswith("std"):
        return "Standard"
    if any(token in value for token in NEGATIVE_LABEL_TOKENS):
        return "Negative"
    if "control" in value or "calibrator" in value:
        return "Control"
    return "Sample"


def transform_vendor_wells(df: pd.DataFrame) -> Optional[pd.DataFrame]:
    """Recognise common vendor exports and reshape into the canonical schema."""
    df = rename_case_insensitive(
        df,
        ["Plate", "Well", "Well Position", "Sample", "Target", "Task", "Cq", "Omit", "Amp Status"],
    )
    lower_map = {col.lower(): col for col in df.columns}
    if not {"target", "sample", "cq"}.issubset(lower_map.keys()):
        return None
    gene_col = lower_map["target"]
    sample_col = lower_map["sample"]
    cq_col = lower_map["cq"]
    well_col = lower_map.get("well") or lower_map.get("well position")
    plate_col = lower_map.get("plate")
    task_col = lower_map.get("task")
    omit_col = lower_map.get("omit")
    status_col = lower_map.get("amp status")

    base = pd.DataFrame()
    base["Gene"] = df[gene_col].astype(str).str.strip()
    base["Label"] = df[sample_col].astype(str).str.strip()
    base["Cq"] = pd.to_numeric(df[cq_col], errors="coerce")
    base["Plate"] = df[plate_col].astype(str).str.strip() if plate_col else "Plate 1"
    base["Well"] = df[well_col].astype(str).str.strip() if well_col else [f"Well_{idx+1}" for idx in range(len(df))]
    base["Type"] = df[task_col].apply(map_task_to_type) if task_col else base["Label"].apply(derive_type_from_label)
    base["Group"] = df[status_col].astype(str).str.strip() if status_col else ""
    base["Replicate"] = base.groupby(["Gene", "Label", "Type"]).cumcount() + 1
    if omit_col:
        base["Omit"] = df[omit_col].astype(str).str.strip().str.upper().isin({"TRUE", "1", "YES"})
    return base


def prepare_wells(df: pd.DataFrame) -> pd.DataFrame:
    """Standardise column names and derive missing metadata from the wells table."""
    df = rename_case_insensitive(
        df,
        ["Plate", "Well", "Gene", "Type", "Label", "Replicate", "Group", "Cq", "Concentration"],
    )
    vendor = transform_vendor_wells(df)
    if vendor is not None:
        df = vendor
    df = df.reset_index(drop=True)
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].astype(str).str.strip()
    if "Gene" not in df.columns or "Label" not in df.columns:
        st.error("Wells table requires at least 'Gene' and 'Label' columns.")
        return pd.DataFrame()
    if "Cq" not in df.columns:
        st.error("Could not identify a 'Cq' column.")
        return pd.DataFrame()
    df["Cq"] = pd.to_numeric(df["Cq"], errors="coerce")
    if "Plate" not in df.columns:
        df["Plate"] = "Plate 1"
    if "Well" not in df.columns:
        df["Well"] = [f"Well_{idx+1}" for idx in range(len(df))]
    if "Type" not in df.columns:
        df["Type"] = df["Label"].apply(derive_type_from_label)
    if "Group" not in df.columns:
        df["Group"] = ""
    if "Replicate" in df.columns:
        df["Replicate"] = pd.to_numeric(df["Replicate"], errors="coerce")
    else:
        df["Replicate"] = np.nan
    needs_rep = df["Replicate"].isna()
    if needs_rep.any():
        df.loc[needs_rep, "Replicate"] = (
            df.loc[needs_rep].groupby(["Gene", "Label", "Type"]).cumcount().add(1).values
        )
    df["Replicate"] = df["Replicate"].astype(int)
    df["InputOrder"] = np.arange(len(df))
    df["Type_lower"] = df["Type"].str.lower()
    df["IsNegativeControl"] = df["Type_lower"] == "negative"
    return df


# ---------------------------------------------------------------------------
# Aggregation helpers
# ---------------------------------------------------------------------------

def ensure_gene_state(genes: Sequence[str]) -> None:
    """Initialise Streamlit session state for keep/reference toggles."""
    if not genes:
        return
    lower_map = {g.lower(): g for g in genes}
    if "reference_gene" not in st.session_state or st.session_state["reference_gene"].lower() not in lower_map:
        default = lower_map.get(DEFAULT_REFERENCE_GENE.lower(), genes[0])
        st.session_state["reference_gene"] = default
    for gene in genes:
        keep_key = f"keep::{gene}"
        ref_key = f"ref::{gene}"
        if keep_key not in st.session_state:
            st.session_state[keep_key] = True
        if ref_key not in st.session_state:
            st.session_state[ref_key] = gene.lower() == st.session_state["reference_gene"].lower()
    ref_flags = [g for g in genes if st.session_state.get(f"ref::{g}", False)]
    if not ref_flags:
        fallback = st.session_state["reference_gene"]
        st.session_state[f"ref::{fallback}"] = True
        ref_flags = [fallback]
    reference = ref_flags[0]
    st.session_state["reference_gene"] = reference
    for other in ref_flags[1:]:
        st.session_state[f"ref::{other}"] = False


def get_reference_gene(genes: Sequence[str]) -> str:
    """Return the currently selected reference gene, ensuring it is valid."""
    ensure_gene_state(genes)
    reference = st.session_state.get("reference_gene", genes[0] if genes else "")
    if reference not in genes and genes:
        reference = genes[0]
        st.session_state["reference_gene"] = reference
    for gene in genes:
        st.session_state[f"ref::{gene}"] = gene == reference
    return reference


def get_kept_genes(genes: Sequence[str]) -> List[str]:
    """List genes that the user left checked for inclusion."""
    kept = [gene for gene in genes if st.session_state.get(f"keep::{gene}", True)]
    return kept


def reset_session_state_for_dataset(signature: str) -> None:
    """Reset per-dataset Streamlit session state when new data are pasted."""
    if not signature:
        return
    previous = st.session_state.get("dataset_signature")
    if previous == signature:
        return
    keys_to_clear = [
        key
        for key in list(st.session_state.keys())
        if key.startswith(("keep::", "ref::", "include::", "pending::"))
    ]
    for key in keys_to_clear:
        del st.session_state[key]
    for scalar in ["reference_gene", "selected_gene", "calibrator_labels"]:
        if scalar in st.session_state:
            del st.session_state[scalar]
    st.session_state["dataset_signature"] = signature


def compute_replicate_summary(wells: pd.DataFrame) -> pd.DataFrame:
    """Aggregate wells into replicate-level statistics."""
    if wells.empty:
        columns = [
            "Gene",
            "Type",
            "Label",
            "MeanCq",
            "StdCq",
            "Replicates",
            "Plates",
            "Groups",
            "IsNegativeControl",
            "FirstOrder",
            "CqSpread",
            "DeltaRep12",
            "FlagLargeDelta",
        ]
        return pd.DataFrame(columns=columns)
    group_cols = ["Gene", "Type", "Label"]
    wells_sorted = wells.sort_values("InputOrder")
    agg = (
        wells_sorted.groupby(group_cols, dropna=False, sort=False)
        .agg(
            MeanCq=("Cq", "mean"),
            StdCq=("Cq", lambda series: series.std(ddof=1) if series.count() > 1 else np.nan),
            Replicates=("Cq", "count"),
            Plates=("Plate", lambda series: ", ".join(sorted(series.dropna().astype(str).unique()))),
            Groups=("Group", lambda series: ", ".join(sorted(g for g in series.dropna().astype(str).unique() if g))),
            IsNegativeControl=("IsNegativeControl", "max"),
            FirstOrder=("InputOrder", "min"),
            CqMin=("Cq", "min"),
            CqMax=("Cq", "max"),
        )
        .reset_index()
    )
    delta_records: List[Dict[str, object]] = []
    for (gene, sample_type, label), group in wells_sorted.groupby(group_cols, sort=False):
        ordered = group.sort_values(["Replicate", "InputOrder"])
        if len(ordered) >= 2:
            delta = abs(float(ordered.iloc[0]["Cq"]) - float(ordered.iloc[1]["Cq"]))
        else:
            delta = np.nan
        delta_records.append(
            {
                "Gene": gene,
                "Type": sample_type,
                "Label": label,
                "DeltaRep12": delta,
            }
        )
    delta_df = pd.DataFrame(delta_records)
    agg = agg.merge(delta_df, on=group_cols, how="left")
    agg["MeanCq"] = agg["MeanCq"].round(3)
    agg["StdCq"] = agg["StdCq"].round(3)
    agg["CqSpread"] = (agg["CqMax"] - agg["CqMin"]).round(3)
    agg["DeltaRep12"] = agg["DeltaRep12"].round(3)
    agg["IsNegativeControl"] = agg["IsNegativeControl"].astype(bool)
    agg["FlagLargeDelta"] = agg["DeltaRep12"] > 1
    agg = agg.drop(columns=["CqMin", "CqMax"])
    agg = agg.sort_values(["Gene", "FirstOrder"]).reset_index(drop=True)
    return agg


def flag_replicate_variance(replicates: pd.DataFrame) -> pd.DataFrame:
    """Flag replicate groups whose SD lies above the 75th percentile for each gene."""
    if replicates.empty:
        replicates["HighSdThreshold"] = np.nan
        replicates["FlagHighVariance"] = False
        return replicates
    flagged = replicates.copy()
    flagged["HighSdThreshold"] = np.nan
    flagged["FlagHighVariance"] = False
    for gene, gene_df in flagged.groupby("Gene"):
        sd_values = gene_df["StdCq"].dropna()
        if sd_values.empty:
            continue
        threshold = float(sd_values.quantile(0.75))
        idx = gene_df.index
        flagged.loc[idx, "HighSdThreshold"] = threshold
        high_idx = gene_df[gene_df["StdCq"] > threshold].index
        flagged.loc[high_idx, "FlagHighVariance"] = True
    flagged["HighSdThreshold"] = flagged["HighSdThreshold"].round(3)
    return flagged


def extract_standard_concentrations(
    df: pd.DataFrame,
) -> Tuple[Dict[str, float], Dict[str, Dict[str, float]], pd.DataFrame, List[str]]:
    """Collect concentration information for standard wells."""
    standards_columns = ["Gene", "Label", "Concentration", "Source", "LabelOrder"]
    empty_standards = pd.DataFrame(columns=standards_columns)
    if df.empty:
        return {}, {}, empty_standards, []
    conc_col = next((col for col in df.columns if col.lower() in {"concentration", "conc"}), None)
    if conc_col:
        df = df.copy()
        df[conc_col] = pd.to_numeric(df[conc_col], errors="coerce")
    warnings: List[str] = []
    per_gene_map: Dict[str, Dict[str, float]] = {}
    global_map: Dict[str, float] = {}
    records: List[Dict[str, object]] = []
    standards = df[df["Type_lower"] == "standard"]
    if standards.empty:
        return global_map, per_gene_map, empty_standards, ["No standards detected in the wells table."]
    default_per_gene: Dict[str, Dict[str, float]] = {}
    for gene, group in standards.groupby("Gene", sort=False):
        order_series = (
            group.groupby("Label")["InputOrder"].min().sort_values()
        )
        defaults: Dict[str, float] = {}
        for idx, (label, _) in enumerate(order_series.items()):
            defaults[label] = 1000.0 / (10 ** idx)
        default_per_gene[gene.lower()] = defaults
    for (gene, label), group in standards.groupby(["Gene", "Label"], sort=False):
        label_order = int(group["InputOrder"].min())
        if conc_col:
            values = group[conc_col].dropna().unique()
        else:
            values = np.array([])
        if len(values) == 0:
            fallback = default_per_gene.get(gene.lower(), {}).get(label)
            if fallback is None:
                warnings.append(f"{gene} missing concentration for label '{label}'.")
                records.append(
                    {
                        "Gene": gene,
                        "Label": label,
                        "Concentration": np.nan,
                        "Source": "Missing",
                        "LabelOrder": label_order,
                    }
                )
                continue
            value = float(fallback)
            source = "Default 10-fold series (1000→0.01)"
        else:
            value = float(values[0])
            if len(values) > 1:
                warnings.append(f"{gene} label '{label}' has multiple concentration values; using {value}.")
            source = "Provided"
        per_gene_map.setdefault(gene.lower(), {})[label] = value
        global_map.setdefault(label, value)
        records.append(
            {
                "Gene": gene,
                "Label": label,
                "Concentration": value,
                "Source": source,
                "LabelOrder": label_order,
            }
        )
    standards_sheet = pd.DataFrame(records, columns=standards_columns)
    if not standards_sheet.empty:
        standards_sheet = standards_sheet.sort_values(["Gene", "LabelOrder"]).reset_index(drop=True)
    return global_map, per_gene_map, standards_sheet, warnings


def fit_standard_curves(
    standards_df: pd.DataFrame,
    global_map: Dict[str, float],
    per_gene_map: Dict[str, Dict[str, float]],
) -> Tuple[Dict[str, StandardCurve], pd.DataFrame, pd.DataFrame, List[str]]:
    """Fit a linear regression for each gene's standard points."""
    curves: Dict[str, StandardCurve] = {}
    points_records: List[Dict[str, object]] = []
    missing_records: List[Dict[str, object]] = []
    warnings: List[str] = []
    if standards_df.empty:
        return curves, pd.DataFrame(), pd.DataFrame(), warnings
    for gene, group in standards_df.groupby("Gene", sort=False):
        records: List[Dict[str, object]] = []
        missing_labels: List[str] = []
        for _, row in group.iterrows():
            label = row["Label"]
            mean_cq = row.get("MeanCq")
            if pd.isna(mean_cq):
                continue
            conc = per_gene_map.get(gene.lower(), {}).get(label, global_map.get(label))
            if conc is None:
                missing_labels.append(label)
                continue
            if conc <= 0:
                warnings.append(f"{gene} label '{label}' has non-positive concentration; skipped.")
                continue
            records.append(
                {
                    "Gene": gene,
                    "Label": label,
                    "Concentration": conc,
                    "log10_conc": np.log10(conc),
                    "MeanCq": mean_cq,
                }
            )
        if missing_labels:
            warnings.append(f"{gene} missing concentrations for labels: {', '.join(sorted(missing_labels))}.")
        records_df = pd.DataFrame(records)
        if records_df.empty or records_df["log10_conc"].isna().any():
            missing_records.append({"Gene": gene, "Reason": "Missing concentrations or invalid values", "Levels": 0})
            continue
        if records_df["Label"].nunique() < 2:
            missing_records.append(
                {"Gene": gene, "Reason": "Less than 2 standard levels", "Levels": records_df["Label"].nunique()}
            )
            continue
        x = records_df["log10_conc"].values
        y = records_df["MeanCq"].values
        slope, intercept = np.polyfit(x, y, 1)
        predicted = slope * x + intercept
        ss_res = float(np.sum((y - predicted) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        efficiency = (10 ** (-1 / slope) - 1) * 100 if slope != 0 else np.nan
        curves[gene] = StandardCurve(
            gene=gene,
            slope=float(slope),
            intercept=float(intercept),
            r2=float(r2),
            efficiency=float(efficiency),
            n_points=len(records_df),
            points=records_df,
        )
        points_records.extend(records_df.to_dict(orient="records"))
    return curves, pd.DataFrame(points_records), pd.DataFrame(missing_records), warnings


def compute_sample_quantities(samples_df: pd.DataFrame, curves: Dict[str, StandardCurve]) -> pd.DataFrame:
    """Interpolate replicate means against fitted standard curves."""
    if samples_df.empty:
        return pd.DataFrame(columns=list(samples_df.columns) + ["MeanQuantity"])
    samples_sorted = samples_df.sort_values("FirstOrder")
    records: List[Dict[str, object]] = []
    for _, row in samples_sorted.iterrows():
        gene = row["Gene"]
        mean_cq = row["MeanCq"]
        curve = curves.get(gene)
        quantity = np.nan
        if curve and not np.isnan(mean_cq) and curve.slope != 0:
            log10_q = (mean_cq - curve.intercept) / curve.slope
            if not np.isnan(log10_q):
                quantity = 10 ** log10_q
        rec = row.to_dict()
        rec["MeanQuantity"] = quantity
        records.append(rec)
    result = pd.DataFrame(records)
    return result.sort_values("FirstOrder").reset_index(drop=True)


def normalise_to_reference(
    sample_quantities: pd.DataFrame,
    reference_gene: str,
) -> Tuple[pd.DataFrame, Optional[str]]:
    """Normalise sample quantities against the chosen reference gene."""
    if sample_quantities.empty:
        return pd.DataFrame(), "No sample quantities available for normalisation."
    ref_rows = sample_quantities[sample_quantities["Gene"].str.lower() == reference_gene.lower()]
    if ref_rows.empty or ref_rows["MeanQuantity"].isna().all():
        return pd.DataFrame(), f"Reference gene '{reference_gene}' has no quantities."
    ref_map = ref_rows.set_index("Label")["MeanQuantity"].to_dict()
    records: List[Dict[str, object]] = []
    for _, row in sample_quantities.sort_values("FirstOrder").iterrows():
        label = row["Label"]
        ref_quantity = ref_map.get(label)
        target_quantity = row["MeanQuantity"]
        if ref_quantity is None or np.isnan(ref_quantity) or np.isnan(target_quantity):
            normalised = np.nan
        elif ref_quantity == 0:
            normalised = np.nan
        else:
            normalised = target_quantity / ref_quantity
        rec = row.to_dict()
        rec["ReferenceGene"] = reference_gene
        rec["ReferenceMeanQuantity"] = ref_quantity
        rec["NormalisedQuantity"] = normalised
        records.append(rec)
    normalised_df = pd.DataFrame(records)
    if normalised_df.empty:
        return normalised_df, None
    normalised_df = normalised_df.sort_values(["Gene", "FirstOrder"]).reset_index(drop=True)
    return normalised_df, None


def compute_relative_expression_all(
    replicate_summary: pd.DataFrame,
    reference_gene: str,
    calibrator_labels: Sequence[str],
) -> Tuple[pd.DataFrame, Optional[str]]:
    """Calculate 2^-ΔΔCt values for every gene in the replicate summary."""
    if replicate_summary.empty:
        return pd.DataFrame(), "No replicate averages available."
    non_standard = replicate_summary[replicate_summary["Type"].str.lower() != "standard"]
    if non_standard.empty:
        return pd.DataFrame(), "No non-standard samples available for relative expression."
    ref_rows = non_standard[non_standard["Gene"].str.lower() == reference_gene.lower()]
    if ref_rows.empty:
        return pd.DataFrame(), f"Reference gene '{reference_gene}' missing replicate averages."
    ref_map = ref_rows.set_index("Label")["MeanCq"].to_dict()
    records: List[pd.DataFrame] = []
    warnings: List[str] = []
    for gene, group in non_standard.groupby("Gene"):
        gene_records: List[Dict[str, object]] = []
        group_sorted = group.sort_values("FirstOrder")
        for _, row in group_sorted.iterrows():
            label = row["Label"]
            target_mean = row["MeanCq"]
            ref_mean = ref_map.get(label)
            delta_ct = target_mean - ref_mean if ref_mean is not None and not np.isnan(ref_mean) else np.nan
            gene_records.append(
                {
                    "Gene": gene,
                    "Label": label,
                    "TargetMeanCq": target_mean,
                    "ReferenceMeanCq": ref_mean,
                    "DeltaCt": delta_ct,
                    "FirstOrder": row.get("FirstOrder"),
                }
            )
        if not gene_records:
            continue
        gene_df = pd.DataFrame(gene_records)
        calibrator_subset = gene_df[gene_df["Label"].isin(calibrator_labels)]
        calibrator_subset = calibrator_subset[calibrator_subset["DeltaCt"].notna()]
        if calibrator_labels and calibrator_subset.empty:
            warnings.append(f"{gene}: calibrator labels missing or lacking ΔCt values.")
            calibrator_mean = np.nan
        elif calibrator_subset.empty:
            calibrator_mean = np.nan
        else:
            calibrator_mean = calibrator_subset["DeltaCt"].mean()
        gene_df["CalibratorMeanDeltaCt"] = calibrator_mean
        if np.isnan(calibrator_mean):
            gene_df["DeltaDeltaCt"] = np.nan
            gene_df["RelativeExpression"] = np.nan
        else:
            gene_df["DeltaDeltaCt"] = gene_df["DeltaCt"] - calibrator_mean
            gene_df["RelativeExpression"] = np.power(2.0, -gene_df["DeltaDeltaCt"])
        records.append(gene_df)
    result = pd.concat(records, ignore_index=True) if records else pd.DataFrame()
    if not result.empty and "FirstOrder" in result.columns:
        result = result.sort_values(["Gene", "FirstOrder"]).reset_index(drop=True)
    combined_warning = "; ".join(warnings) if warnings else None
    return result, combined_warning


def sanitise_sheet_name(name: str) -> str:
    """Trim sheet names to Excel-safe representations."""
    invalid = set(r"[]:*?/\\")
    cleaned = "".join(ch for ch in name if ch not in invalid).strip()
    if not cleaned:
        cleaned = "Sheet"
    return cleaned[:31]


def build_excel_stream(
    cleaned_wells: pd.DataFrame,
    standards_map_df: pd.DataFrame,
    replicate_averages: pd.DataFrame,
    std_curve_points: pd.DataFrame,
    curves: Dict[str, StandardCurve],
    sample_quantities: pd.DataFrame,
    normalised_df: pd.DataFrame,
    relative_df: pd.DataFrame,
    missing_standards: pd.DataFrame,
) -> io.BytesIO:
    """Create an Excel workbook containing all analysis tables."""
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="xlsxwriter") as writer:
        cleaned_wells.to_excel(writer, sheet_name="Cleaned_Wells", index=False)
        for gene, group in cleaned_wells.groupby("Gene"):
            sheet = sanitise_sheet_name(f"Gene_{gene}")
            group.drop(columns=["Type_lower"], errors="ignore").to_excel(writer, sheet_name=sheet, index=False)
        replicate_averages.to_excel(writer, sheet_name="Replicate_Averages", index=False)
        standards_map_df.to_excel(writer, sheet_name="Standards_Map", index=False)
        std_curve_points.to_excel(writer, sheet_name="StdCurve_Points", index=False)
        fits_records = [
            {
                "Gene": curve.gene,
                "Slope": round(curve.slope, 6),
                "Intercept": round(curve.intercept, 6),
                "R2": round(curve.r2, 4) if not np.isnan(curve.r2) else np.nan,
                "EfficiencyPct": round(curve.efficiency, 2) if not np.isnan(curve.efficiency) else np.nan,
                "PointsUsed": curve.n_points,
            }
            for curve in curves.values()
        ]
        fits_df = pd.DataFrame(fits_records).sort_values("Gene") if fits_records else pd.DataFrame(
            columns=["Gene", "Slope", "Intercept", "R2", "EfficiencyPct", "PointsUsed"]
        )
        fits_df.to_excel(writer, sheet_name="StdCurve_Fits", index=False)
        sample_quantities.to_excel(writer, sheet_name="Sample_Quantities", index=False)
        normalised_df.to_excel(writer, sheet_name="Normalised", index=False)
        relative_df.to_excel(writer, sheet_name="Relative_Expression", index=False)
        missing_standards.to_excel(writer, sheet_name="Missing_Standards", index=False)
        workbook = writer.book
        plots_sheet = workbook.add_worksheet("StdCurve_Plots")
        writer.sheets["StdCurve_Plots"] = plots_sheet
        plots_sheet.write(0, 0, "Standard curve plots for each gene")
        plots_sheet.set_column(0, 0, 2)
        plots_sheet.set_column(1, 5, 2)
        plots_sheet.set_column(6, 8, 45)
        row_offset = 2
        image_buffers: List[io.BytesIO] = []
        for curve in curves.values():
            fig = render_standard_curve_plot(curve)
            img_buffer = io.BytesIO()
            fig.savefig(img_buffer, format="png", dpi=180, bbox_inches="tight")
            plt.close(fig)
            img_buffer.seek(0)
            image_buffers.append(img_buffer)
            plots_sheet.insert_image(
                row_offset,
                1,
                f"{curve.gene}.png",
                {"image_data": img_buffer, "x_scale": 0.9, "y_scale": 0.9},
            )
            plots_sheet.write(
                row_offset,
                6,
                f"Cq = {curve.slope:.3f}·log10(conc) + {curve.intercept:.3f}",
            )
            plots_sheet.write(
                row_offset + 1,
                6,
                f"Quantity = 10^((Cq - {curve.intercept:.3f})/{curve.slope:.3f})",
            )
            plots_sheet.write(row_offset + 2, 6, f"R² = {curve.r2:.3f} | Efficiency = {curve.efficiency:.1f}%")
            row_offset += 28
    buffer.seek(0)
    return buffer


def style_high_variance(row: pd.Series) -> List[str]:
    """Helper for styling replicate tables."""
    variance_flag = bool(row.get("FlagHighVariance") or row.get("High SD flag"))
    delta_flag = bool(row.get("FlagLargeDelta") or row.get("Δ>1 flag"))
    if variance_flag:
        color = "background-color: #fde4e4;"
    elif delta_flag:
        color = "background-color: #fff4ce;"
    else:
        color = ""
    return [color] * len(row)


def render_standard_curve_plot(curve: StandardCurve) -> plt.Figure:
    """Generate a matplotlib plot for a standard curve."""
    points = curve.points.sort_values("log10_conc")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(points["log10_conc"], points["MeanCq"], color="#1f77b4", label="Standards")
    x_vals = np.linspace(points["log10_conc"].min(), points["log10_conc"].max(), 100)
    y_vals = curve.slope * x_vals + curve.intercept
    ax.plot(x_vals, y_vals, color="#ff7f0e", label="Regression")
    ax.set_xlabel("log10(Standard concentration)")
    ax.set_ylabel("Mean Cq")
    ax.set_title(curve.gene)
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend()
    quantity_eq = r"Quantity = 10$^{(Cq - %.3f)/%.3f}$" % (curve.intercept, curve.slope)
    regression_eq = r"Cq = %.3f·log10(conc) + %.3f" % (curve.slope, curve.intercept)
    ax.text(
        0.02,
        0.05,
        f"{regression_eq}\n{quantity_eq}",
        transform=ax.transAxes,
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.7, "pad": 4},
    )
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Streamlit interface
# ---------------------------------------------------------------------------

def main() -> None:
    st.set_page_config(page_title="qPCR replicate analysis", layout="wide")
    st.title("qPCR replicate analysis dashboard")
    st.caption(
        "Paste your wells table, choose the gene of interest, and step through each tab to review replicate checks, "
        "standard curves, quantities, normalisation, and relative expression. All tables are synchronised to the gene selector."
    )

    with st.sidebar:
        st.header("Step 0 – Paste wells")
        st.write(
            "Copy the full wells table from your qPCR instrument (including standards) and paste it below. "
            "CSV, TSV, and Excel-style tables are detected automatically."
        )
        wells_text = st.text_area("Wells table", height=240, key="wells_text")
        st.caption(
            "Required columns: Gene, Label, Cq. Optional columns such as Plate, Well, Type, Replicate, Concentration, "
            "and Amp Status will be used when available."
        )

    if not wells_text or not wells_text.strip():
        st.info("Paste the wells table in the sidebar to begin the analysis.")
        return

    dataset_signature = hashlib.sha256(wells_text.strip().encode("utf-8")).hexdigest()
    reset_session_state_for_dataset(dataset_signature)

    raw_df = read_pasted_dataframe("wells", wells_text)
    if raw_df is None or raw_df.empty:
        st.stop()

    wells = prepare_wells(raw_df)
    if wells.empty:
        st.stop()
    wells["StateKey"] = wells["InputOrder"].apply(lambda idx: f"{dataset_signature}:{int(idx)}")

    genes = sorted(wells["Gene"].unique())
    ensure_gene_state(genes)
    reference_gene = get_reference_gene(genes)
    keep_genes = get_kept_genes(genes)
    if reference_gene not in keep_genes:
        st.warning(
            f"Reference gene '{reference_gene}' was unchecked for analysis. It has been re-enabled so normalisation can run."
        )
        st.session_state[f"keep::{reference_gene}"] = True
        keep_genes = get_kept_genes(genes)

    if not keep_genes:
        st.error("All genes are currently deselected. Re-enable at least one gene using the controls in Step 1.")
        return

    filtered_wells = wells[wells["Gene"].isin(keep_genes)].copy()
    filtered_wells["KeepGene"] = filtered_wells["Gene"].apply(lambda g: st.session_state.get(f"keep::{g}", True))
    filtered_wells["IsReferenceGene"] = filtered_wells["Gene"].str.lower() == reference_gene.lower()
    for state_key in filtered_wells["StateKey"]:
        include_key = f"include::{state_key}"
        if include_key not in st.session_state:
            st.session_state[include_key] = True
    filtered_wells["Include"] = filtered_wells["StateKey"].apply(
        lambda key: st.session_state.get(f"include::{key}", True)
    )
    active_wells = filtered_wells[filtered_wells["Include"]].copy()

    replicate_summary = compute_replicate_summary(active_wells)
    replicate_summary = flag_replicate_variance(replicate_summary)

    global_map, per_gene_map, standards_map_df, standards_warnings = extract_standard_concentrations(active_wells)
    standards_replicates = replicate_summary[replicate_summary["Type"].str.lower() == "standard"].copy()
    standards_replicates = standards_replicates.sort_values("FirstOrder")
    curves, std_curve_points, missing_standards_df, curve_warnings = fit_standard_curves(
        standards_replicates, global_map, per_gene_map
    )

    sample_replicates = replicate_summary[replicate_summary["Type"].str.lower() != "standard"].copy()
    sample_replicates = sample_replicates.sort_values("FirstOrder")
    sample_quantities = compute_sample_quantities(sample_replicates, curves)
    normalised_df, normalisation_warning = normalise_to_reference(sample_quantities, reference_gene)

    # Default calibrators: Control or Calibrator types, otherwise first label.
    control_labels: List[str] = []
    seen_controls: Set[str] = set()
    for label in sample_replicates[sample_replicates["Type"].str.lower().isin({"control"})]["Label"]:
        if label not in seen_controls:
            control_labels.append(label)
            seen_controls.add(label)
    all_labels: List[str] = []
    seen_all: Set[str] = set()
    for label in sample_replicates["Label"]:
        if label not in seen_all:
            all_labels.append(label)
            seen_all.add(label)
    default_calibrators = control_labels if control_labels else (all_labels[:1] if all_labels else [])
    if "calibrator_labels" not in st.session_state:
        st.session_state["calibrator_labels"] = default_calibrators
    calibrator_labels = st.session_state.get("calibrator_labels", [])

    relative_all, relative_warning = compute_relative_expression_all(
        replicate_summary, reference_gene, calibrator_labels
    )

    tabs = st.tabs(
        [
            "Step 1 – Selection",
            "Step 2 – Standard curve",
            "Step 3 – Quantities",
            "Step 4 – Normalisation",
            "Step 5 – Relative expression",
            "Step 6 – Missing standards",
        ]
    )

    with tabs[0]:
        st.subheader("Step 1 – Selection")
        st.write(
            "Use the gene-level checkboxes to keep or exclude genes from the downstream analysis. Mark exactly one gene as "
            "the reference (feeds normalisation and ΔΔCt). The editable table lets you untick individual wells/replicates "
            "before they influence any calculations. Toggle as needed, then click **Apply replicate selections** to update the analysis."
        )
        for gene in genes:
            cols = st.columns([3, 1, 1])
            cols[0].markdown(f"**{gene}**")
            cols[1].checkbox("Keep", key=f"keep::{gene}")
            cols[2].checkbox("Reference", key=f"ref::{gene}")
        # Enforce a single reference gene flag.
        reference_flags = [g for g in genes if st.session_state.get(f"ref::{g}", False)]
        if not reference_flags:
            st.session_state[f"ref::{reference_gene}"] = True
            reference_flags = [reference_gene]
        primary_ref = reference_flags[0]
        st.session_state["reference_gene"] = primary_ref
        for other in reference_flags[1:]:
            st.session_state[f"ref::{other}"] = False
        keep_genes = get_kept_genes(genes)
        if primary_ref not in keep_genes:
            st.error("The reference gene must remain in the analysis. Re-enable it before proceeding.")
            keep_genes.append(primary_ref)
            st.session_state[f"keep::{primary_ref}"] = True
        selectable_genes = sorted(set(keep_genes))
        if "selected_gene" not in st.session_state or st.session_state["selected_gene"] not in selectable_genes:
            st.session_state["selected_gene"] = selectable_genes[0]
        selected_gene = st.selectbox("Gene to display", selectable_genes, key="selected_gene")
        gene_wells = filtered_wells[filtered_wells["Gene"] == selected_gene].copy()
        gene_wells = gene_wells.sort_values("InputOrder")
        if gene_wells.empty:
            st.info("No wells are currently available for this gene.")
        else:
            delta_lookup: Dict[Tuple[str, str], float] = {}
            if not replicate_summary.empty:
                delta_lookup = (
                    replicate_summary[replicate_summary["Gene"] == selected_gene]
                    .set_index(["Type", "Label"])["DeltaRep12"]
                    .to_dict()
                )
            gene_wells["DeltaRep12"] = gene_wells.apply(
                lambda row: delta_lookup.get((row["Type"], row["Label"])), axis=1
            )
            display_df = gene_wells[
                ["Include", "Type", "Label", "Replicate", "Cq", "DeltaRep12", "Plate", "Well"]
            ].rename(
                columns={
                    "Include": "Use",
                    "DeltaRep12": "|Δ Rep1-Rep2|",
                }
            )
            display_df.index = gene_wells["StateKey"].tolist()
            display_df.index.name = "RowKey"
            display_df["Use"] = display_df["Use"].astype(bool)
            display_df["Cq"] = display_df["Cq"].round(3)
            display_df["|Δ Rep1-Rep2|"] = display_df["|Δ Rep1-Rep2|"].round(3)
            warning_flags = display_df["|Δ Rep1-Rep2|"].abs() > 1
            display_df["Warning"] = ""
            display_df.loc[warning_flags, "Warning"] = "⚠️ Δ>1"
            column_config = {
                "Use": st.column_config.CheckboxColumn(
                    "Use",
                    help="Untick to exclude this replicate or standard from downstream analysis.",
                ),
                "|Δ Rep1-Rep2|": st.column_config.NumberColumn("|Δ Rep1-Rep2|", format="%.3f", disabled=True),
                "Warning": st.column_config.TextColumn(
                    "Warning",
                    help="⚠️ flags replicate pairs where |Δ Rep1-Rep2| exceeds 1 Cq.",
                    disabled=True,
                ),
                "Type": st.column_config.TextColumn("Type", disabled=True),
                "Label": st.column_config.TextColumn("Label", disabled=True),
                "Replicate": st.column_config.NumberColumn("Replicate", disabled=True, format="%.0f"),
                "Cq": st.column_config.NumberColumn("Cq", disabled=True, format="%.3f"),
                "Plate": st.column_config.TextColumn("Plate", disabled=True),
                "Well": st.column_config.TextColumn("Well", disabled=True),
            }
            with st.form(key=f"form::{selected_gene}"):
                edited = st.data_editor(
                    display_df,
                    hide_index=True,
                    column_config=column_config,
                    disabled=["Type", "Label", "Replicate", "Cq", "|Δ Rep1-Rep2|", "Warning", "Plate", "Well"],
                    key=f"editor::{selected_gene}",
                    use_container_width=True,
                    height=420,
                )
                st.caption(
                    "Toggle the Use column to include or exclude individual replicates and standards. "
                    "Rows marked with ⚠️ have |Δ Rep1-Rep2| > 1 Cq so both replicates merit review."
                )
                submitted = st.form_submit_button("Apply replicate selections")
            if submitted:
                for row_key, use_flag in edited["Use"].items():
                    st.session_state[f"include::{row_key}"] = bool(use_flag)
                st.rerun()

            selected_replicates = replicate_summary[replicate_summary["Gene"] == selected_gene].copy()
            selected_replicates = selected_replicates.sort_values("FirstOrder")
            if not selected_replicates.empty:
                summary_display = selected_replicates.rename(
                    columns={
                        "MeanCq": "Average Cq",
                        "StdCq": "Std Cq",
                        "Replicates": "n",
                        "IsNegativeControl": "Negative control",
                        "DeltaRep12": "|Δ Rep1-Rep2|",
                        "CqSpread": "Cq range",
                        "HighSdThreshold": "75th pct Std Cq",
                        "FlagHighVariance": "High SD flag",
                        "FlagLargeDelta": "Δ>1 flag",
                    }
                )
                summary_display = summary_display[
                    [
                        "Type",
                        "Label",
                        "Average Cq",
                        "Std Cq",
                        "|Δ Rep1-Rep2|",
                        "Cq range",
                        "n",
                        "Plates",
                        "Groups",
                        "Negative control",
                        "75th pct Std Cq",
                        "High SD flag",
                        "Δ>1 flag",
                    ]
                ]
                styled = summary_display.style.apply(style_high_variance, axis=1)
                with st.expander("Replicate summary", expanded=False):
                    st.dataframe(styled, use_container_width=True, height=240)
                    st.caption(
                        "Average Cq is the replicate mean. Std Cq is the standard deviation across wells in the Type × Label group; "
                        "use it to spot noisy replicates alongside the |Δ Rep1-Rep2| metric."
                    )

    with tabs[1]:
        st.subheader("Step 2 – Standard curve")
        gene = st.session_state["selected_gene"]
        gene_standards = standards_map_df[standards_map_df["Gene"] == gene]
        if not gene_standards.empty:
            if "LabelOrder" in gene_standards.columns:
                gene_standards = gene_standards.sort_values("LabelOrder")
            display = gene_standards.drop(columns=["LabelOrder"], errors="ignore")
            st.dataframe(display, use_container_width=True, height=220)
            st.caption(
                "Standard labels detected for this gene along with the concentration used. The Source column notes whether "
                "the value came from the wells table or was filled by a 10-fold default (1000→0.01 series)."
            )
        else:
            st.info("No standards detected for this gene.")
        if gene in curves:
            curve = curves[gene]
            fig = render_standard_curve_plot(curve)
            st.pyplot(fig, use_container_width=True)
            st.caption(
                f"Cq = {curve.slope:.3f}·log10(conc) + {curve.intercept:.3f} (R² = {curve.r2:.3f}). "
                f"Amplification efficiency = {curve.efficiency:.1f}%. "
                "Healthy qPCR assays typically sit between 90% and 110% efficiency."
            )
        else:
            st.warning("A standard curve could not be fitted for this gene.")
        for warning_msg in standards_warnings + curve_warnings:
            st.warning(warning_msg)

    with tabs[2]:
        st.subheader("Step 3 – Quantities")
        gene = st.session_state["selected_gene"]
        gene_quant = sample_quantities[sample_quantities["Gene"] == gene].copy()
        if not gene_quant.empty and "FirstOrder" in gene_quant.columns:
            gene_quant = gene_quant.sort_values("FirstOrder")
        if not gene_quant.empty:
            gene_quant["MeanQuantity"] = gene_quant["MeanQuantity"].round(4)
            display_cols = ["Type", "Label", "MeanCq", "StdCq", "Replicates", "MeanQuantity"]
            rename_map = {
                "MeanCq": "Average Cq",
                "StdCq": "Std Cq",
                "Replicates": "n",
                "MeanQuantity": "Mean Quantity",
            }
            st.dataframe(gene_quant[display_cols].rename(columns=rename_map), use_container_width=True, height=320)
            st.caption(
                "Mean quantities are interpolated from the standard curve using Mean Cq. Values are in the same units as "
                "the supplied standards."
            )
        else:
            st.info("No sample quantities could be calculated for this gene.")

    with tabs[3]:
        st.subheader("Step 4 – Normalisation")
        gene = st.session_state["selected_gene"]
        if normalisation_warning:
            st.warning(normalisation_warning)
        gene_norm = pd.DataFrame()
        if not normalised_df.empty:
            gene_norm = normalised_df[normalised_df["Gene"] == gene].copy()
        if not gene_norm.empty and "FirstOrder" in gene_norm.columns:
            gene_norm = gene_norm.sort_values("FirstOrder")
        if not gene_norm.empty:
            gene_norm["NormalisedQuantity"] = gene_norm["NormalisedQuantity"].round(4)
            gene_norm["MeanQuantity"] = gene_norm["MeanQuantity"].round(4)
            gene_norm["ReferenceMeanQuantity"] = gene_norm["ReferenceMeanQuantity"].round(4)
            display_cols = ["Label", "MeanQuantity", "ReferenceMeanQuantity", "NormalisedQuantity", "Type"]
            if gene.lower() == reference_gene.lower():
                target_quantity_label = f"Mean Quantity ({gene}) – target"
                reference_quantity_label = f"Mean Quantity ({reference_gene}) – reference"
            else:
                target_quantity_label = f"Mean Quantity ({gene})"
                reference_quantity_label = f"Mean Quantity ({reference_gene})"
            rename_map = {
                "MeanQuantity": target_quantity_label,
                "ReferenceMeanQuantity": reference_quantity_label,
                "NormalisedQuantity": "Normalised Quantity",
            }
            st.dataframe(gene_norm[display_cols].rename(columns=rename_map), use_container_width=True, height=320)
            st.caption(
                f"Normalised Quantity = Mean Quantity ({gene}) ÷ Mean Quantity ({reference_gene}). "
                "Rows with missing reference values show NaN."
            )
        else:
            st.info("Normalised quantities are unavailable. Check that the reference gene has valid quantities.")

    with tabs[4]:
        st.subheader("Step 5 – Relative expression (2^-ΔΔCt)")
        gene = st.session_state["selected_gene"]
        gene_samples = sample_replicates[sample_replicates["Gene"] == gene].sort_values("FirstOrder")
        ordered_labels = gene_samples["Label"].tolist()
        options: List[str] = []
        seen_labels: Set[str] = set()
        for label in ordered_labels:
            if label not in seen_labels:
                options.append(label)
                seen_labels.add(label)
        # Ensure session state only stores available labels.
        st.session_state["calibrator_labels"] = [label for label in calibrator_labels if label in options]
        calibrator_selection = st.multiselect(
            "Calibrator labels",
            options=options,
            default=st.session_state["calibrator_labels"],
            key="calibrator_labels",
            help="Pick control samples that represent your baseline condition. The ΔΔCt is referenced to the mean ΔCt of this pool.",
        )
        calibrator_labels = st.session_state.get("calibrator_labels", [])
        relative_all, relative_warning = compute_relative_expression_all(
            replicate_summary, reference_gene, calibrator_labels
        )
        relative_gene = relative_all[relative_all["Gene"] == gene].copy() if not relative_all.empty else pd.DataFrame()
        if not relative_gene.empty and "FirstOrder" in relative_gene.columns:
            relative_gene = relative_gene.sort_values("FirstOrder")
        if not relative_gene.empty:
            relative_gene["RelativeExpression"] = relative_gene["RelativeExpression"].round(4)
            display_cols = [
                "Label",
                "TargetMeanCq",
                "ReferenceMeanCq",
                "DeltaCt",
                "DeltaDeltaCt",
                "RelativeExpression",
            ]
            if gene.lower() == reference_gene.lower():
                target_label = f"{gene} Mean Cq (target)"
                reference_label = f"{reference_gene} Mean Cq (reference)"
            else:
                target_label = f"{gene} Mean Cq"
                reference_label = f"{reference_gene} Mean Cq"
            rename_map = {
                "TargetMeanCq": target_label,
                "ReferenceMeanCq": reference_label,
                "DeltaCt": "ΔCt",
                "DeltaDeltaCt": "ΔΔCt",
                "RelativeExpression": "2^-ΔΔCt",
            }
            st.dataframe(relative_gene[display_cols].rename(columns=rename_map), use_container_width=True, height=320)
            st.caption(
                "ΔCt = Target Mean Cq − Reference Mean Cq. ΔΔCt subtracts the calibrator pool mean ΔCt. "
                "Relative expression is 2^-ΔΔCt."
            )
        else:
            st.info("Relative expression could not be computed for this gene. Check calibrator selections and replicate data.")
        if relative_warning:
            st.warning(relative_warning)

    with tabs[5]:
        st.subheader("Step 6 – Missing standards audit")
        if not missing_standards_df.empty:
            st.dataframe(missing_standards_df, use_container_width=True, height=220)
            st.caption(
                "Genes listed here lack usable standard points. Without at least two standard levels with concentrations, "
                "quantities and curves cannot be fitted."
            )
        else:
            st.success("All genes have the minimum standard coverage required for curve fitting.")

    export_wells = filtered_wells.drop(columns=["IsNegativeControl", "Type_lower", "StateKey"], errors="ignore")
    if not export_wells.empty and "InputOrder" in export_wells.columns:
        export_wells = export_wells.sort_values("InputOrder").reset_index(drop=True)
    export_buffer = build_excel_stream(
        cleaned_wells=export_wells,
        standards_map_df=standards_map_df,
        replicate_averages=replicate_summary,
        std_curve_points=std_curve_points,
        curves=curves,
        sample_quantities=sample_quantities,
        normalised_df=normalised_df,
        relative_df=relative_all,
        missing_standards=missing_standards_df,
    )
    st.download_button(
        "Download qpcr_results.xlsx",
        data=export_buffer.getvalue(),
        file_name="qpcr_results.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        help="Export all tables into a single workbook for record keeping.",
    )

    st.caption(
        "Export includes: Cleaned wells, per-gene wells, replicate averages, standards map, standard curve points and fits, "
        "sample quantities, normalised results, relative expression, and the missing standards audit."
    )


if __name__ == "__main__":
    main()
