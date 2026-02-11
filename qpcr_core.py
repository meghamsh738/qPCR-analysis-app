from __future__ import annotations

import io
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

VALID_COLS = ["Plate", "Well", "Gene", "Type", "Label", "Replicate", "Cq"]
DEFAULT_TOP_CONC = 256.0
DEFAULT_DILUTION_FACTOR = 4.0
PLOT_BG = "#0b1224"
STANDARD_RE = re.compile(r"(?:^|\b)(std\d*|standard|curve|calib|calibrator)(?:\b|$)", re.IGNORECASE)
SAMPLE_RE = re.compile(r"(?:^|\b)(sample|unknown|unk)(?:\b|$)", re.IGNORECASE)
NEGATIVE_RE = re.compile(r"(?:^|\b)(ntc|no\s*template|negative|neg)(?:\b|$)", re.IGNORECASE)
POSITIVE_RE = re.compile(r"(?:^|\b)(positive|pos)(?:\b|$)", re.IGNORECASE)
BLANK_RE = re.compile(r"(?:^|\b)(blank|water)(?:\b|$)", re.IGNORECASE)


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


def _normalize_type_value(val: str) -> str:
    if val is None:
        return ""
    s = str(val).strip()
    if not s or s.lower() in {"nan", "none", "na"}:
        return ""
    if STANDARD_RE.search(s):
        return "Standard"
    if SAMPLE_RE.search(s):
        return "Sample"
    if NEGATIVE_RE.search(s):
        return "Negative"
    if POSITIVE_RE.search(s):
        return "Positive"
    if BLANK_RE.search(s):
        return "Blank"
    return s


def _infer_type_from_label(label: str) -> str:
    if label is None:
        return ""
    s = str(label).strip()
    if not s or s.lower() in {"nan", "none", "na"}:
        return ""
    if STANDARD_RE.search(s):
        return "Standard"
    if NEGATIVE_RE.search(s) or BLANK_RE.search(s):
        return "Blank"
    return ""


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

    df["Type"] = df["Type"].apply(_normalize_type_value)
    inferred = df["Label"].apply(_infer_type_from_label)
    df.loc[df["Type"] == "", "Type"] = inferred
    df.loc[df["Type"] == "", "Type"] = "Sample"

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
    fig.patch.set_facecolor("#0a1020")
    ax.set_facecolor("#111a2d")
    ax.grid(True, color=(1, 1, 1, 0.14), linewidth=0.9)

    ax.scatter(x, y, s=56, color="#38bdf8", edgecolor="#070b14", linewidth=0.6, label="Std mean Cq", zorder=3)
    ax.plot(xp, yp, linestyle="--", linewidth=2.2, color="#a78bfa", label=f"Fit (R²={r2:.3f}, Eff={eff:.1f}%)")

    ax.set_xlabel("log10(concentration)", color="#f8fafc")
    ax.set_ylabel("Cq", color="#f8fafc")
    ax.set_title(f"Standard curve — {gene}", color="#ffffff", pad=10, fontweight="bold")

    for spine in ax.spines.values():
        spine.set_color("#94a3b8")
    ax.tick_params(colors="#e2e8f0")
    legend = ax.legend(frameon=True, fontsize=9, facecolor="#0b1224", edgecolor="#475569", framealpha=0.95)
    if legend is not None:
        for text in legend.get_texts():
            text.set_color("#f8fafc")
    fig.tight_layout()
    return fig


def fit_standard_curve(std_df: pd.DataFrame, mapping: pd.DataFrame):
    m = mapping.dropna(subset=["Label", "Concentration"]).copy()
    m["Label"] = m["Label"].astype(str)
    m["Concentration"] = m["Concentration"].astype(float)
    m = m[(~m["Concentration"].isna()) & (m["Concentration"] > 0)].copy()
    m["log10_conc"] = np.log10(m["Concentration"])
    df = std_df.merge(m, on="Label", how="inner")
    g = df.groupby(["Gene", "Label", "log10_conc"], as_index=False)["Cq"].mean().rename(columns={"Cq": "meanCq"})
    curves = []
    for gene, sub in g.groupby("Gene"):
        sub = sub.dropna(subset=["log10_conc", "meanCq"]).copy()
        if sub.shape[0] < 2:
            curves.append({"Gene": gene, "slope": np.nan, "intercept": np.nan, "R2": np.nan, "Efficiency_%": np.nan, "n_points": sub.shape[0]})
            continue
        x = sub["log10_conc"].values
        y = sub["meanCq"].values
        slope, intercept = np.polyfit(x, y, 1)
        pred = slope * x + intercept
        R2 = r2_score(y, pred)
        eff = (10 ** (-1.0 / slope) - 1.0) * 100.0 if np.isfinite(slope) and slope != 0 else np.nan
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


def _std_curve_fit_metrics(points: pd.DataFrame):
    """Fit a standard curve from pre-aggregated points (log10_conc, meanCq)."""
    if points.empty:
        return {"slope": np.nan, "intercept": np.nan, "R2": np.nan, "Efficiency_%": np.nan, "n_points": 0}
    sub = points.dropna(subset=["log10_conc", "meanCq"]).copy()
    if sub.shape[0] < 2:
        return {"slope": np.nan, "intercept": np.nan, "R2": np.nan, "Efficiency_%": np.nan, "n_points": int(sub.shape[0])}
    x = sub["log10_conc"].to_numpy(dtype=float)
    y = sub["meanCq"].to_numpy(dtype=float)
    if np.allclose(x, x[0]):
        return {"slope": np.nan, "intercept": np.nan, "R2": np.nan, "Efficiency_%": np.nan, "n_points": int(sub.shape[0])}
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    R2 = r2_score(y, pred)
    eff = (10 ** (-1.0 / slope) - 1.0) * 100.0 if np.isfinite(slope) and slope != 0 else np.nan
    return {"slope": float(slope), "intercept": float(intercept), "R2": float(R2), "Efficiency_%": float(eff), "n_points": int(sub.shape[0])}


def suggest_standard_curve_exclusions(
    std_df: pd.DataFrame,
    mapping: pd.DataFrame,
    *,
    outlier_threshold: float = 0.75,
    min_levels_to_keep: int = 6,
    max_levels_to_drop: int = 1,
    r2_target: float = 0.98,
    eff_min: float = 90.0,
    eff_max: float = 110.0,
    max_combinations: int = 200_000,
    optimize_for: str = "balanced",
    replicate_drop_mode: str = "outlier_only",
):
    """Suggest which *standard wells* to exclude to get a more reliable standard curve.

    Hybrid approach:
    - Option 1 (replicate QC): only consider dropping a replicate if it's "outlier-like"
      within its standard level (Label) using `outlier_threshold`.
      * If 2 replicates: allow dropping either only if |ΔCq| > threshold.
      * If 3+ replicates: allow dropping any replicate whose |Cq - median(Cq)| > threshold.
    - Option 4 (global optimization): choose the best combination across all levels,
      with constraints (min kept levels, max dropped levels) and a scoring function
      that prefers meeting R²/efficiency targets with minimal exclusions.
    - Optional aggressive mode:
      * optimize_for="best_r2" prioritizes maximum R² (then fewer drops as tie-breakers).
      * replicate_drop_mode="any_single_replicate" allows considering dropping any one
        replicate in each level (not just outlier-like ones).

    Returns:
      exclusions_df: One row per suggested excluded well (Plate/Well/Gene/Label).
      summary_df: Per-curve (Gene) baseline vs suggested metrics + drop counts.
    """
    if std_df is None or std_df.empty:
        empty_excl = pd.DataFrame(columns=["Gene", "Label", "Concentration", "log10_conc", "Plate", "Well", "Cq", "reason"])
        empty_sum = pd.DataFrame(
            columns=[
                "Gene",
                "baseline_R2",
                "baseline_Efficiency_%",
                "baseline_n_points",
                "suggested_R2",
                "suggested_Efficiency_%",
                "suggested_n_points",
                "dropped_levels",
                "dropped_replicates",
                "excluded_wells",
            ]
        )
        return empty_excl, empty_sum

    m = mapping.dropna(subset=["Label", "Concentration"]).copy()
    m["Label"] = m["Label"].astype(str)
    m["Concentration"] = m["Concentration"].astype(float)
    m = m[(~m["Concentration"].isna()) & (m["Concentration"] > 0)].copy()
    if m.empty:
        empty_excl = pd.DataFrame(columns=["Gene", "Label", "Concentration", "log10_conc", "Plate", "Well", "Cq", "reason"])
        empty_sum = pd.DataFrame(
            columns=[
                "Gene",
                "baseline_R2",
                "baseline_Efficiency_%",
                "baseline_n_points",
                "suggested_R2",
                "suggested_Efficiency_%",
                "suggested_n_points",
                "dropped_levels",
                "dropped_replicates",
                "excluded_wells",
            ]
        )
        return empty_excl, empty_sum

    m["log10_conc"] = np.log10(m["Concentration"])

    optimize_for = str(optimize_for or "balanced").strip().lower()
    if optimize_for not in {"balanced", "best_r2"}:
        optimize_for = "balanced"

    replicate_drop_mode = str(replicate_drop_mode or "outlier_only").strip().lower()
    if replicate_drop_mode not in {"outlier_only", "any_single_replicate"}:
        replicate_drop_mode = "outlier_only"

    df = std_df.copy()
    for col in ["Plate", "Well"]:
        if col not in df.columns:
            df[col] = ""
    # Only rows that can participate in the curve.
    df = df.merge(m[["Label", "Concentration", "log10_conc"]], on="Label", how="inner")
    df = df.dropna(subset=["Gene", "Label", "Cq", "log10_conc"]).copy()
    if df.empty:
        empty_excl = pd.DataFrame(columns=["Gene", "Label", "Concentration", "log10_conc", "Plate", "Well", "Cq", "reason"])
        empty_sum = pd.DataFrame(
            columns=[
                "Gene",
                "baseline_R2",
                "baseline_Efficiency_%",
                "baseline_n_points",
                "suggested_R2",
                "suggested_Efficiency_%",
                "suggested_n_points",
                "dropped_levels",
                "dropped_replicates",
                "excluded_wells",
            ]
        )
        return empty_excl, empty_sum

    baseline_curves, _ = fit_standard_curve(std_df, m)
    baseline_by_gene = baseline_curves.set_index("Gene").to_dict(orient="index")

    exclusions_rows = []
    summary_rows = []

    def build_level_options(level_df: pd.DataFrame):
        level_df = level_df.dropna(subset=["Cq"]).copy()
        level_df = level_df.sort_values(["Plate", "Well"]).reset_index(drop=True)
        if level_df.empty:
            return [], []
        cqs = level_df["Cq"].to_numpy(dtype=float)
        n = int(level_df.shape[0])

        candidate_drop_idx = []
        if replicate_drop_mode == "any_single_replicate":
            if n >= 2:
                candidate_drop_idx = list(range(n))
        else:
            if n == 2:
                if abs(float(cqs[0]) - float(cqs[1])) > float(outlier_threshold):
                    candidate_drop_idx = [0, 1]
            elif n >= 3:
                med = float(np.median(cqs))
                deltas = np.abs(cqs - med)
                candidate_drop_idx = [int(i) for i in np.where(deltas > float(outlier_threshold))[0].tolist()]

        # Option: keep all replicates in the level.
        options = [{
            "drop_level": 0,
            "drop_reps": 0,
            "excluded_idx": [],
            "meanCq": float(np.mean(cqs)),
        }]

        # Option: drop one replicate, but only if it looks outlier-ish.
        for idx in candidate_drop_idx:
            keep_mask = np.ones(n, dtype=bool)
            keep_mask[idx] = False
            kept = cqs[keep_mask]
            if kept.size == 0:
                continue
            options.append({
                "drop_level": 0,
                "drop_reps": 1,
                "excluded_idx": [idx],
                "meanCq": float(np.mean(kept)),
            })

        # Option: drop the whole standard level.
        options.append({
            "drop_level": 1,
            "drop_reps": 0,
            "excluded_idx": list(range(n)),
            "meanCq": None,
        })

        return options, candidate_drop_idx

    def score_candidate(metrics: dict, drop_levels: int, drop_reps: int):
        slope = metrics.get("slope", np.nan)
        r2 = metrics.get("R2", np.nan)
        eff = metrics.get("Efficiency_%", np.nan)
        if not np.isfinite(slope) or slope >= 0:
            return None  # invalid for qPCR (Cq should decrease with log10 conc)

        if optimize_for == "best_r2":
            if not np.isfinite(r2):
                return None
            # Aggressive mode: maximize R² first, then keep fewer exclusions.
            # Lower tuple is better because the optimizer minimizes score.
            r2_term = -float(r2)
            eff_term = abs(float(eff) - 100.0) if np.isfinite(eff) else float("inf")
            return (r2_term, int(drop_levels), int(drop_reps), eff_term)

        violations = 0
        if (not np.isfinite(r2)) or (float(r2) < float(r2_target)):
            violations += 1
        if (not np.isfinite(eff)) or (float(eff) < float(eff_min)) or (float(eff) > float(eff_max)):
            violations += 1

        # Prefer: meeting targets > fewer drops > higher R2 > efficiency closer to 100%.
        r2_term = -float(r2) if np.isfinite(r2) else 0.0
        eff_term = abs(float(eff) - 100.0) if np.isfinite(eff) else float("inf")
        return (violations, int(drop_levels), int(drop_reps), r2_term, eff_term)

    for gene, gdf in df.groupby("Gene"):
        # Prepare levels in a deterministic order (by log10 conc then label).
        levels = []
        for (label, log10_conc), sub in gdf.groupby(["Label", "log10_conc"]):
            sub = sub.dropna(subset=["Cq"]).copy()
            if sub.empty:
                continue
            options, _cands = build_level_options(sub)
            levels.append({
                "Label": str(label),
                "log10_conc": float(log10_conc),
                "Concentration": float(sub["Concentration"].iloc[0]),
                "rows": sub.reset_index(drop=True),
                "options": options,
            })
        levels.sort(key=lambda d: (d["log10_conc"], d["Label"]))

        total_levels = len(levels)
        effective_min_keep = max(2, min(int(min_levels_to_keep), total_levels))
        effective_max_drop = max(0, min(int(max_levels_to_drop), max(0, total_levels - effective_min_keep)))

        best = None
        best_exclusions = []
        best_drop_levels = 0
        best_drop_reps = 0
        best_metrics = {"slope": np.nan, "intercept": np.nan, "R2": np.nan, "Efficiency_%": np.nan, "n_points": 0}

        # Quick exit: if we don't have enough levels to even fit, skip suggestions.
        if total_levels < 2:
            base = baseline_by_gene.get(gene, {})
            summary_rows.append({
                "Gene": gene,
                "baseline_R2": base.get("R2", np.nan),
                "baseline_Efficiency_%": base.get("Efficiency_%", np.nan),
                "baseline_n_points": base.get("n_points", 0),
                "suggested_R2": base.get("R2", np.nan),
                "suggested_Efficiency_%": base.get("Efficiency_%", np.nan),
                "suggested_n_points": base.get("n_points", 0),
                "dropped_levels": 0,
                "dropped_replicates": 0,
                "excluded_wells": 0,
            })
            continue

        # Estimate search space; fall back to greedy if too large.
        combos = 1
        for lv in levels:
            combos *= max(1, len(lv["options"]))

        def evaluate_selection(selected_options):
            points = []
            drop_levels = 0
            drop_reps = 0
            excluded = []
            for lv, opt in zip(levels, selected_options):
                if opt["drop_level"]:
                    drop_levels += 1
                drop_reps += int(opt["drop_reps"])
                rows = lv["rows"]
                for idx in opt["excluded_idx"]:
                    if opt["drop_level"]:
                        reason = "drop_level_fit"
                    else:
                        reason = "replicate_fit_optimization" if replicate_drop_mode == "any_single_replicate" else "replicate_outlier"
                    excluded.append({
                        "Gene": gene,
                        "Label": lv["Label"],
                        "Concentration": lv["Concentration"],
                        "log10_conc": lv["log10_conc"],
                        "Plate": rows.loc[idx, "Plate"],
                        "Well": rows.loc[idx, "Well"],
                        "Cq": float(rows.loc[idx, "Cq"]),
                        "reason": reason,
                    })
                if opt["meanCq"] is not None:
                    points.append({"log10_conc": lv["log10_conc"], "meanCq": float(opt["meanCq"])})

            if len(points) < effective_min_keep or len(points) < 2:
                return None
            pts_df = pd.DataFrame(points)
            metrics = _std_curve_fit_metrics(pts_df)
            score = score_candidate(metrics, drop_levels=drop_levels, drop_reps=drop_reps)
            if score is None:
                return None
            return score, metrics, excluded, drop_levels, drop_reps

        if combos <= int(max_combinations):
            # Exhaustive search with pruning.
            cur = [None] * total_levels

            def rec(i: int, drop_levels: int, drop_reps: int):
                nonlocal best, best_exclusions, best_drop_levels, best_drop_reps, best_metrics
                if drop_levels > effective_max_drop:
                    return
                # Even if we keep all remaining, can we reach min kept levels?
                max_keep_possible = (total_levels - drop_levels)
                if max_keep_possible < effective_min_keep:
                    return
                if i >= total_levels:
                    res = evaluate_selection(cur)
                    if res is None:
                        return
                    score, metrics, excluded, dlevels, dreps = res
                    if best is None or score < best:
                        best = score
                        best_metrics = metrics
                        best_exclusions = excluded
                        best_drop_levels = dlevels
                        best_drop_reps = dreps
                    return

                lv = levels[i]
                for opt in lv["options"]:
                    cur[i] = opt
                    rec(i + 1, drop_levels + int(opt["drop_level"]), drop_reps + int(opt["drop_reps"]))

            rec(0, 0, 0)
        else:
            # Greedy fallback: apply the single best local improvement step-by-step.
            selected = [lv["options"][0] for lv in levels]  # keep-all baseline
            base_res = evaluate_selection(selected)
            if base_res is not None:
                best, best_metrics, best_exclusions, best_drop_levels, best_drop_reps = base_res

            improved = True
            while improved:
                improved = False
                current_res = evaluate_selection(selected)
                if current_res is None:
                    break
                current_score, _m, _e, current_drop_levels, current_drop_reps = current_res
                best_step = None
                best_step_score = current_score
                for i, lv in enumerate(levels):
                    for opt in lv["options"]:
                        if opt is selected[i]:
                            continue
                        tmp = list(selected)
                        tmp[i] = opt
                        tmp_res = evaluate_selection(tmp)
                        if tmp_res is None:
                            continue
                        tmp_score, _tmp_m, _tmp_e, _dl, _dr = tmp_res
                        if tmp_score < best_step_score:
                            best_step_score = tmp_score
                            best_step = (i, opt, tmp_res)
                if best_step is not None:
                    idx, opt, tmp_res = best_step
                    selected[idx] = opt
                    best, best_metrics, best_exclusions, best_drop_levels, best_drop_reps = tmp_res
                    improved = True

        base = baseline_by_gene.get(gene, {})
        excluded_wells = len(best_exclusions) if best_exclusions else 0

        # Count dropped levels / replicate drops based on what's excluded (more intuitive than action counts).
        baseline_counts = gdf.groupby(["Label", "log10_conc"]).size().reset_index(name="n")
        excluded_df = pd.DataFrame(best_exclusions) if best_exclusions else pd.DataFrame(columns=["Label", "Plate", "Well"])
        excl_counts = excluded_df.groupby(["Label", "log10_conc"]).size().reset_index(name="n_excl") if not excluded_df.empty else pd.DataFrame(columns=["Label", "log10_conc", "n_excl"])
        merged_counts = baseline_counts.merge(excl_counts, on=["Label", "log10_conc"], how="left")
        merged_counts["n_excl"] = pd.to_numeric(merged_counts["n_excl"], errors="coerce").fillna(0).astype(int)
        merged_counts["n_remain"] = merged_counts["n"] - merged_counts["n_excl"]
        dropped_levels = int((merged_counts["n_remain"] <= 0).sum())
        dropped_rep_wells = int(merged_counts.loc[merged_counts["n_remain"] > 0, "n_excl"].sum())

        summary_rows.append({
            "Gene": gene,
            "baseline_R2": base.get("R2", np.nan),
            "baseline_Efficiency_%": base.get("Efficiency_%", np.nan),
            "baseline_n_points": base.get("n_points", 0),
            "suggested_R2": best_metrics.get("R2", np.nan),
            "suggested_Efficiency_%": best_metrics.get("Efficiency_%", np.nan),
            "suggested_n_points": best_metrics.get("n_points", 0),
            "dropped_levels": dropped_levels,
            "dropped_replicates": dropped_rep_wells,
            "excluded_wells": excluded_wells,
        })

        exclusions_rows.extend(best_exclusions)

    exclusions_df = pd.DataFrame(exclusions_rows)
    if exclusions_df.empty:
        exclusions_df = pd.DataFrame(columns=["Gene", "Label", "Concentration", "log10_conc", "Plate", "Well", "Cq", "reason"])
    else:
        exclusions_df = exclusions_df.drop_duplicates(subset=["Gene", "Label", "Plate", "Well"]).sort_values(
            ["Gene", "log10_conc", "Label", "Plate", "Well"]
        )

    summary_df = pd.DataFrame(summary_rows)
    if summary_df.empty:
        summary_df = pd.DataFrame(
            columns=[
                "Gene",
                "baseline_R2",
                "baseline_Efficiency_%",
                "baseline_n_points",
                "suggested_R2",
                "suggested_Efficiency_%",
                "suggested_n_points",
                "dropped_levels",
                "dropped_replicates",
                "excluded_wells",
            ]
        )
    else:
        summary_df = summary_df.sort_values(["Gene"])

    return exclusions_df, summary_df


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
