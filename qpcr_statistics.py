"""Small-sample statistical comparisons on biological qPCR endpoints.

Functions consume a single endpoint with one row per biological sample. They
never transform values or infer a pairing from row order. Pairwise tests in a
requested comparison family receive Holm-adjusted p-values.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats


def _truthy(value: Any) -> bool:
    if value is None:
        return False
    try:
        if bool(pd.isna(value)):
            return False
    except (TypeError, ValueError):
        pass
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        token = value.strip().casefold()
        if token in {"true", "yes", "y", "1", "eligible"}:
            return True
        if token in {"false", "no", "n", "0", "ineligible", ""}:
            return False
    raise ValueError("Eligible values must be true/false")


def _finite_values(values: pd.Series, label: str) -> np.ndarray:
    try:
        array = pd.to_numeric(values, errors="raise").to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("{} must contain numeric values".format(label)) from exc
    if not np.isfinite(array).all():
        raise ValueError("{} contains missing or non-finite values".format(label))
    return array


def _holm(p_values: Sequence[float]) -> List[float]:
    """Holm step-down adjusted p-values, in the input order."""
    if not p_values:
        return []
    indexed = sorted(enumerate(float(p) for p in p_values), key=lambda item: item[1])
    adjusted_sorted: List[float] = []
    running = 0.0
    count = len(indexed)
    for rank, (_, p_value) in enumerate(indexed):
        candidate = min(1.0, (count - rank) * p_value)
        running = max(running, candidate)
        adjusted_sorted.append(running)
    adjusted = [0.0] * count
    for sorted_index, (original_index, _) in enumerate(indexed):
        adjusted[original_index] = adjusted_sorted[sorted_index]
    return adjusted


def _has_usable_variation(values: np.ndarray) -> bool:
    scale = max(float(np.max(np.abs(values))) ** 2, 1e-300)
    variation = float(np.sum((values - np.mean(values)) ** 2))
    return variation > np.finfo(float).eps * scale * 100.0


def _ci(effect: float, standard_error: float, degrees_freedom: float) -> Tuple[float, float]:
    if standard_error == 0:
        return float(effect), float(effect)
    critical = float(stats.t.ppf(0.975, degrees_freedom))
    return float(effect - critical * standard_error), float(effect + critical * standard_error)


def _welch_test(a: np.ndarray, b: np.ndarray) -> Dict[str, Any]:
    n_a, n_b = len(a), len(b)
    mean_a, mean_b = float(np.mean(a)), float(np.mean(b))
    effect = mean_a - mean_b
    var_a = float(np.var(a, ddof=1))
    var_b = float(np.var(b, ddof=1))
    if not (_has_usable_variation(a) or _has_usable_variation(b)):
        raise ValueError("Welch test cannot estimate uncertainty when both groups have zero or near-zero variability")
    component_a, component_b = var_a / n_a, var_b / n_b
    se_squared = component_a + component_b
    if se_squared <= 0:
        raise ValueError("Welch test cannot estimate uncertainty from zero variability")
    else:
        statistic = effect / math.sqrt(se_squared)
        denominator = component_a ** 2 / (n_a - 1) + component_b ** 2 / (n_b - 1)
        degrees = se_squared ** 2 / denominator if denominator > 0 else math.inf
        p_value = float(2 * stats.t.sf(abs(statistic), degrees))
        lower, upper = _ci(effect, math.sqrt(se_squared), degrees)
    return {
        "n1": int(n_a), "n2": int(n_b), "N": int(n_a + n_b),
        "EffectMeanDifference": float(effect), "CI95Lower": lower, "CI95Upper": upper,
        "StatisticName": "t", "Statistic": float(statistic), "Df": float(degrees) if math.isfinite(degrees) else None,
        "PValue": float(p_value),
    }


def _paired_test(a: np.ndarray, b: np.ndarray) -> Dict[str, Any]:
    differences = a - b
    n = len(differences)
    effect = float(np.mean(differences))
    sd = float(np.std(differences, ddof=1))
    if not _has_usable_variation(differences):
        raise ValueError("Paired test cannot estimate uncertainty when paired differences have zero or near-zero variability")
    standard_error = sd / math.sqrt(n)
    statistic = effect / standard_error
    if not math.isfinite(statistic):
        raise ValueError("Paired test statistic is non-finite; review the endpoint scale and variability")
    p_value = float(2 * stats.t.sf(abs(statistic), n - 1))
    lower, upper = _ci(effect, standard_error, n - 1)
    return {
        "n1": int(n), "n2": int(n), "N": int(n),
        "EffectMeanDifference": effect, "CI95Lower": lower, "CI95Upper": upper,
        "StatisticName": "paired t", "Statistic": float(statistic), "Df": float(n - 1),
        "PValue": float(p_value),
    }


def _one_way_anova(frame: pd.DataFrame, groups: Sequence[str]) -> Dict[str, Any]:
    arrays = [frame.loc[frame["Group"] == group, "Value"].to_numpy(dtype=float) for group in groups]
    if any(len(values) < 2 for values in arrays):
        counts = {group: int(len(values)) for group, values in zip(groups, arrays)}
        raise ValueError("One-way ANOVA needs at least two eligible biological samples per group: {}".format(counts))
    total_n = sum(len(values) for values in arrays)
    grand_mean = float(np.mean(np.concatenate(arrays)))
    ss_between = sum(len(values) * (float(np.mean(values)) - grand_mean) ** 2 for values in arrays)
    ss_error = sum(float(np.sum((values - np.mean(values)) ** 2)) for values in arrays)
    df_between = len(groups) - 1
    df_error = total_n - len(groups)
    if df_error <= 0:
        raise ValueError("One-way ANOVA has no residual degrees of freedom")
    ms_between = ss_between / df_between
    mse = ss_error / df_error
    all_values = np.concatenate(arrays)
    if not _has_usable_variation(all_values) or ss_error <= np.finfo(float).eps * max(float(np.sum(all_values ** 2)), 1e-300) * 100.0:
        raise ValueError("One-way ANOVA cannot estimate residual uncertainty from zero or near-zero within-group variability")
    f_value = ms_between / mse
    if not math.isfinite(f_value):
        raise ValueError("One-way ANOVA F statistic is non-finite; review residual variability")
    p_value = float(stats.f.sf(f_value, df_between, df_error))
    return {
        "Method": "anova", "Test": "One-way ANOVA",
        "StatisticName": "F", "Statistic": float(f_value), "PValue": float(p_value),
        "DfBetween": int(df_between), "DfError": int(df_error), "N": int(total_n),
        "GroupNs": {group: int(len(values)) for group, values in zip(groups, arrays)},
        "MSE": float(mse), "SSBetween": float(ss_between), "SSError": float(ss_error),
    }


def _pairwise_anova(a: np.ndarray, b: np.ndarray, mse: float, df_error: int) -> Dict[str, Any]:
    mean_a, mean_b = float(np.mean(a)), float(np.mean(b))
    effect = mean_a - mean_b
    se_squared = mse * (1.0 / len(a) + 1.0 / len(b))
    if se_squared <= 0:
        raise ValueError("ANOVA contrast cannot estimate uncertainty from zero variability")
    standard_error = math.sqrt(se_squared)
    statistic = effect / standard_error
    if not math.isfinite(statistic):
        raise ValueError("ANOVA contrast statistic is non-finite; review residual variability")
    p_value = float(2 * stats.t.sf(abs(statistic), df_error))
    lower, upper = _ci(effect, standard_error, df_error)
    return {
        "n1": int(len(a)), "n2": int(len(b)), "N": int(len(a) + len(b)),
        "EffectMeanDifference": float(effect), "CI95Lower": lower, "CI95Upper": upper,
        "StatisticName": "pooled-error t contrast", "Statistic": float(statistic), "Df": float(df_error),
        "PValue": float(p_value),
    }


def run_statistics(
    results: pd.DataFrame,
    method: str,
    comparisons: Sequence[Tuple[str, str]],
    group_order: Optional[Sequence[str]] = None,
) -> Tuple[Dict[str, Any], pd.DataFrame]:
    """Run tests for one biological endpoint and a user-selected comparison family.

    ``welch`` and ``paired`` perform the requested two-group contrasts. For
    paired tests, ``PairId`` is the sole matching key. ``anova`` runs a
    one-way independent-groups omnibus test and selected contrasts using its
    pooled residual error. Rows marked ineligible are reported as exclusions;
    they do not silently contribute to sample sizes.
    """
    if not isinstance(results, pd.DataFrame):
        raise TypeError("results must be a pandas DataFrame")
    required = {"SampleId", "Group", "Value", "Eligible"}
    missing = sorted(required - set(results.columns))
    if missing:
        raise ValueError("results is missing required columns: {}".format(", ".join(missing)))
    if results.empty:
        raise ValueError("No biological-sample results were supplied")
    if "Gene" in results and results["Gene"].dropna().astype(str).nunique() > 1:
        raise ValueError("Select one gene before running statistics")
    if "Endpoint" in results and results["Endpoint"].dropna().astype(str).nunique() > 1:
        raise ValueError("Select one endpoint before running statistics")
    sample_ids = results["SampleId"].astype(str).str.strip()
    if sample_ids.eq("").any() or results["SampleId"].isna().any():
        raise ValueError("Every biological result needs a SampleId")
    if sample_ids.duplicated().any():
        duplicate_ids = sorted(sample_ids[sample_ids.duplicated(keep=False)].unique().tolist())
        raise ValueError("Results must have one row per biological SampleId; duplicates: {}".format(", ".join(duplicate_ids)))
    if "Group" in results:
        blank_eligible = results["Group"].map(lambda value: value is None or (isinstance(value, float) and math.isnan(value)) or not str(value).strip())
        eligible_flags = results["Eligible"].map(_truthy)
        if (blank_eligible & eligible_flags).any():
            raise ValueError("Eligible samples need a Group for comparison")
    eligible = results["Eligible"].map(_truthy)
    values = results.copy()
    values["SampleId"] = sample_ids
    values["Group"] = values["Group"].map(lambda value: str(value).strip() if not (value is None or (isinstance(value, float) and math.isnan(value))) else "")
    numeric = pd.to_numeric(values["Value"], errors="coerce")
    if (eligible & (~np.isfinite(numeric))).any():
        raise ValueError("Eligible biological results must have finite numeric values")
    excluded = values.loc[~eligible]
    values = values.loc[eligible].copy()
    values["Value"] = numeric.loc[eligible].astype(float)
    if values.empty:
        raise ValueError("No eligible biological samples are available for testing")
    if values["SampleId"].duplicated().any():
        raise ValueError("Eligible results contain duplicate biological SampleIds")
    if values["Group"].eq("").any():
        raise ValueError("Eligible samples need a Group for comparison")
    if "PairId" not in values:
        values["PairId"] = None

    selected_groups = list(group_order) if group_order is not None else list(dict.fromkeys(values["Group"].tolist()))
    selected_groups = [str(group).strip() for group in selected_groups if str(group).strip()]
    if len(set(selected_groups)) != len(selected_groups):
        raise ValueError("group_order contains duplicate group names")
    unknown_order = sorted(set(values["Group"].unique()) - set(selected_groups))
    if unknown_order:
        raise ValueError("group_order omits eligible groups: {}".format(", ".join(unknown_order)))
    if len(selected_groups) < 2:
        raise ValueError("At least two groups are required")
    if not comparisons:
        raise ValueError("Choose at least one planned comparison")
    normalized_comparisons: List[Tuple[str, str]] = []
    for pair in comparisons:
        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            raise ValueError("Each comparison must name exactly two groups")
        group_a, group_b = str(pair[0]).strip(), str(pair[1]).strip()
        if group_a == group_b:
            raise ValueError("A comparison must use two different groups")
        if group_a not in selected_groups or group_b not in selected_groups:
            raise ValueError("Comparison names a group outside group_order: {} vs {}".format(group_a, group_b))
        if (group_a, group_b) in normalized_comparisons or (group_b, group_a) in normalized_comparisons:
            raise ValueError("Comparison family contains a duplicate group pair")
        normalized_comparisons.append((group_a, group_b))

    code = str(method).strip().casefold()
    if code not in {"welch", "paired", "anova"}:
        raise ValueError("method must be 'welch', 'paired', or 'anova'")
    if code == "anova":
        omnibus = _one_way_anova(values, selected_groups)
        omnibus["ExcludedN"] = int(len(excluded))
        omnibus["ExcludedSampleIds"] = excluded["SampleId"].astype(str).tolist()
        mse, df_error = omnibus["MSE"], int(omnibus["DfError"])
    else:
        omnibus = {
            "Method": code,
            "Test": "Pairwise tests only; no omnibus test was requested",
            "StatisticName": None, "Statistic": None, "PValue": None,
            "DfBetween": None, "DfError": None,
            "N": int(len(values)),
            "GroupNs": {group: int((values["Group"] == group).sum()) for group in selected_groups},
            "ExcludedN": int(len(excluded)),
            "ExcludedSampleIds": excluded["SampleId"].astype(str).tolist(),
        }
        mse, df_error = np.nan, 0

    comparison_rows: List[Dict[str, Any]] = []
    for group_a, group_b in normalized_comparisons:
        a_frame = values.loc[values["Group"] == group_a]
        b_frame = values.loc[values["Group"] == group_b]
        if code == "paired":
            if a_frame["PairId"].map(lambda value: value is None or (isinstance(value, float) and math.isnan(value)) or not str(value).strip()).any() or b_frame["PairId"].map(lambda value: value is None or (isinstance(value, float) and math.isnan(value)) or not str(value).strip()).any():
                missing_ids = values.loc[
                    (values["Group"].isin([group_a, group_b]))
                    & values["PairId"].map(lambda value: value is None or (isinstance(value, float) and math.isnan(value)) or not str(value).strip()),
                    "SampleId",
                ].astype(str).tolist()
            else:
                missing_ids = []
            a_frame = a_frame.loc[~a_frame["PairId"].map(lambda value: value is None or (isinstance(value, float) and math.isnan(value)) or not str(value).strip())].copy()
            b_frame = b_frame.loc[~b_frame["PairId"].map(lambda value: value is None or (isinstance(value, float) and math.isnan(value)) or not str(value).strip())].copy()
            a_frame["PairId"] = a_frame["PairId"].astype(str).str.strip()
            b_frame["PairId"] = b_frame["PairId"].astype(str).str.strip()
            if a_frame["PairId"].duplicated().any() or b_frame["PairId"].duplicated().any():
                raise ValueError("PairId must identify exactly one biological sample per group")
            map_a = a_frame.set_index("PairId")
            map_b = b_frame.set_index("PairId")
            paired_ids = [pair_id for pair_id in map_a.index.tolist() if pair_id in map_b.index]
            unmatched_pair_ids = sorted((set(map_a.index) ^ set(map_b.index)))
            if len(paired_ids) < 2:
                raise ValueError("Paired comparison {} vs {} needs at least two explicit matched PairIds; found {}".format(group_a, group_b, len(paired_ids)))
            a = map_a.loc[paired_ids, "Value"].to_numpy(dtype=float)
            b = map_b.loc[paired_ids, "Value"].to_numpy(dtype=float)
            result = _paired_test(a, b)
            result["PairIds"] = list(paired_ids)
            result["UnmatchedPairIds"] = unmatched_pair_ids
            result["MissingPairIdSampleIds"] = missing_ids
            result["ExcludedSampleIds"] = sorted(
                a_frame.loc[~a_frame["PairId"].isin(paired_ids), "SampleId"].astype(str).tolist()
                + b_frame.loc[~b_frame["PairId"].isin(paired_ids), "SampleId"].astype(str).tolist()
                + missing_ids
            )
            result["Method"] = "paired"
        else:
            a = _finite_values(a_frame["Value"], "group {} values".format(group_a))
            b = _finite_values(b_frame["Value"], "group {} values".format(group_b))
            if len(a) < 2 or len(b) < 2:
                raise ValueError("Comparison {} vs {} needs at least two eligible biological samples per group".format(group_a, group_b))
            result = _welch_test(a, b) if code == "welch" else _pairwise_anova(a, b, mse, df_error)
            result["Method"] = code
            result["PairIds"] = []
            result["UnmatchedPairIds"] = []
            result["MissingPairIdSampleIds"] = []
            result["ExcludedSampleIds"] = []
        result["Group1"] = group_a
        result["Group2"] = group_b
        comparison_rows.append(result)

    adjusted = _holm([row["PValue"] for row in comparison_rows])
    for row, adjusted_p in zip(comparison_rows, adjusted):
        row["PAdjustedHolm"] = float(adjusted_p)
    if code == "welch" and len(comparison_rows) == 1:
        row = comparison_rows[0]
        omnibus.update({
            "Test": "Welch two-sample t-test",
            "StatisticName": row["StatisticName"],
            "Statistic": row["Statistic"],
            "PValue": row["PValue"],
            "DfBetween": 1,
            "DfError": row["Df"],
            "GroupNs": {row["Group1"]: row["n1"], row["Group2"]: row["n2"]},
        })
    if code == "paired" and len(comparison_rows) == 1:
        row = comparison_rows[0]
        omnibus.update({
            "Test": "Paired two-sample t-test",
            "StatisticName": row["StatisticName"],
            "Statistic": row["Statistic"],
            "PValue": row["PValue"],
            "DfBetween": 1,
            "DfError": row["Df"],
            "GroupNs": {row["Group1"]: row["N"], row["Group2"]: row["N"]},
        })
    return omnibus, pd.DataFrame(comparison_rows)
