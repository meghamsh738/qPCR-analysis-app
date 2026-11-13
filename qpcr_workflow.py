"""Scientific qPCR workflow calculations with explicit review states.

This module deliberately keeps raw well rows separate from biological-sample
summaries. A detected Cq is never silently removed, a standard curve is scoped
to one plate and gene, and any extrapolated quantity remains labelled.
"""

from __future__ import annotations

import json
import math
import re
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


WELL_COLUMNS = [
    "Plate", "Well", "Gene", "Type", "Label", "Cq", "OriginalCq",
    "RowId", "SampleId", "Group", "PairId", "Tissue", "Batch", "Keep", "Reason",
]
SUMMARY_COLUMNS = [
    "SampleId", "Label", "Group", "PairId", "Tissue", "Batch", "Gene", "Plate",
    "PlateTrace", "PlateCount", "MeanCq", "NumericCount",
    "RetainedWellCount", "TotalWellCount", "CqDelta", "Disagreement",
    "Eligible", "Status", "Reason",
]
CURVE_COLUMNS = [
    "Plate", "Gene", "Scope", "Slope", "Intercept", "R2",
    "EfficiencyPercent", "NLevels", "MinCq", "MaxCq",
    "MinConcentration", "MaxConcentration", "Eligible", "Status", "Issue",
]
STANDARD_POINT_COLUMNS = [
    "Plate", "Gene", "Label", "Concentration", "Log10Concentration",
    "Include", "MeanCq", "NumericCount", "WellCount", "Reason", "Issue",
]
RESULT_COLUMNS = [
    "SampleId", "Label", "Group", "PairId", "Gene", "Value",
    "Eligible", "Status", "Plate", "PlateTrace",
]
ISSUE_COLUMNS = ["Severity", "Issue", "Plate", "Gene", "SampleId", "RowId"]

_MISSING_TOKENS = {"", "na", "n/a", "nan", "none", "null", "missing"}
_UNDETERMINED_TOKENS = {
    "undetermined", "undet", "undetected", "no cq", "no ct", "no amplification",
}
_SAMPLE_ALIASES = {"sample", "unknown", "unk", "biological sample"}
_CANONICAL_TYPES = {"sample", "standard", "negative", "positive", "blank", "control"}
_CONTROL_TYPES = {
    "standard", "blank", "negative", "positive", "ntc", "no template control", "control",
}


def _is_blank(value: Any) -> bool:
    if value is None:
        return True
    try:
        if bool(pd.isna(value)):
            return True
    except (TypeError, ValueError):
        pass
    return isinstance(value, str) and not value.strip()


def _identifier(value: Any) -> Optional[str]:
    """Return a comparison key without changing the displayed source value."""
    if _is_blank(value):
        return None
    return str(value).strip()


def _append_message(existing: Any, message: str) -> str:
    prior = "" if _is_blank(existing) else str(existing).strip()
    if not prior:
        return message
    parts = [part.strip() for part in prior.split("; ") if part.strip()]
    if message not in parts:
        parts.append(message)
    return "; ".join(parts)


def _parse_cq(value: Any) -> Tuple[float, str, Optional[str]]:
    if _is_blank(value):
        return np.nan, "Missing", None
    if isinstance(value, str):
        text = value.strip()
        token = text.casefold()
        if token in _UNDETERMINED_TOKENS:
            return np.nan, "Undetermined", None
        if token in _MISSING_TOKENS:
            return np.nan, "Missing", None
    try:
        parsed = float(value)
    except (TypeError, ValueError, OverflowError):
        return np.nan, "Invalid", "Cq is not a numeric value"
    if not math.isfinite(parsed) or parsed < 0:
        return np.nan, "Invalid", "Cq must be a finite non-negative number"
    return parsed, "Detected", None


def _normalise_type(value: Any, label: Any) -> str:
    if _is_blank(value):
        label_key = (_identifier(label) or "").casefold()
        if re.search(r"(^|\b)(std\d*|standard|curve|calib(?:rator)?)(\b|$)", label_key):
            return "Standard"
        if re.search(r"(^|\b)(ntc|no\s*template|negative|neg)(\b|$)", label_key):
            return "Negative"
        if re.search(r"(^|\b)(blank|water)(\b|$)", label_key):
            return "Blank"
        return "Sample"
    text = str(value).strip()
    key = text.casefold()
    if key in _SAMPLE_ALIASES:
        return "Sample"
    if key in {"ntc", "no template", "no template control", "negative", "neg"}:
        return "Negative"
    if key in {"positive", "pos", "positive control"}:
        return "Positive"
    if key in {"standard", "std"} or re.match(r"^(std|standard)\b", key):
        return "Standard"
    if key in {"blank", "water"}:
        return "Blank"
    return text


def _parse_keep(value: Any) -> Tuple[bool, Optional[str]]:
    if _is_blank(value):
        return True, None
    if isinstance(value, (bool, np.bool_)):
        return bool(value), None
    if isinstance(value, (int, np.integer)) and value in (0, 1):
        return bool(value), None
    if isinstance(value, (float, np.floating)) and value in (0.0, 1.0):
        return bool(value), None
    if isinstance(value, str):
        token = value.strip().casefold()
        if token in {"true", "yes", "y", "keep", "1"}:
            return True, None
        if token in {"false", "no", "n", "exclude", "excluded", "0"}:
            return False, None
    return False, "Keep value is invalid; review the inclusion choice"


def _generated_sample_id(
    plate: Any, label: Any, row_id: Any, tissue: Any = None, batch: Any = None
) -> str:
    plate_key = _identifier(plate)
    label_key = _identifier(label)
    if plate_key is None or label_key is None:
        return "unmapped:" + str(row_id)
    # JSON tuple encoding avoids delimiter collisions while leaving Label intact.
    tissue_key, batch_key = _identifier(tissue), _identifier(batch)
    if tissue_key is None and batch_key is None:
        # Preserve the established key for existing files that have no context.
        key = [plate_key, label_key]
    else:
        key = [plate_key, label_key, tissue_key, batch_key]
    return "auto:" + json.dumps(key, ensure_ascii=False, separators=(",", ":"))


def _well_identity(row: pd.Series) -> Optional[Tuple[str, str, str]]:
    values = tuple(_identifier(row.get(col)) for col in ("Plate", "Well", "Gene"))
    if any(value is None for value in values):
        return None
    return values  # type: ignore[return-value]


def prepare_wells(raw: pd.DataFrame) -> pd.DataFrame:
    """Normalize mapped input while retaining every input row and source value.

    Required inputs are ``Plate, Well, Gene, Type, Label, Cq``. ``Type`` may
    be empty and is inferred conservatively from the label. Invalid Cq text is
    retained in ``OriginalCq`` and identified separately from missing or
    undetermined readings in ``CqStatus``.
    """
    if not isinstance(raw, pd.DataFrame):
        raise TypeError("raw must be a pandas DataFrame")
    result = raw.copy(deep=True).reset_index(drop=True)
    for col in ["Plate", "Well", "Gene", "Type", "Label", "Cq", "Tissue", "Batch"]:
        if col not in result.columns:
            result[col] = pd.Series([None] * len(result), dtype=object)
    if "OriginalCq" not in result.columns:
        result["OriginalCq"] = result["Cq"].copy()
    if "RowId" not in result.columns:
        result["RowId"] = ["row:{}".format(i + 1) for i in range(len(result))]
    if "SampleId" not in result.columns:
        result["SampleId"] = None
    if "Group" not in result.columns:
        result["Group"] = None
    if "PairId" not in result.columns:
        result["PairId"] = None
    if "Keep" not in result.columns:
        result["Keep"] = True
    if "Reason" not in result.columns:
        result["Reason"] = ""
    result["Issue"] = ""
    result["IdentityStatus"] = "Unique"
    result["CqStatus"] = "Missing"
    result["TypeStatus"] = "Recognized"

    for index, row in result.iterrows():
        row_id = row["RowId"]
        if _is_blank(row_id):
            row_id = "row:{}".format(index + 1)
            result.at[index, "RowId"] = row_id
            result.at[index, "Issue"] = _append_message(result.at[index, "Issue"], "Row ID is missing")
        cq, cq_status, cq_issue = _parse_cq(row["OriginalCq"])
        result.at[index, "Cq"] = cq
        result.at[index, "CqStatus"] = cq_status
        if cq_issue:
            result.at[index, "Issue"] = _append_message(result.at[index, "Issue"], cq_issue)
        normalized_type = _normalise_type(row["Type"], row["Label"])
        result.at[index, "Type"] = normalized_type
        if normalized_type.casefold() not in _CANONICAL_TYPES:
            result.at[index, "TypeStatus"] = "Unmapped type"
            result.at[index, "Issue"] = _append_message(
                result.at[index, "Issue"], "Unrecognized Type '{}'; map it explicitly".format(normalized_type)
            )
        keep, keep_issue = _parse_keep(row["Keep"])
        result.at[index, "Keep"] = keep
        if keep_issue:
            result.at[index, "Reason"] = _append_message(result.at[index, "Reason"], keep_issue)
            result.at[index, "Issue"] = _append_message(result.at[index, "Issue"], keep_issue)
        if not keep and _is_blank(row["Reason"]):
            result.at[index, "Issue"] = _append_message(result.at[index, "Issue"], "Exclusion reason is missing")
        sample_id = row["SampleId"]
        if _is_blank(sample_id):
            sample_id = _generated_sample_id(
                row["Plate"], row["Label"], row_id,
                tissue=row["Tissue"], batch=row["Batch"],
            )
        else:
            sample_id = str(sample_id).strip()
        result.at[index, "SampleId"] = sample_id

        missing_identity = [col for col in ("Plate", "Well", "Gene") if _identifier(row[col]) is None]
        if missing_identity:
            result.at[index, "IdentityStatus"] = "Missing identity"
            result.at[index, "Issue"] = _append_message(
                result.at[index, "Issue"], "Missing " + ", ".join(missing_identity)
            )
        if _identifier(row["Label"]) is None:
            result.at[index, "Issue"] = _append_message(result.at[index, "Issue"], "Sample/standard label is missing")

    identity_rows: Dict[Tuple[str, str, str], List[int]] = {}
    row_ids: Dict[str, List[int]] = {}
    for index, row in result.iterrows():
        identity = _well_identity(row)
        if identity:
            identity_rows.setdefault(identity, []).append(index)
        row_id = _identifier(row["RowId"])
        if row_id is not None:
            row_ids.setdefault(row_id, []).append(index)
    for indices in identity_rows.values():
        if len(indices) > 1:
            for index in indices:
                result.at[index, "IdentityStatus"] = "Duplicate well identity"
                result.at[index, "Issue"] = _append_message(
                    result.at[index, "Issue"], "Duplicate Plate/Well/Gene identity"
                )
    for indices in row_ids.values():
        if len(indices) > 1:
            for index in indices:
                result.at[index, "Issue"] = _append_message(result.at[index, "Issue"], "Duplicate RowId")

    # Existing identifiers/metadata remain in their source columns. Use
    # normalized keys only for conflict detection.
    by_sample: Dict[str, List[int]] = {}
    for index, sample_id in enumerate(result["SampleId"].tolist()):
        key = _identifier(sample_id)
        if key is not None:
            by_sample.setdefault(key, []).append(index)
    for indices in by_sample.values():
        for field in ("Group", "PairId", "Tissue", "Batch"):
            values = {_identifier(result.at[i, field]) for i in indices}
            values.discard(None)
            if len(values) > 1:
                message = "Conflicting {} metadata for SampleId".format(field)
                for index in indices:
                    result.at[index, "Issue"] = _append_message(result.at[index, "Issue"], message)
    return result


def _sample_groups(wells: pd.DataFrame) -> Iterable[Tuple[Tuple[str, str], List[int]]]:
    grouped: Dict[Tuple[str, str], List[int]] = {}
    for index, row in wells.iterrows():
        type_key = (_identifier(row.get("Type")) or "").casefold()
        type_unmapped = str(row.get("TypeStatus", "Recognized")) == "Unmapped type"
        if type_key != "sample" and not type_unmapped:
            continue
        sample_id = _identifier(row.get("SampleId"))
        gene = _identifier(row.get("Gene"))
        if sample_id is None:
            sample_id = "unmapped:" + str(_identifier(row.get("RowId")) or index)
        if gene is None:
            gene = "unmapped-gene:" + str(_identifier(row.get("RowId")) or index)
        grouped.setdefault((sample_id, gene), []).append(index)
    return grouped.items()


def summarise_samples(wells: pd.DataFrame, threshold: float = 0.75) -> pd.DataFrame:
    """Summarise technical Cq replicates to one row per biological sample/gene.

    ``CqDelta`` is the range (maximum minus minimum) among retained numeric
    technical readings. A large delta remains a review flag; once the well
    review step is confirmed, it does not by itself remove readings or block
    the biological sample.
    """
    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("threshold must be a finite non-negative Cq value")
    required = {"SampleId", "Gene", "Cq", "Keep", "Type", "Plate", "Label", "Group", "PairId"}
    missing = sorted(required - set(wells.columns))
    if missing:
        raise ValueError("wells is missing required columns: {}".format(", ".join(missing)))
    rows: List[Dict[str, Any]] = []
    for (sample_id, gene), indices in _sample_groups(wells):
        part = wells.loc[indices]
        plate_trace_list: List[str] = []
        for value in part["Plate"].tolist():
            key = _identifier(value)
            if key is not None and key not in plate_trace_list:
                plate_trace_list.append(key)
        plate_trace = " | ".join(plate_trace_list)
        plates_count = len(plate_trace_list)
        plate = part.iloc[0]["Plate"] if plates_count == 1 else None

        group_values = {_identifier(value) for value in part["Group"].tolist()}
        group_values.discard(None)
        pair_values = {_identifier(value) for value in part["PairId"].tolist()}
        pair_values.discard(None)
        context_values: Dict[str, set] = {}
        for field in ("Tissue", "Batch"):
            values = {
                _identifier(value)
                for value in part.get(field, pd.Series(None, index=part.index)).tolist()
            }
            values.discard(None)
            context_values[field] = values
        metadata_conflicts: List[str] = []
        if len(group_values) > 1:
            metadata_conflicts.append("Conflicting Group metadata")
        if len(pair_values) > 1:
            metadata_conflicts.append("Conflicting PairId metadata")
        for field, values in context_values.items():
            if len(values) > 1:
                metadata_conflicts.append("Conflicting {} metadata".format(field))
        source_issues = [
            str(value) for value in part.get("Issue", pd.Series("", index=part.index)).tolist()
            if not _is_blank(value)
        ]
        for field in ("Group", "PairId", "Tissue", "Batch"):
            message = "Conflicting {} metadata for SampleId".format(field)
            if any(message in issue for issue in source_issues):
                metadata_conflicts.append(message)
        metadata_conflicts = list(dict.fromkeys(metadata_conflicts))
        group = next(iter(group_values)) if len(group_values) == 1 else None
        pair_id = next(iter(pair_values)) if len(pair_values) == 1 else None
        tissue = next(iter(context_values["Tissue"])) if len(context_values["Tissue"]) == 1 else None
        batch = next(iter(context_values["Batch"])) if len(context_values["Batch"]) == 1 else None

        keep_mask = part["Keep"].map(lambda value: bool(value) if not _is_blank(value) else True)
        numeric_mask = part["CqStatus"].eq("Detected") if "CqStatus" in part.columns else part["Cq"].notna()
        valid_values = part.loc[keep_mask & numeric_mask, "Cq"].dropna().astype(float).to_numpy()
        numeric_count = int(len(valid_values))
        delta = float(np.max(valid_values) - np.min(valid_values)) if numeric_count >= 2 else np.nan
        disagreement = bool(numeric_count >= 2 and delta > threshold)
        mean_cq = float(np.mean(valid_values)) if numeric_count else np.nan

        issues = [str(v).strip() for v in part.get("Issue", pd.Series("", index=part.index)).tolist() if not _is_blank(v)]
        duplicate_identity = any(
            str(value) == "Duplicate well identity"
            for value in part.get("IdentityStatus", pd.Series("Unique", index=part.index)).tolist()
        )
        identity_missing = any(
            str(value) == "Missing identity"
            for value in part.get("IdentityStatus", pd.Series("Unique", index=part.index)).tolist()
        )
        invalid_cq = any(
            (part.get("CqStatus", pd.Series("Detected", index=part.index)) == "Invalid")
            & keep_mask
        )
        reasons: List[str] = []
        reviewer_reasons = [
            str(value).strip() for value in part.get("Reason", pd.Series("", index=part.index)).tolist()
            if not _is_blank(value)
        ]
        reasons.extend(reviewer_reasons)
        if plates_count > 1:
            reasons.append("Sample/Gene has numeric candidates on multiple plates: {}".format(plate_trace))
        if duplicate_identity:
            reasons.append("Duplicate Plate/Well/Gene identity requires review")
        if identity_missing:
            reasons.append("Missing Plate/Well/Gene identity requires review")
        reasons.extend(metadata_conflicts)
        if invalid_cq:
            reasons.append("Invalid Cq text requires review")
        if any(str(value) == "Unmapped type" for value in part.get("TypeStatus", pd.Series("Recognized", index=part.index)).tolist()):
            reasons.append("Unrecognized Type requires explicit mapping")
        if disagreement:
            reasons.append("Technical replicate Cq range exceeds {:.2f}; review flag".format(threshold))
        exclusion_reason_missing = any(
            (not bool(keep_mask.iloc[position])) and _is_blank(part.iloc[position].get("Reason"))
            for position in range(len(part))
        )
        if exclusion_reason_missing:
            reasons.append("Excluded well is missing a review reason")
        if numeric_count == 0:
            statuses = set(part.get("CqStatus", pd.Series("Missing", index=part.index)).tolist())
            if "Undetermined" in statuses and "Detected" not in statuses:
                reasons.append("Both/all replicates undetermined")
            elif bool((~keep_mask).all()):
                reasons.append("All replicates excluded by reviewer")
            else:
                reasons.append("No retained numeric Cq")
        elif numeric_count == 1:
            reasons.append("Single numeric replicate")

        unmapped_type = any(str(value) == "Unmapped type" for value in part.get("TypeStatus", pd.Series("Recognized", index=part.index)).tolist())
        unresolved = plates_count != 1 or duplicate_identity or identity_missing or bool(metadata_conflicts) or invalid_cq or exclusion_reason_missing or unmapped_type or numeric_count == 0
        eligible = not unresolved
        if plates_count > 1:
            status = "Ambiguous plate mapping"
        elif duplicate_identity:
            status = "Duplicate well identity requires review"
        elif identity_missing:
            status = "Missing well identity requires review"
        elif metadata_conflicts:
            status = "Conflicting sample metadata"
        elif invalid_cq:
            status = "Invalid Cq text requires review"
        elif unmapped_type:
            status = "Unrecognized Type requires mapping"
        elif numeric_count == 0:
            statuses = set(part.get("CqStatus", pd.Series("Missing", index=part.index)).tolist())
            if "Undetermined" in statuses and "Detected" not in statuses:
                status = "Both/all replicates undetermined"
            elif bool((~keep_mask).all()):
                status = "All replicates excluded by reviewer"
            else:
                status = "No numeric Cq"
        elif disagreement:
            status = "Replicate disagreement flagged"
        elif numeric_count == 1:
            status = "Single numeric replicate"
        else:
            status = "Ready"
        reason_text = "; ".join(dict.fromkeys(reasons))
        rows.append({
            "SampleId": part.iloc[0]["SampleId"],
            "Label": part.iloc[0]["Label"],
            "Group": group,
            "PairId": pair_id,
            "Tissue": tissue,
            "Batch": batch,
            "Gene": part.iloc[0]["Gene"],
            "Plate": plate,
            "PlateTrace": plate_trace,
            "PlateCount": plates_count,
            "MeanCq": mean_cq if plates_count == 1 and not duplicate_identity and not identity_missing else np.nan,
            "NumericCount": numeric_count,
            "RetainedWellCount": int(keep_mask.sum()),
            "TotalWellCount": int(len(part)),
            "CqDelta": delta,
            "Disagreement": disagreement,
            "Eligible": bool(eligible),
            "Status": status,
            "Reason": reason_text,
            "_Issue": "; ".join(dict.fromkeys(issues)),
        })
    output = pd.DataFrame(rows)
    if output.empty:
        output = pd.DataFrame(columns=SUMMARY_COLUMNS + ["_Issue"])
    else:
        output = output[SUMMARY_COLUMNS + ["_Issue"]]
    return output.reset_index(drop=True)


def _truthy(value: Any) -> bool:
    return _parse_keep(value)[0]


def make_standard_map(
    wells: pd.DataFrame,
    top: float = 16384.0,
    factor: float = 4.0,
) -> pd.DataFrame:
    """Map observed standard labels in plate order onto a serial dilution.

    The generated map is a starting point for reviewer confirmation; the
    concentrations and Include choices can be edited before curve fitting.
    """
    if not np.isfinite(top) or top <= 0:
        raise ValueError("top concentration must be finite and greater than zero")
    if not np.isfinite(factor) or factor <= 1:
        raise ValueError("dilution factor must be finite and greater than one")
    required = {"Plate", "Gene", "Label", "Type"}
    missing = sorted(required - set(wells.columns))
    if missing:
        raise ValueError("wells is missing required columns: {}".format(", ".join(missing)))
    selected = wells.loc[wells["Type"].map(lambda value: (_identifier(value) or "").casefold() == "standard")]
    rows: List[Dict[str, Any]] = []
    seen: Dict[Tuple[Optional[str], Optional[str]], List[str]] = {}
    for _, row in selected.iterrows():
        key = (_identifier(row["Plate"]), _identifier(row["Gene"]))
        label = _identifier(row["Label"])
        if key[0] is None or key[1] is None or label is None:
            continue
        labels = seen.setdefault(key, [])
        if label not in labels:
            labels.append(label)
    for (plate, gene), labels in seen.items():
        for order, label in enumerate(labels):
            concentration = float(top) / (float(factor) ** order)
            rows.append({
                "Plate": plate,
                "Gene": gene,
                "Label": label,
                "Concentration": concentration,
                "Include": True,
                # The generated dilution series is described in the UI. Keep
                # the reason blank so an excluded level requires a real reason.
                "Reason": "",
                "ObservedOrder": order + 1,
            })
    result = pd.DataFrame(rows)
    if result.empty:
        result = pd.DataFrame(columns=["Plate", "Gene", "Label", "Concentration", "Include", "Reason", "ObservedOrder"])
    return result.reset_index(drop=True)


def _fit_line(x: np.ndarray, y: np.ndarray) -> Tuple[float, float, float]:
    design = np.column_stack([x, np.ones(len(x), dtype=float)])
    slope, intercept = np.linalg.lstsq(design, y, rcond=None)[0]
    predicted = slope * x + intercept
    total = float(np.sum((y - np.mean(y)) ** 2))
    residual = float(np.sum((y - predicted) ** 2))
    r2 = float(1.0 - residual / total) if total > 0 else np.nan
    return float(slope), float(intercept), r2


def fit_curves(
    wells: pd.DataFrame,
    standards: pd.DataFrame,
    scope: str = "plate",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Fit Cq against log10 concentration per plate and gene.

    Only levels marked ``Include`` are used. Selected missing readings or
    duplicate concentrations invalidate the curve instead of being skipped
    without notice. A single numeric well at a selected level is allowed but
    called out in the standard-points table.
    """
    if str(scope).casefold() != "plate":
        raise ValueError("Only plate-scoped curves are supported; pooling plates must be explicit")
    required_wells = {"Plate", "Well", "Gene", "Type", "Label", "Cq", "Keep"}
    missing = sorted(required_wells - set(wells.columns))
    if missing:
        raise ValueError("wells is missing required columns: {}".format(", ".join(missing)))
    required_standards = {"Plate", "Gene", "Label", "Concentration", "Include"}
    missing = sorted(required_standards - set(standards.columns))
    if missing:
        raise ValueError("standards is missing required columns: {}".format(", ".join(missing)))

    std_wells = wells.loc[wells["Type"].map(lambda value: (_identifier(value) or "").casefold() == "standard")]
    map_groups: Dict[Tuple[str, str], List[pd.Series]] = {}
    for _, row in standards.iterrows():
        plate, gene = _identifier(row["Plate"]), _identifier(row["Gene"])
        if plate is not None and gene is not None:
            map_groups.setdefault((plate, gene), []).append(row)
    for _, row in std_wells.iterrows():
        plate, gene = _identifier(row["Plate"]), _identifier(row["Gene"])
        if plate is not None and gene is not None:
            map_groups.setdefault((plate, gene), [])

    point_rows: List[Dict[str, Any]] = []
    curve_rows: List[Dict[str, Any]] = []
    for (plate, gene), mappings in map_groups.items():
        mapping_by_label: Dict[str, List[pd.Series]] = {}
        for mapping in mappings:
            label = _identifier(mapping.get("Label"))
            if label is not None:
                mapping_by_label.setdefault(label, []).append(mapping)
        observed_labels: List[str] = []
        for _, row in std_wells.iterrows():
            if _identifier(row["Plate"]) == plate and _identifier(row["Gene"]) == gene:
                label = _identifier(row["Label"])
                if label is not None and label not in observed_labels:
                    observed_labels.append(label)
        labels = list(dict.fromkeys(list(mapping_by_label.keys()) + observed_labels))
        curve_issues: List[str] = []
        selected_points: List[Dict[str, Any]] = []
        selected_missing = False
        missing_exclusion_reason = False
        duplicate_mapping = False
        for label in labels:
            matching_maps = mapping_by_label.get(label, [])
            mapping_issue = ""
            concentration = np.nan
            include = False
            reason = ""
            if len(matching_maps) > 1:
                duplicate_mapping = True
                mapping_issue = "Duplicate Plate/Gene/Label standard mapping"
                curve_issues.append(mapping_issue)
            elif len(matching_maps) == 1:
                mapping = matching_maps[0]
                try:
                    concentration = float(mapping["Concentration"])
                except (TypeError, ValueError, OverflowError):
                    concentration = np.nan
                include = _truthy(mapping["Include"])
                reason = "" if _is_blank(mapping.get("Reason")) else str(mapping.get("Reason"))
                if not np.isfinite(concentration) or concentration <= 0:
                    mapping_issue = "Standard concentration must be finite and positive"
                    curve_issues.append(mapping_issue)
                    if include:
                        selected_missing = True
                if not include and not reason:
                    missing_exclusion_reason = True
                    mapping_issue = _append_message(
                        mapping_issue, "Excluded standard level is missing a review reason"
                    )
                    curve_issues.append(
                        "Excluded standard level {} is missing a review reason".format(label)
                    )
            else:
                mapping_issue = "Observed standard has no concentration mapping"
                curve_issues.append(mapping_issue)
                selected_missing = True

            level = std_wells.loc[
                std_wells["Plate"].map(_identifier).eq(plate)
                & std_wells["Gene"].map(_identifier).eq(gene)
                & std_wells["Label"].map(_identifier).eq(label)
            ]
            excluded_level_wells = ~level["Keep"].map(_truthy)
            level_reasons = level.get("Reason", pd.Series("", index=level.index))
            missing_well_exclusion_reason = bool(
                (excluded_level_wells & level_reasons.map(_is_blank)).any()
            )
            duplicate_well = False
            if "IdentityStatus" in level.columns:
                duplicate_well = level["IdentityStatus"].eq("Duplicate well identity").any()
            keep = level["Keep"].map(_truthy)
            detected = level["CqStatus"].eq("Detected") if "CqStatus" in level.columns else level["Cq"].notna()
            values = level.loc[keep & detected, "Cq"].dropna().astype(float).to_numpy()
            invalid_selected = False
            if "CqStatus" in level.columns:
                invalid_selected = bool((keep & level["CqStatus"].isin(["Invalid", "Undetermined", "Missing"])).any())
            point_issue = mapping_issue
            if missing_well_exclusion_reason:
                point_issue = _append_message(
                    point_issue, "Excluded standard well is missing a review reason"
                )
                curve_issues.append(
                    "Excluded standard well is missing a review reason at {}".format(label)
                )
                missing_exclusion_reason = True
            if duplicate_well:
                point_issue = _append_message(point_issue, "Duplicate well identity in standard level")
                curve_issues.append("Duplicate well identity in standard level {}".format(label))
                if include:
                    selected_missing = True
            if include and 'CqStatus' in level.columns and (keep & level['CqStatus'].eq('Invalid')).any():
                selected_missing = True
                curve_issues.append('Invalid Cq text in a selected standard; correct or explicitly exclude it')
            if include and invalid_selected:
                point_issue = _append_message(point_issue, "One or more selected standard wells have no numeric Cq")
                curve_issues.append("One or more selected standard wells have no numeric Cq at {}".format(label))
            if include and len(values) == 0:
                point_issue = _append_message(point_issue, "Selected standard level has no retained numeric Cq")
                if point_issue not in curve_issues:
                    curve_issues.append("Selected standard level has no retained numeric Cq at {}".format(label))
                selected_missing = True
            mean_cq = float(np.mean(values)) if len(values) else np.nan
            if include and len(values) == 1:
                point_issue = _append_message(point_issue, "Only one retained numeric standard reading")
            point = {
                "Plate": plate,
                "Gene": gene,
                "Label": label,
                "Concentration": concentration,
                "Log10Concentration": float(np.log10(concentration)) if np.isfinite(concentration) and concentration > 0 else np.nan,
                "Include": bool(include),
                "MeanCq": mean_cq,
                "NumericCount": int(len(values)),
                "WellCount": int(len(level)),
                "Reason": reason,
                "Issue": point_issue,
            }
            point_rows.append(point)
            if include and np.isfinite(mean_cq) and np.isfinite(concentration) and concentration > 0 and not duplicate_well:
                selected_points.append(point)

        included_concentrations = [point["Concentration"] for point in selected_points]
        duplicate_concentration = len(included_concentrations) != len(set(included_concentrations))
        if duplicate_concentration:
            curve_issues.append("Duplicate selected standard concentrations")
        x = np.array([point["Log10Concentration"] for point in selected_points], dtype=float)
        y = np.array([point["MeanCq"] for point in selected_points], dtype=float)
        distinct_levels = len(set(included_concentrations))
        slope = intercept = r2 = efficiency = np.nan
        min_cq = max_cq = min_conc = max_conc = np.nan
        eligible = False
        status = "Invalid standard curve"
        if duplicate_mapping:
            status = "Duplicate standard mapping"
        elif duplicate_concentration:
            status = "Duplicate standard concentration"
        elif selected_missing or missing_exclusion_reason:
            status = "Selected standard mapping or reading requires review"
        elif distinct_levels < 2:
            curve_issues.append("At least two distinct selected concentrations are required")
            status = "Insufficient distinct standard levels"
        elif not np.isfinite(x).all() or not np.isfinite(y).all():
            curve_issues.append("Selected standards contain missing numeric values")
            status = "Selected standard values are incomplete"
        else:
            slope, intercept, r2 = _fit_line(x, y)
            min_cq, max_cq = float(np.min(y)), float(np.max(y))
            min_conc, max_conc = float(min(included_concentrations)), float(max(included_concentrations))
            if not np.isfinite(slope) or slope >= 0:
                curve_issues.append("Cq versus log10 concentration must have a negative slope")
                status = "Invalid non-negative slope"
            elif not np.isfinite(r2):
                curve_issues.append("R-squared cannot be calculated from constant standard Cq values")
                status = "Invalid standard fit"
            else:
                with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                    efficiency = float((np.power(10.0, -1.0 / slope) - 1.0) * 100.0)
                eligible = bool(np.isfinite(efficiency))
                status = "Ready" if eligible else "Invalid efficiency calculation"
                if not eligible:
                    curve_issues.append("Efficiency is not finite")
        curve_rows.append({
            "Plate": plate,
            "Gene": gene,
            "Scope": "plate",
            "Slope": slope,
            "Intercept": intercept,
            "R2": r2,
            "EfficiencyPercent": efficiency,
            "NLevels": int(distinct_levels),
            "MinCq": min_cq,
            "MaxCq": max_cq,
            "MinConcentration": min_conc,
            "MaxConcentration": max_conc,
            "Eligible": eligible,
            "Status": status,
            "Issue": "; ".join(dict.fromkeys(curve_issues)),
        })
    curves = pd.DataFrame(curve_rows, columns=CURVE_COLUMNS)
    points = pd.DataFrame(point_rows, columns=STANDARD_POINT_COLUMNS)
    return curves.reset_index(drop=True), points.reset_index(drop=True)


def quantify(
    summary: pd.DataFrame,
    curves: pd.DataFrame,
    allow_extrapolation: bool = False,
) -> pd.DataFrame:
    """Back-calculate quantity per sample/gene from its plate-specific curve."""
    required_summary = {"SampleId", "Gene", "Plate", "MeanCq", "Eligible", "Status", "PlateTrace"}
    missing = sorted(required_summary - set(summary.columns))
    if missing:
        raise ValueError("summary is missing required columns: {}".format(", ".join(missing)))
    required_curves = {"Plate", "Gene", "Slope", "Intercept", "MinCq", "MaxCq", "MinConcentration", "MaxConcentration", "Eligible", "Status"}
    missing = sorted(required_curves - set(curves.columns))
    if missing:
        raise ValueError("curves is missing required columns: {}".format(", ".join(missing)))
    curve_map: Dict[Tuple[str, str], pd.Series] = {}
    for _, curve in curves.iterrows():
        key = (_identifier(curve["Plate"]), _identifier(curve["Gene"]))
        if None not in key:
            if key in curve_map:
                raise ValueError("Duplicate plate/gene curves are ambiguous")
            curve_map[key] = curve
    records: List[Dict[str, Any]] = []
    for _, row in summary.iterrows():
        record = row.to_dict()
        plate, gene = _identifier(row["Plate"]), _identifier(row["Gene"])
        curve = curve_map.get((plate, gene)) if plate is not None and gene is not None else None
        record["Quantity"] = np.nan
        record["RangeStatus"] = "Unavailable"
        record["ObservedRangeStatus"] = "Unavailable"
        record["ConcentrationRangeStatus"] = "Unavailable"
        record["Eligible"] = False
        if curve is None:
            record["Status"] = "No plate/gene standard curve"
            record["Reason"] = _append_message(record.get("Reason", ""), "No plate/gene standard curve")
            records.append(record)
            continue
        if not _truthy(curve["Eligible"]):
            record["Status"] = "Standard curve requires review: {}".format(curve["Status"])
            record["Reason"] = _append_message(record.get("Reason", ""), curve.get("Issue", curve["Status"]))
            records.append(record)
            continue
        mean_cq = row["MeanCq"]
        if not np.isfinite(mean_cq):
            record["Status"] = row["Status"]
            record["RangeStatus"] = "No numeric Cq"
            records.append(record)
            continue
        mean_cq = float(mean_cq)
        min_cq, max_cq = float(curve["MinCq"]), float(curve["MaxCq"])
        slope, intercept = float(curve["Slope"]), float(curve["Intercept"])
        if mean_cq < min_cq:
            observed_range_status = "AboveRange"
        elif mean_cq > max_cq:
            observed_range_status = "BelowRange"
        else:
            observed_range_status = "InRange"
        with np.errstate(over="ignore", under="ignore", invalid="ignore", divide="ignore"):
            quantity = float(np.power(10.0, (mean_cq - intercept) / slope))
        if not np.isfinite(quantity) or quantity <= 0:
            record["Status"] = "Calculated quantity is invalid"
            record["Reason"] = _append_message(record.get("Reason", ""), "Calculated quantity is not finite and positive")
            records.append(record)
            continue
        min_concentration = float(curve["MinConcentration"])
        max_concentration = float(curve["MaxConcentration"])
        tolerance = 1e-9 * max(abs(min_concentration), abs(max_concentration), abs(quantity), 1e-300)
        if quantity < min_concentration - tolerance:
            concentration_range_status = "BelowRange"
        elif quantity > max_concentration + tolerance:
            concentration_range_status = "AboveRange"
        else:
            concentration_range_status = "InRange"
        if observed_range_status != "InRange" and concentration_range_status != "InRange" and observed_range_status != concentration_range_status:
            range_status = "RangeMismatch"
        elif concentration_range_status != "InRange":
            range_status = concentration_range_status
        else:
            range_status = observed_range_status
        record["ObservedRangeStatus"] = observed_range_status
        record["ConcentrationRangeStatus"] = concentration_range_status
        record["RangeStatus"] = range_status
        if range_status == "RangeMismatch":
            record["Status"] = "Cq and fitted concentration ranges disagree; curve review required"
            record["Reason"] = _append_message(record.get("Reason", ""), "Observed Cq range conflicts with back-calculated concentration range")
            records.append(record)
            continue
        if range_status != "InRange" and not allow_extrapolation:
            record["Status"] = "Outside calibrated range; final quantity withheld"
            record["Reason"] = _append_message(record.get("Reason", ""), "Observed Cq or fitted quantity is outside selected standards")
            records.append(record)
            continue
        record["Quantity"] = quantity
        sample_eligible = _truthy(row["Eligible"])
        if range_status != "InRange":
            record["Status"] = "Provisional extrapolation"
            record["Reason"] = _append_message(record.get("Reason", ""), "Outside selected standard range; calculated by explicit choice")
            record["Eligible"] = sample_eligible
        else:
            record["Status"] = row["Status"] if not sample_eligible else "Ready"
            record["Eligible"] = sample_eligible
        records.append(record)
    output = pd.DataFrame(records)
    if output.empty:
        output = summary.copy()
        output["Quantity"] = pd.Series(dtype=float)
        output["RangeStatus"] = pd.Series(dtype=object)
    return output.reset_index(drop=True)


def _as_names(values: Optional[Sequence[Any]]) -> List[str]:
    if values is None:
        return []
    return [str(value).strip() for value in values if not _is_blank(value)]


def _add_issue(
    issues: List[Dict[str, Any]],
    issue: str,
    severity: str = "Review",
    plate: Any = None,
    gene: Any = None,
    sample_id: Any = None,
    row_id: Any = None,
) -> None:
    issues.append({
        "Severity": severity,
        "Issue": issue,
        "Plate": plate,
        "Gene": gene,
        "SampleId": sample_id,
        "RowId": row_id,
    })


def analyze(
    wells: pd.DataFrame,
    standards: Optional[pd.DataFrame],
    references: Sequence[str],
    calibrators: Sequence[str],
    allow_extrapolation: bool = False,
    efficiency_confirmed: bool = False,
    threshold: float = 0.75,
) -> Dict[str, pd.DataFrame]:
    """Run reviewed standard-curve and ΔΔCq endpoints from prepared wells.

    ``efficiency_confirmed`` is an explicit researcher acceptance that the
    2^-ΔΔCq assumptions are suitable for the selected target/reference assays;
    it is not an estimate or claim that measured efficiency is exactly 100%.
    """
    if not isinstance(wells, pd.DataFrame):
        raise TypeError("wells must be a pandas DataFrame")
    prepared = wells.copy(deep=True) if {"CqStatus", "IdentityStatus", "OriginalCq"}.issubset(wells.columns) else prepare_wells(wells)
    references_list = _as_names(references)
    calibrator_list = _as_names(calibrators)
    if len(set(references_list)) != len(references_list):
        raise ValueError("references must not contain duplicates")
    if len(set(calibrator_list)) != len(calibrator_list):
        raise ValueError("calibrators must not contain duplicate SampleIds")

    issues: List[Dict[str, Any]] = []
    if not references_list:
        _add_issue(issues, "No reference genes selected; normalized endpoints are unavailable")
    if not calibrator_list:
        _add_issue(issues, "No ΔΔCq calibrator SampleIds selected; relative expression is unavailable")

    for _, row in prepared.iterrows():
        if not _is_blank(row.get("Issue")):
            for issue in str(row["Issue"]).split("; "):
                _add_issue(issues, issue, plate=row.get("Plate"), gene=row.get("Gene"), sample_id=row.get("SampleId"), row_id=row.get("RowId"))

    summary = summarise_samples(prepared, threshold=threshold)
    for _, row in summary.iterrows():
        if not _is_blank(row.get("_Issue")):
            for issue in str(row["_Issue"]).split("; "):
                _add_issue(issues, issue, plate=row.get("PlateTrace"), gene=row.get("Gene"), sample_id=row.get("SampleId"))
        if not _truthy(row["Eligible"]):
            _add_issue(issues, str(row["Status"]), plate=row.get("PlateTrace"), gene=row.get("Gene"), sample_id=row.get("SampleId"))
    for sample_id, group in summary.groupby(summary["SampleId"].map(_identifier), sort=False, dropna=False):
        for field in ("Group", "PairId", "Tissue", "Batch"):
            values = {_identifier(value) for value in group[field].tolist()}
            values.discard(None)
            if len(values) > 1:
                _add_issue(issues, "Conflicting {} metadata across genes".format(field), sample_id=sample_id)

    if standards is None:
        standards = make_standard_map(prepared)
    curves, standard_points = fit_curves(prepared, standards, scope="plate")
    for _, curve in curves.iterrows():
        if not _is_blank(curve.get("Issue")):
            _add_issue(issues, str(curve["Issue"]), plate=curve.get("Plate"), gene=curve.get("Gene"))
        if not _truthy(curve["Eligible"]):
            _add_issue(issues, str(curve["Status"]), plate=curve.get("Plate"), gene=curve.get("Gene"))
    quantities = quantify(summary, curves, allow_extrapolation=allow_extrapolation)
    reference_set = set(references_list)
    quantity_rows: Dict[Tuple[str, str], pd.Series] = {}
    summary_rows: Dict[Tuple[str, str], pd.Series] = {}
    for _, row in quantities.iterrows():
        key = (_identifier(row["SampleId"]), _identifier(row["Gene"]))
        if None not in key:
            quantity_rows[key] = row
            summary_rows[key] = summary.loc[
                summary["SampleId"].map(_identifier).eq(key[0]) & summary["Gene"].map(_identifier).eq(key[1])
            ].iloc[0]

    standard_results: List[Dict[str, Any]] = []
    ddcq_results: List[Dict[str, Any]] = []
    target_genes = [
        str(gene) for gene in summary["Gene"].dropna().drop_duplicates().tolist()
        if str(gene) not in reference_set
    ]
    sample_ids = list(dict.fromkeys(str(value) for value in summary["SampleId"].dropna().tolist()))
    for sample_id in sample_ids:
        sample_targets = [row for (sid, gene), row in quantity_rows.items() if sid == sample_id and gene not in reference_set]
        for target in sample_targets:
            target_gene = _identifier(target["Gene"]) or ""
            ref_rows = [quantity_rows.get((sample_id, reference)) for reference in references_list]
            missing_refs = [reference for reference, ref_row in zip(references_list, ref_rows) if ref_row is None]
            invalid_refs = [
                reference for reference, ref_row in zip(references_list, ref_rows)
                if ref_row is not None and (not _truthy(ref_row.get("Eligible")) or not np.isfinite(ref_row.get("Quantity", np.nan)) or float(ref_row.get("Quantity", np.nan)) <= 0)
            ]
            target_value = target.get("Quantity", np.nan)
            target_ok = _truthy(target.get("Eligible")) and np.isfinite(target_value) and float(target_value) > 0
            value = np.nan
            eligible = False
            statuses: List[str] = []
            if target.get("Status") == "Provisional extrapolation":
                statuses.append("Target extrapolated")
            if not target_ok:
                statuses.append("Target quantity unavailable or ineligible: {}".format(target.get("Status", "unknown")))
            if missing_refs:
                statuses.append("Missing reference gene(s): {}".format(", ".join(missing_refs)))
                for reference in missing_refs:
                    _add_issue(issues, "Missing reference measurement", plate=target.get("Plate"), gene=reference, sample_id=sample_id)
            if invalid_refs:
                statuses.append("Reference quantity unavailable or ineligible: {}".format(", ".join(invalid_refs)))
                for reference in invalid_refs:
                    ref_row = quantity_rows.get((sample_id, reference))
                    _add_issue(issues, "Reference quantity unavailable or ineligible", plate=ref_row.get("Plate") if ref_row is not None else None, gene=reference, sample_id=sample_id)
            ref_quantities = [float(row["Quantity"]) for row in ref_rows if row is not None and np.isfinite(row.get("Quantity", np.nan)) and float(row.get("Quantity", np.nan)) > 0]
            if references_list and not missing_refs and not invalid_refs and len(ref_quantities) == len(references_list) and target_ok:
                reference_quantity = float(np.exp(np.mean(np.log(ref_quantities))))
                value = float(target_value) / reference_quantity
                eligible = bool(np.isfinite(value) and value > 0)
                ref_extrapolated = any(row.get("Status") == "Provisional extrapolation" for row in ref_rows if row is not None)
                if ref_extrapolated:
                    statuses.append("Reference quantity extrapolated")
                if statuses:
                    statuses.append("Provisional standard-curve ratio")
                else:
                    statuses = ["Ready"]
            else:
                reference_quantity = np.nan
                if not references_list:
                    statuses.append("No reference genes selected")
            if not eligible and not statuses:
                statuses = ["Normalized standard-curve result unavailable"]
            reference_plates = ["{}:{}".format(reference, _identifier(row.get("Plate")) or "unknown") for reference, row in zip(references_list, ref_rows) if row is not None]
            plate_trace = "Target {} | References {}".format(
                _identifier(target.get("PlateTrace")) or "unknown",
                ", ".join(reference_plates) or "missing",
            )
            standard_results.append({
                "Endpoint": "Target/reference ratio",
                "SampleId": target["SampleId"],
                "Label": target.get("Label"),
                "Group": target.get("Group"),
                "PairId": target.get("PairId"),
                "Gene": target_gene,
                "Value": value,
                "Eligible": bool(eligible),
                "Status": "; ".join(dict.fromkeys(statuses)),
                "Plate": target.get("Plate"),
                "PlateTrace": plate_trace,
                "TargetQuantity": target_value,
                "ReferenceQuantity": reference_quantity,
                "ReferenceGenes": ", ".join(references_list),
                "RangeStatus": target.get("RangeStatus"),
            })

    # Give the curve-based ratio a calibrator-relative expression endpoint as
    # a separate row. Each selected biological calibrator contributes one
    # ratio, regardless of how many technical wells it had.
    standard_ratio_rows = list(standard_results)
    for target_gene in target_genes:
        gene_rows = [row for row in standard_ratio_rows if _identifier(row.get("Gene")) == target_gene]
        by_id = {_identifier(row.get("SampleId")): row for row in gene_rows}
        missing_calibrators = [sample_id for sample_id in calibrator_list if by_id.get(sample_id) is None]
        invalid_calibrators = [
            sample_id for sample_id in calibrator_list
            if by_id.get(sample_id) is not None and not _truthy(by_id[sample_id].get("Eligible"))
        ]
        valid_calibrator_values = [
            float(by_id[sample_id]["Value"])
            for sample_id in calibrator_list
            if by_id.get(sample_id) is not None
            and _truthy(by_id[sample_id].get("Eligible"))
            and np.isfinite(by_id[sample_id].get("Value", np.nan))
            and float(by_id[sample_id]["Value"]) > 0
        ]
        calibrator_mean = float(np.mean(valid_calibrator_values)) if valid_calibrator_values else np.nan
        complete_calibrators = bool(calibrator_list) and not missing_calibrators and not invalid_calibrators and len(valid_calibrator_values) == len(calibrator_list)
        if not calibrator_list:
            calibrator_status = "No calibrator SampleIds selected"
        elif missing_calibrators or invalid_calibrators:
            affected = list(dict.fromkeys(missing_calibrators + invalid_calibrators))
            calibrator_status = "Selected calibrator ratio(s) unavailable: {}".format(", ".join(affected))
        elif not valid_calibrator_values:
            calibrator_status = "No eligible calibrator ratios"
        else:
            calibrator_status = "Ready"
        if calibrator_status != "Ready":
            _add_issue(issues, calibrator_status, gene=target_gene)
        for row in gene_rows:
            fold = dict(row)
            fold["Endpoint"] = "Calibrator-relative fold change"
            fold["CalibratorMeanValue"] = calibrator_mean
            fold["CalibratorN"] = len(valid_calibrator_values)
            fold["Calibrators"] = ", ".join(calibrator_list)
            if complete_calibrators and np.isfinite(calibrator_mean) and calibrator_mean > 0 and _truthy(row.get("Eligible")) and np.isfinite(row.get("Value", np.nan)):
                fold["Value"] = float(row["Value"]) / calibrator_mean
                fold["Eligible"] = bool(np.isfinite(fold["Value"]))
                if fold["Status"] == "Ready":
                    fold["Status"] = "Ready"
                else:
                    fold["Status"] = _append_message(fold["Status"], "Calibrator-relative value is provisional")
            else:
                fold["Value"] = np.nan
                fold["Eligible"] = False
                fold["Status"] = _append_message(fold.get("Status", ""), calibrator_status)
            standard_results.append(fold)

    # Compute sample-level ΔCq values from biological sample means, never wells.
    delta_by_gene: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for target_gene in target_genes:
        delta_by_gene[target_gene] = {}
        for sample_id in sample_ids:
            target_summary = summary_rows.get((sample_id, target_gene))
            ref_summaries = [summary_rows.get((sample_id, reference)) for reference in references_list]
            missing_refs = [reference for reference, row in zip(references_list, ref_summaries) if row is None]
            bad_refs = [
                reference for reference, row in zip(references_list, ref_summaries)
                if row is not None and (not _truthy(row.get("Eligible")) or not np.isfinite(row.get("MeanCq", np.nan)))
            ]
            record: Dict[str, Any] = {
                "SampleId": sample_id,
                "Gene": target_gene,
                "Target": target_summary,
                "Refs": ref_summaries,
                "MissingRefs": missing_refs,
                "BadRefs": bad_refs,
                "DeltaCq": np.nan,
                "ReferenceMeanCq": np.nan,
                "Eligible": False,
                "Status": "",
                "RangeStatus": "Unavailable",
                "PlateTrace": "",
            }
            range_rows = [quantity_rows.get((sample_id, target_gene))] + [quantity_rows.get((sample_id, reference)) for reference in references_list]
            range_details = [
                "{}:{}".format(_identifier(row.get("Gene")) or "unknown", row.get("RangeStatus"))
                for row in range_rows if row is not None and row.get("RangeStatus") in {"AboveRange", "BelowRange", "RangeMismatch"}
            ]
            known_out_of_range = [
                row for row in range_rows if row is not None and row.get("RangeStatus") in {"AboveRange", "BelowRange"}
            ]
            range_mismatch = any(row is not None and row.get("RangeStatus") == "RangeMismatch" for row in range_rows)
            if range_details:
                record["RangeStatus"] = "; ".join(range_details)
            if target_summary is None:
                record["Status"] = "Target measurement missing"
            elif not _truthy(target_summary.get("Eligible")) or not np.isfinite(target_summary.get("MeanCq", np.nan)):
                record["Status"] = "Target Cq unavailable or requires review: {}".format(target_summary.get("Status", "unknown"))
            elif missing_refs:
                record["Status"] = "Missing reference gene(s): {}".format(", ".join(missing_refs))
            elif bad_refs:
                record["Status"] = "Reference Cq unavailable or requires review: {}".format(", ".join(bad_refs))
            elif not references_list:
                record["Status"] = "No reference genes selected"
            else:
                ref_cqs = [float(row["MeanCq"]) for row in ref_summaries if row is not None]
                record["ReferenceMeanCq"] = float(np.mean(ref_cqs))
                record["DeltaCq"] = float(target_summary["MeanCq"]) - record["ReferenceMeanCq"]
                if range_mismatch:
                    record["Status"] = "Observed Cq and fitted concentration ranges disagree; curve review required"
                elif known_out_of_range and not allow_extrapolation:
                    record["Status"] = "Outside calibrated range; final ΔΔCq result withheld pending review"
                else:
                    record["Eligible"] = True
                    if known_out_of_range:
                        record["Status"] = "Provisional outside calibrated range: {}".format(", ".join(range_details))
                    else:
                        record["Status"] = "Ready"
            trace_parts = []
            if target_summary is not None:
                trace_parts.append("Target {}:{}".format(target_gene, target_summary.get("PlateTrace", "unknown")))
            for reference, row in zip(references_list, ref_summaries):
                if row is not None:
                    trace_parts.append("{}:{}".format(reference, row.get("PlateTrace", "unknown")))
            record["PlateTrace"] = " | ".join(trace_parts)
            if record["Status"] and record["Status"] != "Ready":
                _add_issue(issues, record["Status"], gene=target_gene, sample_id=sample_id)
            delta_by_gene[target_gene][sample_id] = record

    for target_gene in target_genes:
        gene_records = delta_by_gene[target_gene]
        selected_calibrators = [gene_records.get(sample_id) for sample_id in calibrator_list]
        missing_calibrators = [sample_id for sample_id, record in zip(calibrator_list, selected_calibrators) if record is None]
        invalid_calibrators = [
            sample_id for sample_id, record in zip(calibrator_list, selected_calibrators)
            if record is not None and not record["Eligible"]
        ]
        valid_values = [record["DeltaCq"] for record in selected_calibrators if record is not None and record["Eligible"]]
        calibrator_mean = float(np.mean(valid_values)) if valid_values else np.nan
        if not calibrator_list:
            calibrator_status = "No calibrator SampleIds selected"
        elif not valid_values:
            calibrator_status = "No eligible calibrator samples"
        elif missing_calibrators or invalid_calibrators:
            affected = list(dict.fromkeys(missing_calibrators + invalid_calibrators))
            calibrator_status = "Selected calibrator(s) need review: {}".format(", ".join(affected))
        else:
            calibrator_status = "Ready"
        calibrator_provisional = any(
            record is not None and "provisional outside calibrated range" in str(record.get("Status", "")).casefold()
            for record in selected_calibrators
        )
        if calibrator_status != "Ready":
            _add_issue(issues, calibrator_status, gene=target_gene)

        for sample_id, record in gene_records.items():
            target_summary = record["Target"]
            value = np.nan
            eligible = bool(record["Eligible"] and np.isfinite(calibrator_mean) and calibrator_status == "Ready")
            status_parts = []
            if not record["Eligible"]:
                status_parts.append(record["Status"])
            if not np.isfinite(calibrator_mean):
                status_parts.append(calibrator_status)
            elif calibrator_status != "Ready":
                status_parts.append(calibrator_status)
            if record["Eligible"] and np.isfinite(calibrator_mean):
                dd_cq = float(record["DeltaCq"] - calibrator_mean)
                value = float(np.exp2(-dd_cq))
                if not np.isfinite(value):
                    value = np.nan
                    eligible = False
                    status_parts.append("Calculated fold change is not finite")
            if not efficiency_confirmed:
                eligible = False
                status_parts.append("Exploratory: 2^-ΔΔCq efficiency assumptions not confirmed")
            if calibrator_provisional and calibrator_status == "Ready":
                status_parts.append("Provisional: selected calibrator is outside calibrated range")
            elif "provisional outside calibrated range" in str(record.get("Status", "")).casefold():
                status_parts.append(record["Status"])
            elif record["Eligible"] and calibrator_status == "Ready" and np.isfinite(value):
                status_parts.append("Ready")
            standard_results_meta = target_summary
            if standard_results_meta is None:
                standard_results_meta = {}
            ddcq_results.append({
                "Endpoint": "Calibrator-relative fold change",
                "SampleId": sample_id,
                "Label": standard_results_meta.get("Label"),
                "Group": standard_results_meta.get("Group"),
                "PairId": standard_results_meta.get("PairId"),
                "Gene": target_gene,
                "Value": value,
                "Eligible": bool(eligible),
                "Status": "; ".join(dict.fromkeys(status_parts)) if status_parts else "Unavailable",
                "Plate": standard_results_meta.get("Plate"),
                "PlateTrace": record["PlateTrace"],
                "DeltaCq": record["DeltaCq"],
                "ReferenceMeanCq": record["ReferenceMeanCq"],
                "RangeStatus": record["RangeStatus"],
                "CalibratorMeanDeltaCq": calibrator_mean,
                "CalibratorN": len(valid_values),
                "Calibrators": ", ".join(calibrator_list),
                "EfficiencyConfirmed": bool(efficiency_confirmed),
            })

    standard_df = pd.DataFrame(standard_results)
    ddcq_df = pd.DataFrame(ddcq_results)
    if standard_df.empty:
        standard_df = pd.DataFrame(columns=["Endpoint"] + RESULT_COLUMNS + ["TargetQuantity", "ReferenceQuantity", "ReferenceGenes", "RangeStatus", "CalibratorMeanValue", "CalibratorN", "Calibrators"])
    if ddcq_df.empty:
        ddcq_df = pd.DataFrame(columns=["Endpoint"] + RESULT_COLUMNS + ["DeltaCq", "ReferenceMeanCq", "CalibratorMeanDeltaCq", "CalibratorN", "Calibrators", "EfficiencyConfirmed"])
    issue_df = pd.DataFrame(issues, columns=ISSUE_COLUMNS).drop_duplicates().reset_index(drop=True)
    return {
        "samples": summary,
        "curves": curves,
        "standard_points": standard_points,
        "quantities": quantities,
        "standard_results": standard_df,
        "ddcq_results": ddcq_df,
        "issues": issue_df,
    }
