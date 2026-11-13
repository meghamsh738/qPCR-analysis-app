"""Create an editable CYC1-style qPCR review workbook.

``cyc1_review_workbook(wells, standards, results, settings)`` returns XLSX
bytes. ``wells`` should be the mapped/prepared well DataFrame,
``standards`` the editable plate/gene/label concentration map, and ``results``
the DataFrame dictionary returned by :mod:`qpcr_workflow.analyze`.

The workbook uses live Excel formulas for technical-replicate means, standard
curve fits, sample quantities, and target/reference ratios. Formula results
are cached at export so the workbook opens with values populated; Excel or a
compatible spreadsheet app recalculates them after edits. XlsxWriter does not
calculate formulas itself. The raw Results block is retained for provenance,
but it is not a live input table: edit the visible Cq replicate cells and
standard Include cells to review the workbook, then re-import reviewed data
and settings into the app to update the app's analysis state. Sample-level
review eligibility/reason flags come from the analysis result and are retained
as review gates; editing them in Excel is not a supported workflow.

Curves are always plate-scoped. Samples with ambiguous multi-plate mappings
remain ineligible for normalization. Out-of-range quantities are blank unless
``settings['allow_extrapolation']`` is true, in which case the displayed
quantity and normalized result remain labelled provisional. This workbook is
a review aid, not a replacement for the original instrument export.
"""

from __future__ import annotations

import io
import math
import re
from collections import defaultdict
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import xlsxwriter
from xlsxwriter.utility import xl_col_to_name, xl_range


def _blank(value: Any) -> bool:
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _text(value: Any) -> str:
    return "" if _blank(value) else str(value).strip()


def _truthy(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)) and not _blank(value):
        return float(value) != 0
    return _text(value).casefold() in {"yes", "y", "true", "1", "keep", "included"}


def _finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError, OverflowError):
        return False


def _safe_sheet_name(name: Any) -> str:
    clean = re.sub(r"[\\/*?:\[\]]", " ", _text(name)).strip().lower()
    clean = re.sub(r"\s+", " ", clean) or "assay"
    return clean[:31]


def _quote_sheet(name: str) -> str:
    return "'{}'".format(name.replace("'", "''"))


def _write_literal(sheet: Any, row: int, col: int, value: Any, fmt: Any = None) -> None:
    """Write source text literally so imported text cannot become a formula."""
    if _blank(value):
        sheet.write_blank(row, col, None, fmt)
    elif isinstance(value, (bool, np.bool_)):
        sheet.write_boolean(row, col, bool(value), fmt)
    elif isinstance(value, (int, float, np.integer, np.floating)) and _finite(value):
        sheet.write_number(row, col, float(value), fmt)
    else:
        sheet.write_string(row, col, str(value), fmt)


def _frame(results: Mapping[str, Any], name: str) -> pd.DataFrame:
    value = results.get(name)
    return value.copy(deep=False) if isinstance(value, pd.DataFrame) else pd.DataFrame()


def _initial_fit(
    standard_rows: Sequence[Mapping[str, Any]],
) -> Dict[str, Any]:
    """Calculate cache values for formulas; Excel remains the live calculator."""
    selected: List[Tuple[float, float]] = []
    missing_map = False
    duplicate_labels = False
    missing_reason = False
    missing_reading = False
    labels: List[str] = []
    for item in standard_rows:
        label = _text(item.get("Label"))
        quantity = item.get("Concentration")
        include = _truthy(item.get("Include"))
        labels.append(label)
        if not label or not item.get("mapped", True):
            missing_map = True
        if include:
            readings = item.get("readings", [])
            if not readings or item.get("blocking_reading_status", False) or item.get("duplicate_well", False):
                missing_reading = True
            if not _finite(quantity) or float(quantity) <= 0:
                missing_map = True
            elif readings and not item.get("duplicate_well", False):
                selected.append((math.log10(float(quantity)), float(np.mean(readings))))
        elif item.get("mapped", True) and not _text(item.get("explicit_reason")):
            missing_reason = True
        if item.get("missing_well_exclusion_reason", False):
            missing_reason = True
    if len(labels) != len(set(labels)):
        duplicate_labels = True
    if missing_map or duplicate_labels:
        return {"status": "Duplicate standard mapping" if duplicate_labels else "Selected standard mapping or reading requires review"}
    if missing_reason:
        return {"status": "Selected standard mapping or reading requires review"}
    if missing_reading:
        return {"status": "Selected standard mapping or reading requires review"}
    if len(selected) < 2:
        return {"status": "Insufficient distinct standard levels"}
    concentrations = [float(item["Concentration"]) for item in standard_rows if _truthy(item.get("Include"))
                      and _finite(item.get("Concentration")) and float(item["Concentration"]) > 0
                      and item.get("readings") and not item.get("duplicate_well", False)]
    if len(concentrations) != len(set(concentrations)):
        return {"status": "Duplicate standard concentration"}
    x = np.array([point[0] for point in selected], dtype=float)
    y = np.array([point[1] for point in selected], dtype=float)
    slope, intercept = np.linalg.lstsq(np.column_stack([x, np.ones(len(x))]), y, rcond=None)[0]
    predicted = slope * x + intercept
    total = float(np.sum((y - np.mean(y)) ** 2))
    residual = float(np.sum((y - predicted) ** 2))
    r2 = float(1.0 - residual / total) if total > 0 else math.nan
    if not math.isfinite(float(slope)) or slope >= 0:
        return {"status": "Invalid non-negative slope"}
    if not math.isfinite(r2):
        return {"status": "Invalid standard fit"}
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        efficiency = float(np.power(10.0, -1.0 / float(slope)) - 1.0) if math.isfinite(float(slope)) else math.nan
    if not math.isfinite(efficiency):
        return {"status": "Invalid efficiency calculation"}
    return {
        "status": "Ready",
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": r2,
        "efficiency": efficiency,
        "min_cq": float(np.min(y)),
        "max_cq": float(np.max(y)),
        "min_quantity": float(min(10 ** point[0] for point in selected)),
        "max_quantity": float(max(10 ** point[0] for point in selected)),
    }


def _assay_groups(wells: pd.DataFrame) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    seen = set()
    for _, row in wells.iterrows():
        gene, plate = _text(row.get("Gene")), _text(row.get("Plate"))
        if not gene or not plate:
            continue
        key = (gene, plate)
        if key not in seen:
            seen.add(key)
            out.append(key)
    return out


def _numeric_retained(group: pd.DataFrame) -> List[float]:
    if group.empty:
        return []
    if "Keep" in group.columns:
        group = group[group["Keep"].map(_truthy)]
    if "CqStatus" in group.columns:
        group = group[group["CqStatus"].astype(str).eq("Detected")]
    values = []
    for value in group.get("Cq", pd.Series(dtype=object)).tolist():
        if _finite(value):
            values.append(float(value))
    return values


def _status_cell(status: Any) -> str:
    value = _text(status)
    return value or "Ready"


def _authoritative_fit_cache(
    curves: pd.DataFrame,
    gene: str,
    plate: str,
    fallback: Dict[str, Any],
) -> Dict[str, Any]:
    """Use the analysis pipeline's fit as the initial formula-cache value."""
    if curves.empty or not {"Gene", "Plate", "Eligible", "Status"}.issubset(curves.columns):
        return fallback
    rows = curves[curves["Gene"].map(_text).eq(gene) & curves["Plate"].map(_text).eq(plate)]
    if len(rows) != 1:
        return fallback
    row = rows.iloc[0]
    cache: Dict[str, Any] = {"status": _status_cell(row.get("Status"))}
    if not _truthy(row.get("Eligible")):
        return cache
    required = {
        "slope": row.get("Slope"),
        "intercept": row.get("Intercept"),
        "r2": row.get("R2"),
        "efficiency": float(row.get("EfficiencyPercent")) / 100.0 if _finite(row.get("EfficiencyPercent")) else math.nan,
        "min_cq": row.get("MinCq"),
        "max_cq": row.get("MaxCq"),
        "min_quantity": row.get("MinConcentration"),
        "max_quantity": row.get("MaxConcentration"),
    }
    if not all(_finite(value) for value in required.values()):
        return fallback
    cache.update({key: float(value) for key, value in required.items()})
    cache["status"] = "Ready"
    return cache


def cyc1_review_workbook(
    wells: pd.DataFrame,
    standards: pd.DataFrame,
    results: Dict[str, pd.DataFrame],
    settings: Dict[str, Any],
) -> bytes:
    """Return an editable CYC1-style qPCR review workbook as XLSX bytes.

    Expected inputs are the mapped/prepared well DataFrame, standard map,
    ``analyze()`` result tables, and analysis settings. The file contains one
    lower-case assay sheet per gene/plate, a compact ``normalised`` sheet, and
    original source/result fields in a raw table beneath each assay review.

    Limitations: Excel formula caches are calculated from the supplied result
    DataFrames and do not update until a spreadsheet engine recalculates.
    XlsxWriter only writes formulas; it does not evaluate them. Visible Cq and
    standard-inclusion edits recalculate in Excel, but raw table edits,
    sample-level eligibility decisions, mappings, reference selections, and
    settings must be re-imported into the app to regenerate the workbook.
    The exporter retains the supplied sample review eligibility/status as a
    gate, does not infer absent biological metadata, and never pools standard
    curves across plates. It does not perform the app's separate ΔΔCq endpoint.
    """
    if not isinstance(wells, pd.DataFrame):
        raise TypeError("wells must be a pandas DataFrame")
    if not isinstance(standards, pd.DataFrame):
        raise TypeError("standards must be a pandas DataFrame")
    if not isinstance(results, Mapping):
        raise TypeError("results must be the result-table mapping from analyze()")
    settings = dict(settings or {})

    # analyze() normally receives prepared wells. Prepare raw mapped inputs if
    # this API is called independently, preserving source fields in the block.
    from qpcr_workflow import prepare_wells

    source_frame = wells.copy(deep=True).reset_index(drop=True)
    prepared = source_frame.copy(deep=True)
    if not {"CqStatus", "IdentityStatus", "OriginalCq", "SampleId"}.issubset(prepared.columns):
        prepared = prepare_wells(prepared)
    if "Type" not in prepared.columns:
        prepared["Type"] = "Sample"

    samples = _frame(results, "samples")
    curves = _frame(results, "curves")
    standard_results = _frame(results, "standard_results")
    allow_extrapolation = bool(settings.get("allow_extrapolation", False))
    references = [str(value).strip() for value in settings.get("references", []) if not _blank(value)]
    ratio_results = standard_results
    if not ratio_results.empty and "Endpoint" in ratio_results.columns:
        ratio_results = ratio_results[ratio_results["Endpoint"].astype(str).eq("Target/reference ratio")]

    groups = _assay_groups(prepared)
    if not groups:
        # An empty/malformed import should still preserve all source rows.
        groups = [("assay", "plate")]
    plate_counts: Dict[str, int] = defaultdict(int)
    for gene, _plate in groups:
        plate_counts[gene.casefold()] += 1

    used_names = {"normalised", "raw results"}
    sheet_for_group: Dict[Tuple[str, str], str] = {}
    for gene, plate in groups:
        base = _safe_sheet_name(gene if plate_counts[gene.casefold()] == 1 else "{} - {}".format(gene, plate))
        sheet_name = base
        suffix = 2
        while sheet_name in used_names:
            ending = " {}".format(suffix)
            sheet_name = base[:31 - len(ending)] + ending
            suffix += 1
        used_names.add(sheet_name)
        sheet_for_group[(gene, plate)] = sheet_name

    output = io.BytesIO()
    workbook = xlsxwriter.Workbook(output, {
        "in_memory": True,
        "strings_to_formulas": False,
        "strings_to_urls": False,
    })
    workbook.set_calc_mode("auto")
    workbook.set_properties({
        "title": "qPCR CYC1-style review workbook",
        "subject": "Editable plate-scoped review export",
        "author": "EasyLab Next",
        "comments": "Formula results are cached at export and recalculate in a spreadsheet application.",
    })
    workbook.set_custom_property("EasyLabWorkbookKind", "qpcr-review", "text")
    formats = {
        "base": workbook.add_format({"font_name": "Aptos", "font_size": 10, "font_color": "#203040", "valign": "vcenter"}),
        "header": workbook.add_format({"font_name": "Aptos", "font_size": 10, "bold": True, "font_color": "#FFFFFF", "bg_color": "#17365D", "border": 1, "border_color": "#17365D", "text_wrap": True, "valign": "vcenter"}),
        "cell": workbook.add_format({"font_name": "Aptos", "font_size": 10, "font_color": "#203040", "border": 1, "border_color": "#D7E2F1", "valign": "vcenter"}),
        "text": workbook.add_format({"font_name": "Aptos", "font_size": 10, "font_color": "#203040", "border": 1, "border_color": "#D7E2F1", "valign": "vcenter", "num_format": "@"}),
        "cq": workbook.add_format({"font_name": "Aptos", "font_size": 10, "font_color": "#203040", "border": 1, "border_color": "#D7E2F1", "num_format": "0.000"}),
        "quantity": workbook.add_format({"font_name": "Aptos", "font_size": 10, "font_color": "#203040", "border": 1, "border_color": "#D7E2F1", "num_format": "0.000"}),
        "ratio": workbook.add_format({"font_name": "Aptos", "font_size": 10, "font_color": "#203040", "border": 1, "border_color": "#D7E2F1", "num_format": "0.000E+00"}),
        "log": workbook.add_format({"font_name": "Aptos", "font_size": 10, "font_color": "#203040", "border": 1, "border_color": "#D7E2F1", "num_format": "0.0000"}),
        "slope": workbook.add_format({"font_name": "Aptos", "font_size": 10, "font_color": "#203040", "border": 1, "border_color": "#D7E2F1", "num_format": "0.0000"}),
        "r2": workbook.add_format({"font_name": "Aptos", "font_size": 10, "font_color": "#203040", "border": 1, "border_color": "#D7E2F1", "num_format": "0.0000"}),
        "efficiency": workbook.add_format({"font_name": "Aptos", "font_size": 10, "font_color": "#203040", "border": 1, "border_color": "#D7E2F1", "num_format": "0.0%"}),
        "label": workbook.add_format({"font_name": "Aptos", "font_size": 10, "bold": True, "font_color": "#203040", "border": 1, "border_color": "#D7E2F1", "valign": "vcenter"}),
        "section": workbook.add_format({"font_name": "Aptos", "font_size": 10, "bold": True, "font_color": "#17365D", "bottom": 1, "bottom_color": "#4472C4"}),
        "warning": workbook.add_format({"font_name": "Aptos", "font_size": 10, "font_color": "#7F6000", "bg_color": "#FFF2CC", "border": 1, "border_color": "#D7E2F1", "text_wrap": True}),
        "raw_header": workbook.add_format({"font_name": "Aptos", "font_size": 10, "bold": True, "font_color": "#FFFFFF", "bg_color": "#44546A", "border": 1, "border_color": "#44546A", "text_wrap": True}),
    }

    assay_details: Dict[Tuple[str, str], Dict[str, Any]] = {}

    for group_index, (gene, plate) in enumerate(groups):
        sheet_name = sheet_for_group[(gene, plate)]
        sheet = workbook.add_worksheet(sheet_name)
        sheet.set_default_row(18)
        sheet.freeze_panes(1, 0)
        sheet.set_zoom(90)
        sheet.set_landscape()
        sheet.set_paper(8)  # A3 keeps the broad review sheet legible when printed.
        sheet.fit_to_pages(1, 0)
        sheet.set_margins(0.25, 0.25, 0.35, 0.35)
        sheet.center_horizontally()
        group_wells = prepared[
            prepared["Gene"].map(_text).eq(gene)
            & prepared["Plate"].map(_text).eq(plate)
        ].copy()
        sample_wells = group_wells[group_wells["Type"].astype(str).str.casefold().eq("sample")]
        std_wells = group_wells[group_wells["Type"].astype(str).str.casefold().eq("standard")]

        sample_records: List[Dict[str, Any]] = []
        seen_ids = set()
        for _, well in sample_wells.iterrows():
            sample_id = _text(well.get("SampleId"))
            if not sample_id or sample_id in seen_ids:
                continue
            seen_ids.add(sample_id)
            rows = sample_wells[sample_wells["SampleId"].map(_text).eq(sample_id)]
            readings = _numeric_retained(rows)
            summary_rows = samples[
                samples.get("SampleId", pd.Series(dtype=object)).map(_text).eq(sample_id)
                & samples.get("Gene", pd.Series(dtype=object)).map(_text).eq(gene)
            ] if not samples.empty else pd.DataFrame()
            summary = summary_rows.iloc[0].to_dict() if not summary_rows.empty else {}
            record = {
                "SampleId": sample_id,
                "Label": well.get("Label"),
                "Group": well.get("Group"),
                "PairId": well.get("PairId"),
                "Plate": plate,
                "readings": readings,
                "base_eligible": _truthy(summary.get("Eligible", False)),
                "base_status": _status_cell(summary.get("Status", "No sample review result")),
                "mean": float(np.mean(readings)) if readings else None,
            }
            sample_records.append(record)

        max_reps = max([len(item["readings"]) for item in sample_records] + [2])
        replicate_headers = ["Cq duplicate {}".format(index + 1) for index in range(max_reps)]
        visible_headers = ["Sample", "Group", "Pair ID", "Plate", "SampleId"] + replicate_headers + [
            "Cq mean", "{} Power (AU)".format(gene), "QC / inclusion"
        ]
        sample_header_row = 0
        for col, header in enumerate(visible_headers):
            sheet.write_string(sample_header_row, col, header, formats["header"])
        sample_first_row = 1
        sample_last_row = sample_first_row + len(sample_records) - 1
        mean_col = 5 + max_reps
        power_col = mean_col + 1
        status_col = mean_col + 2
        chart_col = max(status_col + 4, 11)
        for index, item in enumerate(sample_records):
            row = sample_first_row + index
            item["sheet_row"] = row
            sheet.write(row, 0, item["Label"], formats["text"])
            sheet.write(row, 1, item["Group"], formats["text"])
            sheet.write(row, 2, item["PairId"], formats["text"])
            sheet.write(row, 3, plate, formats["text"])
            sheet.write_string(row, 4, item["SampleId"], formats["text"])
            for rep in range(max_reps):
                value = item["readings"][rep] if rep < len(item["readings"]) else None
                if value is not None:
                    sheet.write_number(row, 5 + rep, float(value), formats["cq"])
                else:
                    sheet.write_blank(row, 5 + rep, None, formats["cq"])
            excel_row = row + 1
            mean_letter = xl_col_to_name(mean_col)
            cq_start = xl_col_to_name(5)
            cq_end = xl_col_to_name(5 + max_reps - 1)
            mean_value = item["mean"]
            sheet.write_formula(
                row, mean_col,
                '=IF(COUNT({0}{1}:{2}{1})=0,"",AVERAGE({0}{1}:{2}{1}))'.format(cq_start, excel_row, cq_end),
                formats["cq"], mean_value if mean_value is not None else "",
            )
            item["base_eligible"] = bool(item["base_eligible"])
            item["base_status"] = _text(item["base_status"])

        # Standard points are shown in plate order, including observed levels
        # without a concentration mapping and mapped levels without wells.
        mappings = standards[
            standards.get("Gene", pd.Series(dtype=object)).map(_text).eq(gene)
            & standards.get("Plate", pd.Series(dtype=object)).map(_text).eq(plate)
        ].copy() if not standards.empty else pd.DataFrame()
        standard_records: List[Dict[str, Any]] = []
        map_labels = mappings.get("Label", pd.Series(dtype=object)).map(_text).tolist() if not mappings.empty else []
        observed = std_wells.get("Label", pd.Series(dtype=object)).map(_text).drop_duplicates().tolist()
        for _, mapping in mappings.iterrows():
            label = _text(mapping.get("Label"))
            level_wells = std_wells[std_wells["Label"].map(_text).eq(label)]
            retained = level_wells[level_wells.get("Keep", pd.Series(True, index=level_wells.index)).map(_truthy)]
            readings = _numeric_retained(level_wells)
            duplicate_well = bool(level_wells.get("IdentityStatus", pd.Series(index=level_wells.index, dtype=object)).eq("Duplicate well identity").any())
            blocking_reading_status = bool(
                (retained.get("CqStatus", pd.Series(index=retained.index, dtype=object)).astype(str) == "Invalid").any()
            )
            excluded = level_wells[~level_wells.get("Keep", pd.Series(True, index=level_wells.index)).map(_truthy)]
            missing_well_exclusion_reason = bool(
                excluded.get("Reason", pd.Series("", index=excluded.index)).map(_blank).any()
            )
            try:
                concentration = float(mapping.get("Concentration"))
            except (TypeError, ValueError, OverflowError):
                concentration = None
            if concentration is not None and not math.isfinite(concentration):
                concentration = None
            standard_records.append({
                "Label": mapping.get("Label"),
                "Concentration": concentration,
                "Include": _truthy(mapping.get("Include")),
                "readings": readings,
                "retained_count": len(retained),
                "explicit_reason": _text(mapping.get("Reason")),
                "mapped": True,
                "duplicate_well": duplicate_well,
                "blocking_reading_status": blocking_reading_status,
                "missing_well_exclusion_reason": missing_well_exclusion_reason,
                "missing_mapping": not label or concentration is None or concentration <= 0,
            })
        for label in observed:
            if label and label not in map_labels:
                level_wells = std_wells[std_wells["Label"].map(_text).eq(label)]
                retained = level_wells[level_wells.get("Keep", pd.Series(True, index=level_wells.index)).map(_truthy)]
                excluded = level_wells[~level_wells.get("Keep", pd.Series(True, index=level_wells.index)).map(_truthy)]
                standard_records.append({
                    "Label": label, "Concentration": None, "Include": False,
                    "readings": _numeric_retained(level_wells), "retained_count": len(retained),
                    "explicit_reason": "", "mapped": False, "duplicate_well": False,
                    "blocking_reading_status": False,
                    "missing_well_exclusion_reason": bool(excluded.get("Reason", pd.Series("", index=excluded.index)).map(_blank).any()),
                    "missing_mapping": True,
                })
        # Preserve all retained replicate positions in the visible panel. If a
        # kept reading is missing/invalid, its cell stays blank and the fit QC
        # helper sees the expected retained count.
        max_std_reps = max([len(item["readings"]) for item in standard_records] + [2])
        std_headers = ["Standard", "Quantity (AU)", "log10 quantity"] + [
            "Cq duplicate {}".format(index + 1) for index in range(max_std_reps)
        ] + ["Cq mean", "Used in curve", "Reason"]
        std_mean_col = 3 + max_std_reps
        include_col = std_mean_col + 1
        reason_col = include_col + 1
        std_start_col = chart_col
        chart_row = 1
        for offset, header in enumerate(std_headers):
            sheet.write_string(19, std_start_col + offset, header, formats["header"])
        std_first_row = 20
        std_last_row = std_first_row + len(standard_records) - 1
        std_helper_start = max(std_start_col + len(std_headers), len(source_frame.columns) + 1, chart_col + 10)
        std_x_col, std_y_col, std_conc_col = std_helper_start, std_helper_start + 1, std_helper_start + 2
        std_missing_col, std_badmap_col, std_dupconc_col = std_helper_start + 3, std_helper_start + 4, std_helper_start + 5
        std_reason_col, std_dupmap_col = std_helper_start + 6, std_helper_start + 7
        initial_standard_rows = []
        for index, item in enumerate(standard_records):
            row = std_first_row + index
            item["sheet_row"] = row
            xlrow = row + 1
            label_col = std_start_col
            conc_col = std_start_col + 1
            log_col = std_start_col + 2
            reps_col = std_start_col + 3
            sheet.write(row, label_col, item["Label"], formats["text"])
            if _finite(item["Concentration"]) and float(item["Concentration"]) > 0:
                sheet.write_number(row, conc_col, float(item["Concentration"]), formats["quantity"])
                conc_cache = float(item["Concentration"])
                log_cache: Any = math.log10(conc_cache)
                log_formula = '=IF(AND(ISNUMBER({0}{1}),{0}{1}>0),LOG10({0}{1}),"")'.format(xl_col_to_name(conc_col), xlrow)
                sheet.write_formula(row, log_col, log_formula, formats["log"], log_cache)
            else:
                sheet.write_blank(row, conc_col, None, formats["quantity"])
                sheet.write_formula(row, log_col, '=IF(AND(ISNUMBER({0}{1}),{0}{1}>0),LOG10({0}{1}),"")'.format(xl_col_to_name(conc_col), xlrow), formats["log"], "")
                conc_cache = None
            for rep in range(max_std_reps):
                value = item["readings"][rep] if rep < len(item["readings"]) else None
                if value is None:
                    sheet.write_blank(row, reps_col + rep, None, formats["cq"])
                else:
                    sheet.write_number(row, reps_col + rep, float(value), formats["cq"])
            rep_start_letter = xl_col_to_name(reps_col)
            rep_end_letter = xl_col_to_name(reps_col + max_std_reps - 1)
            mean_formula = '=IF(COUNT({0}{1}:{2}{1})=0,"",AVERAGE({0}{1}:{2}{1}))'.format(rep_start_letter, xlrow, rep_end_letter)
            mean_val: Any = float(np.mean(item["readings"])) if item["readings"] else ""
            sheet.write_formula(row, std_start_col + std_mean_col, mean_formula, formats["cq"], mean_val)
            used_val = "Yes" if item["Include"] else "No"
            sheet.write_string(row, std_start_col + include_col, used_val, formats["text"])
            reason_text = _text(item.get("explicit_reason"))
            reason_warning = (item.get("mapped", True) and not item["Include"] and not reason_text) or item.get("missing_well_exclusion_reason", False)
            sheet.write_string(row, std_start_col + reason_col, reason_text,
                               formats["warning"] if reason_warning else formats["text"])
            initial_standard_rows.append({**item, "Label": _text(item["Label"]), "Concentration": conc_cache,
                                          "Include": item["Include"], "readings": item["readings"]})
        if standard_records:
            sheet.data_validation(std_first_row, std_start_col + include_col, std_last_row, std_start_col + include_col,
                                  {"validate": "list", "source": ["Yes", "No"]})
        sheet.set_column(std_start_col, std_start_col, 18)
        sheet.set_column(std_start_col + 1, std_start_col + 1, 15)
        sheet.set_column(std_start_col + 2, std_start_col + 2, 15)
        sheet.set_column(std_start_col + 3, std_start_col + 3 + max_std_reps - 1, 13)
        sheet.set_column(std_start_col + std_mean_col, std_start_col + std_mean_col, 13)
        sheet.set_column(std_start_col + include_col, std_start_col + include_col, 15)
        sheet.set_column(std_start_col + reason_col, std_start_col + reason_col, 38)

        fit_cache = _authoritative_fit_cache(
            curves, gene, plate, _initial_fit(initial_standard_rows)
        )
        # Hidden formula helpers filter the full standards panel by its live
        # Include cells. SLOPE, INTERCEPT, RSQ, the chart, and range bounds all
        # consume these same selected point columns.
        std_label_range = xl_range(std_first_row, std_start_col, max(std_last_row, std_first_row), std_start_col)
        std_conc_range = xl_range(std_first_row, std_start_col + 1, max(std_last_row, std_first_row), std_start_col + 1)
        std_include_range = xl_range(std_first_row, std_start_col + include_col, max(std_last_row, std_first_row), std_start_col + include_col)
        fit_x_range = xl_range(std_first_row, std_x_col, max(std_last_row, std_first_row), std_x_col)
        fit_y_range = xl_range(std_first_row, std_y_col, max(std_last_row, std_first_row), std_y_col)
        fit_conc_range = xl_range(std_first_row, std_conc_col, max(std_last_row, std_first_row), std_conc_col)
        missing_range = xl_range(std_first_row, std_missing_col, max(std_last_row, std_first_row), std_missing_col)
        badmap_range = xl_range(std_first_row, std_badmap_col, max(std_last_row, std_first_row), std_badmap_col)
        dupconc_range = xl_range(std_first_row, std_dupconc_col, max(std_last_row, std_first_row), std_dupconc_col)
        reason_range = xl_range(std_first_row, std_reason_col, max(std_last_row, std_first_row), std_reason_col)
        dupmap_range = xl_range(std_first_row, std_dupmap_col, max(std_last_row, std_first_row), std_dupmap_col)
        hidden_cols = [std_x_col, std_y_col, std_conc_col, std_missing_col, std_badmap_col, std_dupconc_col, std_reason_col]
        hidden_titles = ["Selected log10 quantity", "Selected Cq mean", "Selected quantity", "Selected reading issue", "Mapping issue", "Duplicate selected concentration", "Missing exclusion reason"]
        hidden_cols.append(std_dupmap_col)
        hidden_titles.append("Duplicate standard mapping")
        for col, title in zip(hidden_cols, hidden_titles):
            sheet.write(19, col, title)
            sheet.set_column(col, col, None, None, {"hidden": True})
        main_helper_start = max(len(source_frame.columns) + 1, std_dupmap_col + 1, chart_col + 10)
        base_ok_col, eligible_col = main_helper_start, main_helper_start + 1
        sheet.set_column(base_ok_col, eligible_col, None, None, {"hidden": True})
        for item in sample_records:
            sheet.write_boolean(int(item["sheet_row"]), base_ok_col, item["base_eligible"])
        for index, item in enumerate(standard_records):
            row = std_first_row + index
            xlrow = row + 1
            inc_cell = "{}{}".format(xl_col_to_name(std_start_col + include_col), xlrow)
            qty_cell = "{}{}".format(xl_col_to_name(std_start_col + 1), xlrow)
            log_cell = "{}{}".format(xl_col_to_name(std_start_col + 2), xlrow)
            mean_cell = "{}{}".format(xl_col_to_name(std_start_col + std_mean_col), xlrow)
            rep_start = xl_col_to_name(std_start_col + 3)
            rep_end = xl_col_to_name(std_start_col + 3 + max_std_reps - 1)
            duplicate_well = bool(item.get("duplicate_well", False))
            helper_x = '=IF(AND({0}="Yes",ISNUMBER({1}),{1}>0,ISNUMBER({2}),{3}=0),{4},"")'.format(inc_cell, qty_cell, mean_cell, int(duplicate_well), log_cell)
            helper_y = '=IF(AND({0}="Yes",ISNUMBER({1}),{1}>0,ISNUMBER({2}),{3}=0),{2},"")'.format(inc_cell, qty_cell, mean_cell, int(duplicate_well))
            helper_c = '=IF(AND({0}="Yes",ISNUMBER({1}),{1}>0,{2}=0),{1},"")'.format(inc_cell, qty_cell, int(duplicate_well))
            blocker = int(bool(item.get("blocking_reading_status", False) or item.get("duplicate_well", False)))
            helper_m = '=IF({0}="Yes",IF(OR(COUNT({1}{2}:{3}{2})=0,{4}=1),1,0),0)'.format(
                inc_cell, rep_start, xlrow, rep_end, blocker
            )
            mapped = 1 if item.get("mapped", True) else 0
            helper_bad = '=IF(OR({0}="",{1}=0,AND({2}="Yes",OR(NOT(ISNUMBER({3})),{3}<=0))),1,0)'.format(
                "{}{}".format(xl_col_to_name(std_start_col), xlrow), mapped, inc_cell, qty_cell
            )
            helper_dup = '=IF(AND({0}="Yes",ISNUMBER({1}),COUNTIFS({2},{1},{3},"Yes")>1),1,0)'.format(
                inc_cell, qty_cell, std_conc_range, std_include_range
            )
            reason_cell = "{}{}".format(xl_col_to_name(std_start_col + reason_col), xlrow)
            missing_well_reason = 1 if item.get("missing_well_exclusion_reason", False) else 0
            helper_reason = '=IF(OR({0}=1,AND({1}=1,{2}="No",{3}="")),1,0)'.format(
                missing_well_reason, mapped, inc_cell, reason_cell
            )
            helper_dupmap = '=IF(COUNTIF({0},{1})>1,1,0)'.format(
                std_label_range, "{}{}".format(xl_col_to_name(std_start_col), xlrow)
            )
            # cached selected points mirror the currently selected map.
            valid_level = bool(item["Include"] and _finite(item.get("Concentration")) and float(item["Concentration"]) > 0 and item["readings"])
            sheet.write_formula(row, std_x_col, helper_x, None, math.log10(float(item["Concentration"])) if valid_level else "")
            sheet.write_formula(row, std_y_col, helper_y, None, float(np.mean(item["readings"])) if valid_level else "")
            sheet.write_formula(row, std_conc_col, helper_c, None, float(item["Concentration"]) if valid_level else "")
            incomplete = bool(item["Include"] and (not item["readings"] or item.get("blocking_reading_status", False) or item.get("duplicate_well", False)))
            sheet.write_formula(row, std_missing_col, helper_m, None, 1 if incomplete else 0)
            bad_mapping = bool(not item.get("Label") or not item.get("mapped", True) or (
                item["Include"] and (not _finite(item.get("Concentration")) or float(item.get("Concentration") or 0) <= 0)
            ))
            sheet.write_formula(row, std_badmap_col, helper_bad, None, 1 if bad_mapping else 0)
            duplicate_conc = bool(item["Include"] and _finite(item.get("Concentration")) and sum(
                1 for other in standard_records if other["Include"] and _finite(other.get("Concentration")) and float(other["Concentration"]) == float(item["Concentration"])
            ) > 1)
            sheet.write_formula(row, std_dupconc_col, helper_dup, None, 1 if duplicate_conc else 0)
            missing_reason = bool(item.get("missing_well_exclusion_reason", False) or (
                item.get("mapped", True) and not item["Include"] and not _text(item.get("explicit_reason"))
            ))
            sheet.write_formula(row, std_reason_col, helper_reason, None, 1 if missing_reason else 0)
            duplicate_mapping = sum(1 for label in map_labels if label == _text(item.get("Label"))) > 1
            sheet.write_formula(row, std_dupmap_col, helper_dupmap, None, 1 if duplicate_mapping else 0)

        # Fit summary next to the standard table, immediately underneath it.
        fit_summary_row = std_first_row + len(standard_records) + 2
        summary_labels = ["Fitted equation", "Slope", "Intercept", "R-squared", "Efficiency", "Blank wells", "Assay QC"]
        summary_col = std_start_col
        summary_value_col = std_start_col + 1
        summary_rows_by_label = {label: fit_summary_row + i for i, label in enumerate(summary_labels)}
        slope_ref = "{}{}".format(xl_col_to_name(summary_value_col), summary_rows_by_label["Slope"] + 1)
        intercept_ref = "{}{}".format(xl_col_to_name(summary_value_col), summary_rows_by_label["Intercept"] + 1)
        qc_excel_ref = "{}{}".format(xl_col_to_name(summary_value_col), summary_rows_by_label["Assay QC"] + 1)
        equation_body = '"Cq = "&TEXT({0},"0.0000")&" × log10(quantity) "&IF({1}<0,"- ","+ ")&TEXT(ABS({1}),"0.0000")'.format(slope_ref, intercept_ref)
        equation = '=IF({0}="Ready",{1},"")'.format(qc_excel_ref, equation_body)
        # QC is formula-driven so Include changes and edited standard readings
        # invalidate/recover the fit and every dependent quantity together.
        qc_formula = '=IF(SUM({0})>0,"Duplicate standard mapping",IF(SUM({1})>0,"Duplicate standard concentration",IF(OR(SUM({2})>0,SUM({3})>0,SUM({4})>0),"Selected standard mapping or reading requires review",IF(COUNT({5})<2,"Insufficient distinct standard levels",IFERROR(IF(SLOPE({6},{5})>=0,"Invalid non-negative slope",IF(NOT(ISNUMBER(RSQ({6},{5}))),"Invalid standard fit",IFERROR(IF(ISNUMBER(10^(-1/SLOPE({6},{5}))-1),"Ready","Invalid efficiency calculation"),"Invalid efficiency calculation"))),"Invalid standard fit")))))'.format(
            dupmap_range, dupconc_range, badmap_range, missing_range, reason_range, fit_x_range, fit_y_range
        )
        qc_cache = _text(fit_cache.get("status", "Standard fit unavailable"))
        fit_formulas = {
            "Fitted equation": (equation, ("Cq = {:.4f} × log10(quantity) {} {:.4f}".format(fit_cache["slope"], "-" if fit_cache["intercept"] < 0 else "+", abs(fit_cache["intercept"])) if qc_cache == "Ready" else ""), formats["text"]),
            "Slope": ('=IF({0}="Ready",IFERROR(SLOPE({1},{2}),""),"")'.format(xl_col_to_name(summary_value_col) + str(summary_rows_by_label["Assay QC"] + 1), fit_y_range, fit_x_range), fit_cache.get("slope", ""), formats["slope"]),
            "Intercept": ('=IF({0}="Ready",IFERROR(INTERCEPT({1},{2}),""),"")'.format(xl_col_to_name(summary_value_col) + str(summary_rows_by_label["Assay QC"] + 1), fit_y_range, fit_x_range), fit_cache.get("intercept", ""), formats["slope"]),
            "R-squared": ('=IF({0}="Ready",IFERROR(RSQ({1},{2}),""),"")'.format(xl_col_to_name(summary_value_col) + str(summary_rows_by_label["Assay QC"] + 1), fit_y_range, fit_x_range), fit_cache.get("r2", ""), formats["r2"]),
            "Efficiency": ('=IF({0}="Ready",IFERROR(10^(-1/{1})-1,""),"")'.format(xl_col_to_name(summary_value_col) + str(summary_rows_by_label["Assay QC"] + 1), slope_ref), fit_cache.get("efficiency", ""), formats["efficiency"]),
            "Blank wells": ("", int(group_wells[
                group_wells["Type"].astype(str).str.casefold().isin(["blank", "negative"])
            ].shape[0]), formats["cell"]),
            "Assay QC": (qc_formula, qc_cache, formats["warning"] if qc_cache != "Ready" else formats["cell"]),
        }
        for label in summary_labels:
            row = summary_rows_by_label[label]
            sheet.write_string(row, summary_col, label, formats["label"])
            formula, cached, fmt = fit_formulas[label]
            if formula:
                sheet.write_formula(row, summary_value_col, formula, fmt, cached)
            else:
                sheet.write_number(row, summary_value_col, cached, fmt)
        # Equation formula must be gated by current QC, just like other fit
        # diagnostics; rewrite the cached result above with an empty cache for
        # invalid initial fits (already prepared in fit_formulas).

        slope_ref = "{}{}".format(xl_col_to_name(summary_value_col), summary_rows_by_label["Slope"] + 1)
        intercept_ref = "{}{}".format(xl_col_to_name(summary_value_col), summary_rows_by_label["Intercept"] + 1)
        min_y_ref = "MIN({})".format(fit_y_range)
        max_y_ref = "MAX({})".format(fit_y_range)
        min_q_ref = "MIN({})".format(fit_conc_range)
        max_q_ref = "MAX({})".format(fit_conc_range)
        # Live sample formulas use the exact fit/table scoped to this sheet.
        for item in sample_records:
            row = int(item["sheet_row"])
            xlrow = row + 1
            mean_ref = "{}{}".format(xl_col_to_name(mean_col), xlrow)
            power_cell = "{}{}".format(xl_col_to_name(power_col), xlrow)
            base_cell = "{}{}".format(xl_col_to_name(base_ok_col), xlrow)
            slope_ok = 'AND({0}="Ready",ISNUMBER({1}),{1}<0,ISNUMBER({2}))'.format(qc_excel_ref, slope_ref, intercept_ref)
            estimate = '10^(({0}-{1})/{2})'.format(mean_ref, intercept_ref, slope_ref)
            tolerance = '1E-9*MAX(ABS({0}),ABS({1}),ABS({2}),1E-300)'.format(min_q_ref, max_q_ref, estimate)
            out_range = 'OR({0}<({1}-{2}),{0}>({3}+{2}),{4}<({5}-{6}),{4}>({7}+{6}))'.format(
                estimate, min_q_ref, tolerance, max_q_ref, mean_ref, min_y_ref,
                '1E-9*MAX(ABS({}),ABS({}),1E-300)'.format(min_y_ref, max_y_ref), max_y_ref,
            )
            out_decision = 'IF({0},{1},"")'.format(out_range, estimate) if allow_extrapolation else 'IF({0},"",{1})'.format(out_range, estimate)
            power_formula = '=IF(OR(NOT(ISNUMBER({0})),NOT({1})),"",IFERROR({2},""))'.format(mean_ref, slope_ok, out_decision)
            power_cache: Any = ""
            if fit_cache.get("status") == "Ready" and _finite(item.get("mean")):
                estimate_value = 10.0 ** ((float(item["mean"]) - fit_cache["intercept"]) / fit_cache["slope"])
                if math.isfinite(estimate_value) and estimate_value > 0:
                    in_range = (
                        fit_cache["min_cq"] <= float(item["mean"]) <= fit_cache["max_cq"]
                        and fit_cache["min_quantity"] * (1 - 1e-9) <= estimate_value <= fit_cache["max_quantity"] * (1 + 1e-9)
                    )
                    if in_range or allow_extrapolation:
                        power_cache = estimate_value
            sheet.write_formula(row, power_col, power_formula, formats["quantity"], power_cache)
            # QC communicates the source sample review decision, fit state and
            # range handling. It is formula-backed for the dynamic conditions.
            base_status = item["base_status"].replace('"', '""')
            if not item["base_eligible"]:
                status_formula = '="{}"'.format(base_status or "Requires review")
            else:
                if allow_extrapolation:
                    dynamic_status = 'IF({0}="","No retained numeric Cq",IF(NOT({1}),"Standard curve requires review",IF({2},"Provisional extrapolation",IF(OR("{3}"="Ready","{3}"="Single numeric replicate"),"Ready","{3}"))))'.format(
                        mean_ref, slope_ok, out_range, base_status
                    )
                else:
                    dynamic_status = 'IF({0}="","No retained numeric Cq",IF(NOT({1}),"Standard curve requires review",IF({2},"Outside calibrated range; final quantity withheld",IF(OR("{3}"="Ready","{3}"="Single numeric replicate"),"Ready","{3}"))))'.format(
                        mean_ref, slope_ok, out_range, base_status
                    )
                status_formula = "=" + dynamic_status
            status_cache = base_status if not item["base_eligible"] else (
                "No retained numeric Cq" if item["mean"] is None else (
                    "Standard curve requires review" if fit_cache.get("status") != "Ready" else (
                        "Provisional extrapolation" if power_cache != "" and (not fit_cache["min_quantity"] <= float(power_cache) <= fit_cache["max_quantity"]) else
                        "Outside calibrated range; final quantity withheld" if power_cache == "" else
                        ("Ready" if base_status in {"Ready", "Single numeric replicate"} else base_status)
                    )
                )
            )
            sheet.write_formula(row, status_col, status_formula, formats["warning"] if status_cache != "Ready" else formats["text"], status_cache)
            sheet.write_formula(row, eligible_col, '=IF(AND(ISNUMBER({0}),{0}>0,{1}=TRUE),"Yes","No")'.format(power_cell, base_cell), formats["text"],
                                "Yes" if item["base_eligible"] and power_cache != "" else "No")

        if sample_records:
            sheet.autofilter(0, 0, sample_last_row, len(visible_headers) - 1)
        # Main table widths adapt to replicate count while keeping identifiers
        # readable and sample metadata compact.
        for col, width in enumerate([22, 16, 16, 14, 24] + [13] * max_reps + [13, 16, 34]):
            sheet.set_column(col, col, width)
        sheet.set_row(0, 30)

        # Add chart and points using the live selected helper columns.
        chart = workbook.add_chart({"type": "scatter", "subtype": "straight_with_markers"})
        chart.set_title({"name": "{} standard curve".format(gene.lower())})
        chart.set_x_axis({"name": "log10 quantity (AU)", "major_gridlines": {"visible": False}, "line": {"color": "#7F8FA6"}})
        chart.set_y_axis({"name": "Cq", "major_gridlines": {"visible": True, "line": {"color": "#D9E2F3"}}, "line": {"color": "#7F8FA6"}})
        chart.set_legend({"none": True})
        chart.set_size({"width": 520, "height": 300})
        chart.show_hidden_data()
        chart.show_blanks_as("span")  # Non-contiguous included levels remain connected without adding points.
        if standard_records:
            cached_chart_points = []
            if fit_cache.get("status") == "Ready":
                cached_chart_points = [
                    (math.log10(float(item["Concentration"])), float(np.mean(item["readings"])))
                    for item in initial_standard_rows
                    if item["Include"] and _finite(item.get("Concentration"))
                    and float(item["Concentration"]) > 0 and item["readings"]
                ]
            chart.add_series({
                "name": "{} selected standards".format(gene),
                "categories": [sheet_name, std_first_row, std_x_col, max(std_last_row, std_first_row), std_x_col],
                "categories_data": [point[0] for point in cached_chart_points],
                "values": [sheet_name, std_first_row, std_y_col, max(std_last_row, std_first_row), std_y_col],
                "values_data": [point[1] for point in cached_chart_points],
                "marker": {"type": "circle", "size": 6, "border": {"color": "#4472C4"}, "fill": {"color": "#4472C4"}},
                "line": {"color": "#4472C4", "width": 1.5},
            })
        sheet.insert_chart(chart_row, chart_col, chart, {"x_offset": 8, "y_offset": 0})

        # Raw Results rows are kept verbatim beneath this assay/plate review.
        raw_first = max(42, fit_summary_row + len(summary_labels) + 3, sample_first_row + len(sample_records) + 3)
        sheet.write_string(raw_first, 0, "Raw Results", formats["section"])
        raw_header_row = raw_first + 1
        source_columns = list(source_frame.columns)
        for col, header in enumerate(source_columns):
            sheet.write_string(raw_header_row, col, str(header), formats["raw_header"])
        raw_group_mask = prepared["Gene"].map(_text).eq(gene) & prepared["Plate"].map(_text).eq(plate)
        # Keep rows lacking plate/gene identity with the first assay block so
        # every source row remains represented once in the workbook.
        if group_index == 0:
            raw_group_mask = raw_group_mask | prepared["Gene"].map(_text).eq("") | prepared["Plate"].map(_text).eq("")
        raw_positions = np.flatnonzero(raw_group_mask.to_numpy())
        raw_rows = source_frame.iloc[raw_positions]
        for offset, (_, raw_row) in enumerate(raw_rows.iterrows(), start=raw_header_row + 1):
            for col, source_col in enumerate(source_columns):
                _write_literal(sheet, offset, col, raw_row.get(source_col))
        for col, header in enumerate(source_columns):
            width = min(38, max(12, len(str(header)) + 2))
            if col >= len(visible_headers):
                sheet.set_column(col, col, width)
        print_last_col = max(std_start_col + len(std_headers) - 1, len(source_columns) - 1)
        print_last_row = raw_header_row + len(raw_rows)
        sheet.print_area(0, 0, print_last_row, print_last_col)
        assay_details[(gene, plate)] = {
            "sheet": sheet_name,
            "sample_records": sample_records,
            "fit_cache": fit_cache,
            "power_col": power_col,
            "eligible_col": eligible_col,
            "status_col": status_col,
            "sample_rows": {item["SampleId"]: int(item["sheet_row"]) for item in sample_records},
        }

    # A single compact normalization table, linked by exact SampleId. It uses
    # only the Target/reference ratio endpoint, never the calibrator fold row.
    normalized = workbook.add_worksheet("normalised")
    normalized.set_default_row(18)
    normalized.freeze_panes(1, 0)
    normalized.set_landscape()
    normalized.set_paper(8)
    normalized.fit_to_pages(1, 0)
    normalized.set_margins(0.25, 0.25, 0.35, 0.35)
    normalized.center_horizontally()
    norm_headers = ["SampleId", "Sample", "Group", "Target", "Target/reference ratio", "QC / inclusion"]
    for col, header in enumerate(norm_headers):
        normalized.write_string(0, col, header, formats["header"])
    norm_rows: List[Dict[str, Any]] = []
    if not ratio_results.empty:
        for _, result in ratio_results.iterrows():
            if _text(result.get("Endpoint")) != "Target/reference ratio":
                continue
            norm_rows.append(result.to_dict())
    for index, item in enumerate(norm_rows, start=1):
        sample_id = _text(item.get("SampleId"))
        gene = _text(item.get("Gene"))
        normalized.write_string(index, 0, sample_id, formats["text"])
        normalized.write(index, 1, item.get("Label"), formats["text"])
        normalized.write(index, 2, item.get("Group"), formats["text"])
        normalized.write_string(index, 3, gene, formats["text"])
        # Exact sample/gene matches are required. If any assay maps to multiple
        # physical plates, expose a blank ratio and a specific review status.
        target_matches = [(g, p, detail) for (g, p), detail in assay_details.items()
                          if g == gene and sample_id in detail["sample_rows"]]
        ref_matches: List[Tuple[str, str, Dict[str, Any]]] = []
        for reference in references:
            matches = [(g, p, detail) for (g, p), detail in assay_details.items()
                       if g == reference and sample_id in detail["sample_rows"]]
            ref_matches.extend(matches)
        cached = item.get("Value") if _truthy(item.get("Eligible")) else ""
        formula: Optional[str] = None
        status_formula: Optional[str] = None
        status_cache = _text(item.get("Status")) or "Normalized standard-curve result unavailable"
        if len(target_matches) == 1 and len(ref_matches) == len(references) and references and all(
            sum(1 for match in ref_matches if match[0] == reference) == 1 for reference in references
        ):
            _g, _p, target_detail = target_matches[0]
            target_row = target_detail["sample_rows"][sample_id]
            target_power = "{}!{}{}".format(
                _quote_sheet(target_detail["sheet"]), xl_col_to_name(target_detail["power_col"]), target_row + 1
            )
            target_ok = "{}!{}{}".format(
                _quote_sheet(target_detail["sheet"]), xl_col_to_name(target_detail["eligible_col"]), target_row + 1
            )
            ref_powers = []
            ref_oks = []
            ref_statuses = []
            for reference in references:
                match = next((match for match in ref_matches if match[0] == reference), None)
                if match is None:
                    continue
                _rg, _rp, ref_detail = match
                ref_row = ref_detail["sample_rows"][sample_id]
                ref_powers.append("{}!{}{}".format(_quote_sheet(ref_detail["sheet"]), xl_col_to_name(ref_detail["power_col"]), ref_row + 1))
                ref_oks.append("{}!{}{}".format(_quote_sheet(ref_detail["sheet"]), xl_col_to_name(ref_detail["eligible_col"]), ref_row + 1))
                ref_statuses.append("{}!{}{}".format(_quote_sheet(ref_detail["sheet"]), xl_col_to_name(ref_detail["status_col"]), ref_row + 1))
            conditions = ["{}=\"Yes\"".format(target_ok), "ISNUMBER({})".format(target_power), "{}>0".format(target_power)]
            for power, ok in zip(ref_powers, ref_oks):
                conditions.extend(["{}=\"Yes\"".format(ok), "ISNUMBER({})".format(power), "{}>0".format(power)])
            formula = "=IF(AND({}),{}/GEOMEAN({}),\"\")".format(",".join(conditions), target_power, ",".join(ref_powers))
            # The status follows the same live eligibility gate and recognizes
            # an allowed extrapolation on either target or reference assay.
            provisional_checks = ["ISNUMBER(SEARCH(\"Provisional extrapolation\",{}))".format(
                "{}!{}{}".format(_quote_sheet(target_detail["sheet"]), xl_col_to_name(target_detail["status_col"]), target_row + 1)
            )] + ["ISNUMBER(SEARCH(\"Provisional extrapolation\",{}))".format(value) for value in ref_statuses]
            status_formula = '=IF(AND({0}),IF(OR({1}),"Provisional standard-curve ratio","Ready"),"Target/reference quantity unavailable or ineligible")'.format(
                ",".join(conditions), ",".join(provisional_checks) or "FALSE"
            )
            if _truthy(item.get("Eligible")) and _finite(item.get("Value")):
                cached = float(item["Value"])
                if "Provisional" in status_cache:
                    status_cache = "Provisional standard-curve ratio"
                else:
                    status_cache = "Ready"
            else:
                cached = ""
                status_cache = "Target/reference quantity unavailable or ineligible"
        elif len(target_matches) > 1 or any(sum(1 for match in ref_matches if match[0] == reference) > 1 for reference in references):
            status_cache = "Ambiguous plate mapping; review required"
        elif not references:
            status_cache = "No reference genes selected"
        else:
            status_cache = "Target or reference measurement is missing"
        if formula:
            normalized.write_formula(index, 4, formula, formats["ratio"], cached)
            normalized.write_formula(index, 5, status_formula or '="{}"'.format(status_cache.replace('"', '""')),
                                     formats["warning"] if status_cache != "Ready" else formats["text"], status_cache)
        else:
            normalized.write_blank(index, 4, None, formats["ratio"])
            normalized.write_string(index, 5, status_cache, formats["warning"] if status_cache != "Ready" else formats["text"])
    normalized.set_column(0, 0, 24)
    normalized.set_column(1, 1, 22)
    normalized.set_column(2, 3, 18)
    normalized.set_column(4, 4, 24)
    normalized.set_column(5, 5, 44)
    if norm_rows:
        normalized.autofilter(0, 0, len(norm_rows), len(norm_headers) - 1)
    normalized.print_area(0, 0, max(1, len(norm_rows)), len(norm_headers) - 1)

    workbook.close()
    return output.getvalue()


__all__ = ["cyc1_review_workbook"]
