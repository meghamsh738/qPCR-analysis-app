"""Synthetic-only contract checks for the editable qPCR review export."""

from io import BytesIO
from pathlib import Path
import shutil
import subprocess
import tempfile
from zipfile import ZipFile

import numpy as np
import pandas as pd
from openpyxl import load_workbook

from qpcr_review_export import cyc1_review_workbook
from qpcr_workflow import analyze, make_standard_map, prepare_wells


def _fixture(plates=("P1",), allow_extrapolation=False):
    rows = []
    concentrations = [("S1", 1000.0), ("S2", 100.0), ("S3", 10.0), ("S4", 1.0)]
    slope, intercept = -3.32, 35.0
    for plate in plates:
        for gene, offset in (("CYC1", 0.0), ("Ref", 0.2)):
            for label, quantity in concentrations:
                cq = intercept + offset + slope * np.log10(quantity)
                if label == "S2":
                    cq += 0.4  # excluded initially; inclusion should change the live fit
                for replicate, shift in ((1, -0.05), (2, 0.05)):
                    rows.append({
                        "Instrument well": "{}-{}-{}-{}".format(plate, gene, label, replicate),
                        "Plate": plate, "Well": "{}{}".format(label, replicate),
                        "Gene": gene, "Type": "Standard", "Label": label,
                        "Cq": cq + shift, "Original comment": "synthetic",
                    })
            sample_rows = [("A", 30.0), ("Out", 40.0)] if gene == "CYC1" else [("A", 30.0), ("Out", 30.0)]
            for sample_id, cq in sample_rows:
                rows.append({
                    "Instrument well": "{}-{}-{}-1".format(plate, gene, sample_id),
                    "Plate": plate, "Well": "{}1".format(sample_id), "Gene": gene,
                    "Type": "Sample", "Label": "Mouse {}".format(sample_id),
                    "Cq": cq, "SampleId": sample_id, "Group": "Control",
                    "PairId": "pair-{}".format(sample_id), "Original comment": "synthetic",
                })
                if gene == "CYC1" and sample_id == "A":
                    rows.append({
                        "Instrument well": "{}-{}-A-2".format(plate, gene),
                        "Plate": plate, "Well": "A2", "Gene": gene, "Type": "Sample",
                        "Label": "Mouse A", "Cq": 33.0, "SampleId": "A", "Group": "Control",
                        "PairId": "pair-A", "Keep": False, "Reason": "reviewer excluded",
                        "Original comment": "excluded source row",
                    })
    wells = prepare_wells(pd.DataFrame(rows))
    standards = make_standard_map(wells, top=1000.0, factor=10.0)
    standards["Include"] = ~standards["Label"].eq("S2")
    standards.loc[standards["Label"].eq("S2"), "Reason"] = "synthetic review exclusion"
    results = analyze(
        wells,
        standards,
        references=["Ref"],
        calibrators=["A"],
        allow_extrapolation=allow_extrapolation,
        efficiency_confirmed=True,
    )
    settings = {"references": ["Ref"], "allow_extrapolation": allow_extrapolation}
    return wells, standards, results, settings


def _bytes_workbook(wells, standards, results, settings):
    content = cyc1_review_workbook(wells, standards, results, settings)
    assert isinstance(content, bytes) and content[:2] == b"PK"
    return content


def test_export_has_live_fit_mean_quantity_and_id_matched_ratio_with_formula_caches():
    wells, standards, results, settings = _fixture()
    content = _bytes_workbook(wells, standards, results, settings)
    formulas = load_workbook(BytesIO(content), data_only=False)
    values = load_workbook(BytesIO(content), data_only=True)

    assert "cyc1" in formulas.sheetnames
    assert "ref" in formulas.sheetnames
    assert "normalised" in formulas.sheetnames
    assay = formulas["cyc1"]
    assay_values = values["cyc1"]
    assert assay["A1"].value == "Sample"
    assert assay["C1"].value == "Pair ID"
    assert not assay.sheet_format.zeroHeight  # Blank review rows must leave room for the live chart.
    # Main assay review starts at A1 and uses formula cells for the mean/AU.
    assert assay["H2"].data_type == "f"  # Cq mean
    assert assay["I2"].data_type == "f"  # CYC1 Power (AU)
    assert np.isclose(assay_values["H2"].value, 30.0)
    assert assay_values["I2"].value > 0
    assert assay["T21"].value == "Yes"  # S1 selected
    assert assay["T22"].value == "No"  # S2 excluded
    assert assay["T23"].value == "Yes"  # S3 selected: non-contiguous set
    assert assay["T24"].value == "Yes"  # S4 selected
    helper_columns = {cell.value: cell.column for cell in assay[20] if cell.value}
    selected_log_col = helper_columns["Selected log10 quantity"]
    selected_y_col = helper_columns["Selected Cq mean"]
    selected_quantity_col = helper_columns["Selected quantity"]
    assert assay.cell(21, selected_log_col).data_type == "f"
    assert assay.cell(21, selected_quantity_col).data_type == "f"
    summary_rows = {cell.value: cell.row for row in assay.iter_rows() for cell in row if cell.value in {
        "Fitted equation", "Slope", "Intercept", "R-squared", "Efficiency", "Blank wells", "Assay QC"
    }}
    qc_row = summary_rows["Assay QC"]
    assert assay.cell(qc_row, 15).data_type == "f"
    assert assay_values.cell(qc_row, 15).value == "Ready"

    normalized = formulas["normalised"]
    normalized_values = values["normalised"]
    a_row = next(row for row in range(2, normalized.max_row + 1) if normalized.cell(row, 1).value == "A")
    assert normalized.cell(a_row, 5).data_type == "f"
    assert normalized_values.cell(a_row, 5).value > 0
    assert normalized_values.cell(a_row, 6).value == "Ready"

    # Native chart uses hidden selected-point helper ranges rather than the
    # visible standard panel's contiguous rows (which include excluded S2).
    assert len(formulas["cyc1"]._charts) == 1
    chart_series = formulas["cyc1"]._charts[0].ser[0]
    chart_formula = chart_series.yVal.numRef.f
    assert "cyc1" in chart_formula.casefold()
    y_letter = assay.cell(20, selected_y_col).column_letter
    x_letter = assay.cell(20, selected_log_col).column_letter
    assert "{}21:{}24".format(y_letter, y_letter) in chart_formula.replace("$", "")
    x_formula = chart_series.xVal.numRef.f.replace("$", "")
    assert "{}21:{}24".format(x_letter, x_letter) in x_formula
    with ZipFile(BytesIO(content)) as archive:
        chart_xml = archive.read("xl/charts/chart1.xml").decode("utf-8")
    assert "<c:ptCount val=\"3\"/>" in chart_xml
    assert "<c:plotVisOnly" not in chart_xml  # chart intentionally reads hidden live helpers
    assert "SLOPE" in assay.cell(summary_rows["Slope"], 15).value
    assert "INTERCEPT" in assay.cell(summary_rows["Intercept"], 15).value
    assert "RSQ" in assay.cell(summary_rows["R-squared"], 15).value


def test_excluded_replicate_is_blank_in_main_table_but_original_row_is_retained_raw():
    wells, standards, results, settings = _fixture()
    content = _bytes_workbook(wells, standards, results, settings)
    formulas = load_workbook(BytesIO(content), data_only=False)
    values = load_workbook(BytesIO(content), data_only=True)
    sheet = formulas["cyc1"]
    value_sheet = values["cyc1"]
    # Excluded numeric well is omitted from the mean inputs; the remaining
    # detected replicate is shown once, with no fabricated second value.
    assert value_sheet["F2"].value == 30.0
    assert value_sheet["G2"].value is None
    raw_header = next(row for row in range(1, sheet.max_row + 1) if sheet.cell(row, 1).value == "Instrument well")
    raw_values = [sheet.cell(row, 1).value for row in range(raw_header + 1, sheet.max_row + 1)]
    assert "P1-CYC1-A-2" in raw_values
    assert "Original comment" in [sheet.cell(raw_header, col).value for col in range(1, sheet.max_column + 1)]


def test_out_of_range_is_withheld_by_default_and_provisional_when_explicitly_allowed():
    wells, standards, results, settings = _fixture(allow_extrapolation=False)
    content = _bytes_workbook(wells, standards, results, settings)
    formula_wb = load_workbook(BytesIO(content), data_only=False)
    value_wb = load_workbook(BytesIO(content), data_only=True)
    assay, assay_values = formula_wb["cyc1"], value_wb["cyc1"]
    out_row = next(row for row in range(2, assay.max_row + 1) if assay.cell(row, 5).value == "Out")
    assert assay.cell(out_row, 9).data_type == "f"
    assert assay_values.cell(out_row, 9).value is None
    norm = value_wb["normalised"]
    out_norm_row = next(row for row in range(2, norm.max_row + 1) if norm.cell(row, 1).value == "Out")
    assert norm.cell(out_norm_row, 5).value is None

    wells2, standards2, results2, settings2 = _fixture(allow_extrapolation=True)
    content2 = _bytes_workbook(wells2, standards2, results2, settings2)
    formula_wb2 = load_workbook(BytesIO(content2), data_only=False)
    value_wb2 = load_workbook(BytesIO(content2), data_only=True)
    assay2, assay_values2 = formula_wb2["cyc1"], value_wb2["cyc1"]
    out_row2 = next(row for row in range(2, assay2.max_row + 1) if assay2.cell(row, 5).value == "Out")
    assert assay_values2.cell(out_row2, 9).value > 0
    assert assay_values2.cell(out_row2, 10).value == "Provisional extrapolation"
    norm2 = value_wb2["normalised"]
    out_norm_row2 = next(row for row in range(2, norm2.max_row + 1) if norm2.cell(row, 1).value == "Out")
    assert norm2.cell(out_norm_row2, 5).value > 0
    assert norm2.cell(out_norm_row2, 6).value == "Provisional standard-curve ratio"


def test_multiple_plates_get_separate_sheets_and_ambiguous_ids_do_not_normalise():
    wells, standards, results, settings = _fixture(plates=("P1", "P2"))
    content = _bytes_workbook(wells, standards, results, settings)
    formula_wb = load_workbook(BytesIO(content), data_only=False)
    value_wb = load_workbook(BytesIO(content), data_only=True)
    assert {"cyc1 - p1", "cyc1 - p2", "ref - p1", "ref - p2"}.issubset(
        {name.casefold() for name in formula_wb.sheetnames}
    )
    assert all(len(formula_wb[name]._charts) == 1 for name in ("cyc1 - p1", "cyc1 - p2"))
    normalized = value_wb["normalised"]
    a_row = next(row for row in range(2, normalized.max_row + 1) if normalized.cell(row, 1).value == "A")
    assert normalized.cell(a_row, 5).value is None
    assert normalized.cell(a_row, 6).value == "Ambiguous plate mapping; review required"


def test_sheet_names_are_lowercase_and_source_strings_are_not_exported_as_formulas():
    wells, standards, results, settings = _fixture()
    wells.loc[wells.index[0], "Original comment"] = "=1+1"
    content = _bytes_workbook(wells, standards, results, settings)
    formula_wb = load_workbook(BytesIO(content), data_only=False)
    with ZipFile(BytesIO(content)) as archive:
        workbook_xml = archive.read("xl/workbook.xml").decode("utf-8")
    assert all(name == name.lower() for name in formula_wb.sheetnames)
    assert formula_wb.custom_doc_props["EasyLabWorkbookKind"].value == "qpcr-review"
    assert "<f>1+1</f>" not in workbook_xml
    raw = formula_wb["cyc1"]
    raw_header = next(row for row in range(1, raw.max_row + 1) if raw.cell(row, 1).value == "Instrument well")
    comment_col = next(col for col in range(1, raw.max_column + 1) if raw.cell(raw_header, col).value == "Original comment")
    assert raw.cell(raw_header + 1, comment_col).value == "=1+1"
    assert raw.cell(raw_header + 1, comment_col).data_type == "s"


def test_export_docstring_states_recalculation_and_reimport_limitations():
    doc = cyc1_review_workbook.__doc__ or ""
    assert "re-import" in doc
    assert "does not evaluate them" in doc
    assert "plate" in doc.lower()


def test_live_standard_inclusion_recalculates_in_headless_spreadsheet_engine():
    soffice = shutil.which("soffice")
    if not soffice:
        import pytest
        pytest.skip("No headless spreadsheet recalculation engine is installed")
    wells, standards, results, settings = _fixture()
    content = _bytes_workbook(wells, standards, results, settings)
    with tempfile.TemporaryDirectory(prefix="qpcr-review-export-") as temp:
        root = Path(temp)
        source = root / "source.xlsx"
        output_dir = root / "recalculated"
        profile = root / "lo-profile"
        output_dir.mkdir()
        source.write_bytes(content)
        editable = load_workbook(source)
        sheet = editable["cyc1"]
        # S2 is the non-contiguous excluded level; switching it on must update
        # the live point set and the slope/quantity formulas after recalculation.
        sheet["T22"] = "Yes"
        editable.save(source)
        run = subprocess.run(
            [
                soffice,
                "--headless",
                "-env:UserInstallation={}".format(profile.as_uri()),
                "--convert-to", "xlsx",
                "--outdir", str(output_dir), str(source),
            ],
            check=False,
            capture_output=True,
            text=True,
            timeout=45,
        )
        recalculated = output_dir / source.name
        assert run.returncode == 0, run.stderr or run.stdout
        assert recalculated.exists(), run.stderr or run.stdout
        values = load_workbook(recalculated, data_only=True)
        assay = values["cyc1"]
        # The fit summary starts beneath the standards table. Inclusion of the
        # shifted S2 point changes the fit and leaves the curve eligible.
        slope_before = results["curves"].loc[
            results["curves"]["Gene"].eq("CYC1"), "Slope"
        ].iloc[0]
        labels = {assay.cell(row, 14).value: row for row in range(1, assay.max_row + 1)}
        slope_after = assay.cell(labels["Slope"], 15).value
        assert np.isfinite(slope_after)
        assert not np.isclose(slope_after, slope_before)
        assert assay.cell(labels["Assay QC"], 15).value == "Ready"

        # The visible Reason cell is live QC input when a level is deselected.
        editable = load_workbook(recalculated)
        editable["cyc1"]["T22"] = "No"
        editable["cyc1"]["U22"] = ""
        source2 = root / "source2.xlsx"
        editable.save(source2)
        output_dir2 = root / "recalculated2"
        output_dir2.mkdir()
        run2 = subprocess.run(
            [soffice, "--headless", "-env:UserInstallation={}".format(profile.as_uri()),
             "--convert-to", "xlsx", "--outdir", str(output_dir2), str(source2)],
            check=False, capture_output=True, text=True, timeout=45,
        )
        recalculated2 = output_dir2 / source2.name
        assert run2.returncode == 0, run2.stderr or run2.stdout
        assert recalculated2.exists(), run2.stderr or run2.stdout
        values2 = load_workbook(recalculated2, data_only=True)
        assay2 = values2["cyc1"]
        labels2 = {assay2.cell(row, 14).value: row for row in range(1, assay2.max_row + 1)}
        assert assay2.cell(labels2["Assay QC"], 15).value == "Selected standard mapping or reading requires review"

        editable2 = load_workbook(recalculated2)
        editable2["cyc1"]["U22"] = "synthetic review reason"
        source3 = root / "source3.xlsx"
        editable2.save(source3)
        output_dir3 = root / "recalculated3"
        output_dir3.mkdir()
        run3 = subprocess.run(
            [soffice, "--headless", "-env:UserInstallation={}".format(profile.as_uri()),
             "--convert-to", "xlsx", "--outdir", str(output_dir3), str(source3)],
            check=False, capture_output=True, text=True, timeout=45,
        )
        recalculated3 = output_dir3 / source3.name
        assert run3.returncode == 0, run3.stderr or run3.stdout
        assert recalculated3.exists(), run3.stderr or run3.stdout
        values3 = load_workbook(recalculated3, data_only=True)
        assay3 = values3["cyc1"]
        labels3 = {assay3.cell(row, 14).value: row for row in range(1, assay3.max_row + 1)}
        assert assay3.cell(labels3["Assay QC"], 15).value == "Ready"


def _recalculate_fixture(wells, standards, allow_extrapolation=False):
    prepared = prepare_wells(wells)
    results = analyze(
        prepared, standards, references=["Ref"], calibrators=["A"],
        allow_extrapolation=allow_extrapolation, efficiency_confirmed=True,
    )
    return prepared, results, {"references": ["Ref"], "allow_extrapolation": allow_extrapolation}


def test_formula_fit_cache_matches_pipeline_for_undetermined_invalid_and_empty_standards():
    wells, standards, _results, settings = _fixture()
    cases = [
        ("one numeric plus undetermined", "Undetermined", "Ready"),
        ("retained invalid text", "not a Cq", "Selected standard mapping or reading requires review"),
        ("no retained numeric Cq", "Undetermined", "Selected standard mapping or reading requires review"),
    ]
    for name, replacement, expected in cases:
        altered = wells.copy(deep=True)
        altered["Cq"] = altered["Cq"].astype(object)
        altered["OriginalCq"] = altered["OriginalCq"].astype(object)
        target = altered.index[(altered["Gene"] == "CYC1") & (altered["Label"] == "S1")].tolist()
        if name == "no retained numeric Cq":
            changed_rows = target
        else:
            changed_rows = [target[-1]]
        for row in changed_rows:
            altered.at[row, "Cq"] = replacement
            altered.at[row, "OriginalCq"] = replacement
        prepared, result, case_settings = _recalculate_fixture(altered, standards.copy(deep=True), settings["allow_extrapolation"])
        curve = result["curves"].loc[result["curves"]["Gene"].eq("CYC1")].iloc[0]
        content = _bytes_workbook(prepared, standards, result, case_settings)
        values = load_workbook(BytesIO(content), data_only=True)
        assay_sheet = values["cyc1"]
        # Initial export QC is taken from the pipeline result and the formula
        # cache must agree with it for this selected-standard state.
        qc_row = next(row for row in range(1, assay_sheet.max_row + 1) if assay_sheet.cell(row, 14).value == "Assay QC")
        qc_value = assay_sheet.cell(qc_row, 15).value
        assert qc_value == curve["Status"], name
        assert qc_value == expected, name


def test_unselected_invalid_concentration_does_not_block_fit_and_overflow_is_exportable():
    wells, standards, _results, settings = _fixture()
    altered_standards = standards.copy(deep=True)
    unselected = altered_standards["Label"].eq("S2")
    altered_standards.loc[unselected, "Concentration"] = np.nan
    altered_standards.loc[unselected, "Include"] = False
    altered_standards.loc[unselected, "Reason"] = "synthetic review exclusion"
    prepared, results, case_settings = _recalculate_fixture(wells, altered_standards, False)
    curve = results["curves"].loc[results["curves"]["Gene"].eq("CYC1")].iloc[0]
    assert bool(curve["Eligible"])
    content = _bytes_workbook(prepared, altered_standards, results, case_settings)
    values = load_workbook(BytesIO(content), data_only=True)
    assay = values["cyc1"]
    qc_row = next(row for row in range(1, assay.max_row + 1) if assay.cell(row, 14).value == "Assay QC")
    assert assay.cell(qc_row, 15).value == "Ready"

    shallow = wells.copy(deep=True)
    concentration_by_label = standards.set_index("Label")["Concentration"].to_dict()
    for row in shallow.index[(shallow["Gene"] == "CYC1") & (shallow["Type"] == "Standard")]:
        quantity = float(concentration_by_label[shallow.at[row, "Label"]])
        cq = 20.0 - 0.001 * np.log10(quantity)
        shallow.at[row, "Cq"] = cq
        shallow.at[row, "OriginalCq"] = cq
    prepared_shallow, shallow_results, shallow_settings = _recalculate_fixture(shallow, standards.copy(deep=True), False)
    shallow_curve = shallow_results["curves"].loc[shallow_results["curves"]["Gene"].eq("CYC1")].iloc[0]
    assert shallow_curve["Status"] == "Invalid efficiency calculation"
    shallow_content = _bytes_workbook(prepared_shallow, standards, shallow_results, shallow_settings)
    shallow_values = load_workbook(BytesIO(shallow_content), data_only=True)
    shallow_assay = shallow_values["cyc1"]
    shallow_qc_row = next(row for row in range(1, shallow_assay.max_row + 1) if shallow_assay.cell(row, 14).value == "Assay QC")
    assert shallow_assay.cell(shallow_qc_row, 15).value == "Invalid efficiency calculation"
