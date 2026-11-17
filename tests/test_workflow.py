import numpy as np
import pandas as pd

from qpcr_workflow import (
    analyze,
    fit_curves,
    make_standard_map,
    prepare_wells,
    quantify,
    summarise_samples,
)


def _base_rows():
    rows = []
    slope = -3.32192809489
    intercepts = {"Target": 24.0, "RefA": 24.0, "RefB": 24.0}
    for gene, intercept in intercepts.items():
        for label, concentration in [("Std1", 16.0), ("Std2", 4.0), ("Std3", 1.0)]:
            cq = intercept + slope * np.log10(concentration)
            for replicate, shift in [(1, -0.05), (2, 0.05)]:
                rows.append({
                    "Plate": "P1", "Well": "S{}-{}-{}".format(gene, label, replicate),
                    "Gene": gene, "Type": "Standard", "Label": label,
                    "Cq": cq + shift, "Replicate": replicate,
                })
    # Two biological samples, deliberately with unequal technical replicate
    # counts so tests can catch well-weighting errors.
    samples = {
        "A": {"Group": "Control", "PairId": "pair-1", "Target": [22.0, 22.0, 22.0, 22.0], "RefA": [23.0, 23.0], "RefB": [24.0, 24.0]},
        "B": {"Group": "Treatment", "PairId": "pair-2", "Target": [20.0], "RefA": [22.0, 22.0, 22.0], "RefB": [22.0]},
    }
    for sample, values in samples.items():
        for gene in ("Target", "RefA", "RefB"):
            for replicate, cq in enumerate(values[gene], start=1):
                rows.append({
                    "Plate": "P1", "Well": "{}-{}-{}".format(sample, gene, replicate),
                    "Gene": gene, "Type": "Sample", "Label": sample,
                    "Cq": cq, "SampleId": sample, "Group": values["Group"],
                    "PairId": values["PairId"],
                })
    return pd.DataFrame(rows)


def test_prepare_wells_preserves_rows_and_distinguishes_undetermined_invalid_and_missing():
    raw = pd.DataFrame([
        {"Plate": "P1", "Well": "A1", "Gene": "G", "Type": "Sample", "Label": "mouse1", "Cq": "22.4"},
        {"Plate": "P1", "Well": "A2", "Gene": "G", "Type": "Sample", "Label": "mouse1", "Cq": "Undetermined"},
        {"Plate": "P1", "Well": "A3", "Gene": "G", "Type": "Sample", "Label": "mouse2", "Cq": "not-a-cq"},
        {"Plate": "P1", "Well": "A4", "Gene": "G", "Type": "Sample", "Label": "mouse2", "Cq": None},
    ])
    wells = prepare_wells(raw)
    assert len(wells) == len(raw)
    assert wells["Well"].tolist() == raw["Well"].tolist()
    assert wells["Cq"].iloc[0] == 22.4
    assert wells["CqStatus"].tolist() == ["Detected", "Undetermined", "Invalid", "Missing"]
    assert wells["OriginalCq"].iloc[2] == "not-a-cq"
    assert "numeric" in wells["Issue"].iloc[2].lower()


def test_automatic_ids_are_scoped_to_plate_and_unknown_type_needs_mapping():
    raw = pd.DataFrame([
        {"Plate": "P1", "Well": "A1", "Gene": "G", "Type": "", "Label": "101", "Cq": 20},
        {"Plate": "P2", "Well": "A1", "Gene": "G", "Type": "", "Label": "101", "Cq": 20},
        {"Plate": "P1", "Well": "A2", "Gene": "G", "Type": "mystery", "Label": "101", "Cq": 21},
    ])
    wells = prepare_wells(raw)
    assert wells.loc[0, "SampleId"] != wells.loc[1, "SampleId"]
    assert wells.loc[0, "SampleId"] == 'auto:["P1","101"]'
    assert wells[["Tissue", "Batch"]].isna().all().all()
    assert wells.loc[2, "TypeStatus"] == "Unmapped type"
    assert "map it explicitly" in wells.loc[2, "Issue"]
    summary = summarise_samples(wells)
    unresolved = summary.loc[summary["TypeStatus"].eq("Unmapped type")] if "TypeStatus" in summary else summary.loc[summary["Status"].str.contains("Unrecognized Type")]
    assert len(unresolved) == 1
    assert not bool(unresolved.iloc[0]["Eligible"])


def test_automatic_ids_separate_same_plate_label_by_tissue_and_batch():
    raw = pd.DataFrame([
        {"Plate": "P1", "Well": "A1", "Gene": "G", "Type": "Sample", "Label": "animal-1", "Cq": 20,
         "Tissue": "hippocampus", "Batch": "run-1"},
        {"Plate": "P1", "Well": "A2", "Gene": "G", "Type": "Sample", "Label": "animal-1", "Cq": 21,
         "Tissue": "cortex", "Batch": "run-1"},
        {"Plate": "P1", "Well": "A3", "Gene": "G", "Type": "Sample", "Label": "animal-1", "Cq": 22,
         "Tissue": "hippocampus", "Batch": "run-2"},
    ])
    wells = prepare_wells(raw)
    assert wells["SampleId"].nunique() == 3
    assert wells["SampleId"].str.startswith("auto:").all()
    assert summarise_samples(wells)["Eligible"].all()


def test_explicit_sample_id_links_assays_across_plates_with_matching_context():
    raw = pd.DataFrame([
        {"Plate": "P1", "Well": "A1", "Gene": "Target", "Type": "Sample", "Label": "animal-7", "Cq": 20,
         "SampleId": "animal-7", "Group": "Control", "Tissue": "hippocampus", "Batch": "run-3"},
        {"Plate": "P2", "Well": "B1", "Gene": "Ref", "Type": "Sample", "Label": "animal-7", "Cq": 22,
         "SampleId": "animal-7", "Group": "Control", "Tissue": "hippocampus", "Batch": "run-3"},
    ])
    wells = prepare_wells(raw)
    summary = summarise_samples(wells)
    assert summary["SampleId"].eq("animal-7").all()
    assert summary["Gene"].tolist() == ["Target", "Ref"]
    assert summary["Eligible"].all()
    assert summary["Tissue"].eq("hippocampus").all()
    assert summary["Batch"].eq("run-3").all()


def test_explicit_sample_id_with_conflicting_tissue_or_batch_is_blocked_across_genes():
    raw = pd.DataFrame([
        {"Plate": "P1", "Well": "A1", "Gene": "Target", "Type": "Sample", "Label": "animal-7", "Cq": 20,
         "SampleId": "animal-7", "Group": "Control", "Tissue": "hippocampus", "Batch": "run-1"},
        {"Plate": "P2", "Well": "B1", "Gene": "RefA", "Type": "Sample", "Label": "animal-7", "Cq": 22,
         "SampleId": "animal-7", "Group": "Control", "Tissue": "cortex", "Batch": "run-1"},
        {"Plate": "P2", "Well": "B2", "Gene": "RefB", "Type": "Sample", "Label": "animal-7", "Cq": 23,
         "SampleId": "animal-7", "Group": "Control", "Tissue": "hippocampus", "Batch": "run-2"},
    ])
    wells = prepare_wells(raw)
    assert wells["Issue"].str.contains("Conflicting Tissue metadata for SampleId").all()
    assert wells["Issue"].str.contains("Conflicting Batch metadata for SampleId").all()
    summary = summarise_samples(wells)
    assert not summary["Eligible"].any()
    assert summary["Status"].eq("Conflicting sample metadata").all()
    assert summary["Reason"].str.contains("Conflicting Tissue metadata for SampleId").all()
    assert summary["Reason"].str.contains("Conflicting Batch metadata for SampleId").all()


def test_summaries_average_technical_wells_once_and_leave_disagreement_as_review_flag():
    raw = pd.DataFrame([
        {"Plate": "P1", "Well": "A1", "Gene": "G", "Type": "Sample", "Label": "mouse1", "Cq": 20.0},
        {"Plate": "P1", "Well": "A2", "Gene": "G", "Type": "Sample", "Label": "mouse1", "Cq": 21.0},
        {"Plate": "P1", "Well": "A3", "Gene": "G", "Type": "Sample", "Label": "mouse2", "Cq": 23.0},
    ])
    wells = prepare_wells(raw)
    summary = summarise_samples(wells, threshold=0.75)
    assert len(summary) == 2
    first = summary.loc[summary["Label"].eq("mouse1")].iloc[0]
    assert first["MeanCq"] == 20.5
    assert first["NumericCount"] == 2
    assert first["CqDelta"] == 1.0
    assert bool(first["Disagreement"])
    assert bool(first["Eligible"]), "the reviewed warning does not silently exclude either retained replicate"
    assert summary["SampleId"].nunique() == 2


def test_duplicate_well_identity_and_missing_exclusion_reason_are_blocking_and_visible():
    raw = pd.DataFrame([
        {"Plate": "P1", "Well": "A1", "Gene": "G", "Type": "Sample", "Label": "mouse1", "Cq": 20.0, "Keep": True},
        {"Plate": "P1", "Well": "A1", "Gene": "G", "Type": "Sample", "Label": "mouse1", "Cq": 20.2, "Keep": True},
        {"Plate": "P1", "Well": "A2", "Gene": "G", "Type": "Sample", "Label": "mouse2", "Cq": 21.0, "Keep": False, "Reason": ""},
    ])
    wells = prepare_wells(raw)
    assert wells["IdentityStatus"].iloc[:2].eq("Duplicate well identity").all()
    assert "Exclusion reason is missing" in wells.loc[2, "Issue"]
    summary = summarise_samples(wells)
    mouse1 = summary.loc[summary["Label"].eq("mouse1")].iloc[0]
    mouse2 = summary.loc[summary["Label"].eq("mouse2")].iloc[0]
    assert not bool(mouse1["Eligible"])
    assert not bool(mouse2["Eligible"])


def test_cross_plate_same_assay_is_not_averaged_into_one_sample():
    raw = pd.DataFrame([
        {"Plate": plate, "Well": "A1", "Gene": "G", "Type": "Sample", "Label": "same", "SampleId": "bio1", "Cq": cq}
        for plate, cq in [("P1", 20.0), ("P2", 22.0)]
    ])
    summary = summarise_samples(prepare_wells(raw))
    assert len(summary) == 1
    row = summary.iloc[0]
    assert row["PlateCount"] == 2
    assert pd.isna(row["MeanCq"])
    assert row["PlateTrace"] == "P1 | P2"
    assert not bool(row["Eligible"])


def test_default_standard_map_uses_observed_order_and_editable_serial_values():
    raw = pd.DataFrame([
        {"Plate": "P1", "Well": "A1", "Gene": "G", "Type": "Standard", "Label": label, "Cq": 20}
        for label in ["Top", "Middle", "Bottom"]
    ])
    standards = make_standard_map(prepare_wells(raw), top=16384, factor=4)
    assert standards["Label"].tolist() == ["Top", "Middle", "Bottom"]
    assert standards["Concentration"].tolist() == [16384.0, 4096.0, 1024.0]
    assert standards["Include"].tolist() == [True, True, True]
    assert standards["Reason"].eq("").all()


def test_fit_curves_uses_plate_gene_and_level_means_and_reports_single_standard_reading():
    raw = _base_rows()
    raw["Cq"] = raw["Cq"].astype(object)
    # Make one standard level have only one numeric well. It remains visible
    # and usable as explicitly included, with a warning on the points table.
    standard_rows = raw.index[(raw.Gene == "Target") & (raw.Label == "Std1")].tolist()
    raw.loc[standard_rows[1], "Cq"] = "Undetermined"
    wells = prepare_wells(raw)
    standards = make_standard_map(wells, top=16, factor=4)
    curves, points = fit_curves(wells, standards)
    curve = curves.loc[(curves.Plate == "P1") & (curves.Gene == "Target")].iloc[0]
    assert bool(curve["Eligible"])
    point = points.loc[(points.Gene == "Target") & (points.Label == "Std1")].iloc[0]
    assert point["NumericCount"] == 1
    assert "Only one" in point["Issue"]
    selected_points = points.loc[(points.Gene == "Target") & points.Include]
    expected_slope, expected_intercept = np.polyfit(
        selected_points["Log10Concentration"], selected_points["MeanCq"], 1
    )
    assert np.isclose(curve["Slope"], expected_slope)
    assert np.isclose(curve["Intercept"], expected_intercept)


def test_selected_standard_without_numeric_cq_and_duplicate_concentration_are_rejected():
    raw = pd.DataFrame([
        {"Plate": "P1", "Well": "A1", "Gene": "G", "Type": "Standard", "Label": "S1", "Cq": 20},
        {"Plate": "P1", "Well": "A2", "Gene": "G", "Type": "Standard", "Label": "S2", "Cq": 22},
        {"Plate": "P1", "Well": "A3", "Gene": "G", "Type": "Standard", "Label": "S3", "Cq": "Undetermined"},
    ])
    wells = prepare_wells(raw)
    standards = make_standard_map(wells, top=16, factor=4)
    curves, _ = fit_curves(wells, standards)
    assert not bool(curves.iloc[0]["Eligible"])
    assert "no retained numeric Cq" in curves.iloc[0]["Issue"]
    standards.loc[1, "Concentration"] = standards.loc[0, "Concentration"]
    duplicated, _ = fit_curves(wells, standards)
    assert not bool(duplicated.iloc[0]["Eligible"])
    assert "Duplicate standard concentration" in duplicated.iloc[0]["Status"]


def test_analysis_uses_geometric_reference_quantity_and_equal_weight_calibrator_samples():
    wells = prepare_wells(_base_rows())
    standards = make_standard_map(wells, top=16, factor=4)
    result = analyze(
        wells, standards, references=["RefA", "RefB"], calibrators=["A", "B"],
        efficiency_confirmed=True,
    )
    assert set(result) == {"samples", "curves", "standard_points", "quantities", "standard_results", "ddcq_results", "issues"}
    normalized = result["standard_results"]
    ratios = normalized.loc[normalized.Endpoint.eq("Target/reference ratio")]
    folds = normalized.loc[normalized.Endpoint.eq("Calibrator-relative fold change")]
    assert len(ratios) == len(folds) == 2
    ratio_a = ratios.loc[ratios.SampleId.eq("A"), "Value"].iloc[0]
    ratio_b = ratios.loc[ratios.SampleId.eq("B"), "Value"].iloc[0]
    assert np.isclose(ratio_a, 4.0 / np.sqrt(2.0), rtol=0.01)
    assert np.isclose(ratio_b, 16.0 / 4.0, rtol=0.01)
    assert np.isclose(folds.Value.mean(), 1.0)
    assert folds.Eligible.all()
    assert folds.CalibratorN.eq(2).all()
    ddcq = result["ddcq_results"]
    assert ddcq.Endpoint.eq("Calibrator-relative fold change").all()
    assert ddcq.Eligible.all()


def test_undetermined_target_stays_visible_and_out_of_range_quantity_is_withheld_by_default():
    raw = _base_rows()
    raw["Cq"] = raw["Cq"].astype(object)
    raw.loc[(raw.SampleId == "B") & (raw.Gene == "Target"), "Cq"] = "Undetermined"
    wells = prepare_wells(raw)
    standards = make_standard_map(wells, top=16, factor=4)
    result = analyze(wells, standards, references=["RefA"], calibrators=["A"], efficiency_confirmed=True)
    row = result["quantities"].loc[(result["quantities"].SampleId == "B") & (result["quantities"].Gene == "Target")].iloc[0]
    assert pd.isna(row["Quantity"])
    assert "undetermined" in row["Status"].lower()
    assert row["NumericCount"] == 0

    # A detected target later than the lowest calibrator remains range-flagged.
    raw = _base_rows()
    raw.loc[(raw.SampleId == "B") & (raw.Gene == "Target"), "Cq"] = 28.0
    wells = prepare_wells(raw)
    result = analyze(wells, make_standard_map(wells, top=16, factor=4), references=["RefA"], calibrators=["A"], efficiency_confirmed=True)
    row = result["quantities"].loc[(result["quantities"].SampleId == "B") & (result["quantities"].Gene == "Target")].iloc[0]
    assert row["RangeStatus"] == "BelowRange"
    assert pd.isna(row["Quantity"])
    ratio = result["standard_results"].loc[(result["standard_results"].SampleId == "B") & (result["standard_results"].Endpoint == "Target/reference ratio")].iloc[0]
    assert not bool(ratio["Eligible"])
    assert "outside" in ratio["Status"].lower()
    ddcq = result["ddcq_results"].loc[result["ddcq_results"].SampleId.eq("B")].iloc[0]
    assert not bool(ddcq["Eligible"])
    assert "outside calibrated range" in ddcq["Status"].lower()
    allowed = analyze(
        wells, make_standard_map(wells, top=16, factor=4), references=["RefA"],
        calibrators=["A"], allow_extrapolation=True, efficiency_confirmed=True,
    )
    provisional_ddcq = allowed["ddcq_results"].loc[allowed["ddcq_results"].SampleId.eq("B")].iloc[0]
    assert bool(provisional_ddcq["Eligible"])
    assert "provisional outside calibrated range" in provisional_ddcq["Status"].lower()


def test_extrapolated_reference_propagates_provisional_status_and_plate_trace():
    raw = _base_rows()
    raw.loc[(raw.SampleId == "A") & (raw.Gene == "RefA"), "Cq"] = 26.0
    wells = prepare_wells(raw)
    standards = make_standard_map(wells, top=16, factor=4)
    withheld = analyze(wells, standards, references=["RefA", "RefB"], calibrators=["B"], efficiency_confirmed=True)
    row = withheld["standard_results"].loc[(withheld["standard_results"].SampleId == "A") & (withheld["standard_results"].Endpoint == "Target/reference ratio")].iloc[0]
    assert not bool(row["Eligible"])
    assert "Reference quantity unavailable" in row["Status"]
    provisional = analyze(wells, standards, references=["RefA", "RefB"], calibrators=["B"], allow_extrapolation=True, efficiency_confirmed=True)
    row = provisional["standard_results"].loc[(provisional["standard_results"].SampleId == "A") & (provisional["standard_results"].Endpoint == "Target/reference ratio")].iloc[0]
    assert bool(row["Eligible"])
    assert "provisional" in row["Status"].lower()
    assert "RefA:P1" in row["PlateTrace"]


def test_missing_reference_and_invalid_calibrator_have_explicit_statuses():
    raw = _base_rows().loc[lambda frame: ~((frame.SampleId == "B") & (frame.Gene == "RefB"))]
    wells = prepare_wells(raw)
    result = analyze(wells, make_standard_map(wells, top=16, factor=4), references=["RefA", "RefB"], calibrators=["A", "missing"], efficiency_confirmed=False)
    missing_ref = result["standard_results"].loc[(result["standard_results"].SampleId == "B") & (result["standard_results"].Endpoint == "Target/reference ratio")].iloc[0]
    assert not bool(missing_ref["Eligible"])
    assert "Missing reference gene(s): RefB" in missing_ref["Status"]
    folds = result["standard_results"].loc[result["standard_results"].Endpoint.eq("Calibrator-relative fold change")]
    assert folds["Value"].isna().all()
    assert (~folds["Eligible"]).all()
    ddcq = result["ddcq_results"]
    assert ddcq.Endpoint.eq("Calibrator-relative fold change").all()
    assert ddcq.Status.str.contains("calibrator", case=False).any()
    assert ddcq.Status.str.contains("efficiency assumptions not confirmed", case=False).all()


def test_analysis_excludes_unmapped_type_but_reports_it():
    raw = _base_rows()
    extra = raw.loc[(raw.SampleId == "A") & (raw.Gene == "Target")].iloc[0].copy()
    extra["Well"] = "A-unmapped"
    extra["Label"] = "mystery"
    extra["SampleId"] = "mystery"
    extra["Type"] = "new assay control"
    raw = pd.concat([raw, pd.DataFrame([extra])], ignore_index=True)
    result = analyze(prepare_wells(raw), make_standard_map(prepare_wells(raw), top=16, factor=4), references=["RefA"], calibrators=["A"], efficiency_confirmed=True)
    issue = result["issues"].loc[result["issues"].Issue.str.contains("Unrecognized Type", regex=False)]
    assert not issue.empty
    row = result["samples"].loc[result["samples"].SampleId.eq("mystery")].iloc[0]
    assert not bool(row["Eligible"])


def test_cross_gene_group_conflict_blocks_every_normalized_endpoint_for_sample():
    raw = _base_rows()
    raw.loc[(raw.SampleId == "A") & (raw.Gene == "RefA"), "Group"] = "Treatment"
    wells = prepare_wells(raw)
    result = analyze(
        wells, make_standard_map(wells, top=16, factor=4), references=["RefA"],
        calibrators=["B"], efficiency_confirmed=True,
    )
    summaries = result["samples"].loc[result["samples"].SampleId.eq("A")]
    assert not summaries["Eligible"].any()
    assert summaries["Status"].str.contains("Conflicting sample metadata").all()
    standard = result["standard_results"].loc[result["standard_results"].SampleId.eq("A")]
    assert not standard["Eligible"].any()
    ddcq = result["ddcq_results"].loc[result["ddcq_results"].SampleId.eq("A")]
    assert not ddcq["Eligible"].any()


def test_empty_reference_selection_never_produces_eligible_nan_endpoints():
    wells = prepare_wells(_base_rows())
    result = analyze(
        wells, make_standard_map(wells, top=16, factor=4), references=[],
        calibrators=["A"], efficiency_confirmed=True,
    )
    standard = result["standard_results"]
    assert not standard["Eligible"].any()
    assert standard["Value"].isna().all()
    assert standard["Status"].str.contains("No reference genes selected").all()
    ddcq = result["ddcq_results"]
    assert not ddcq["Eligible"].any()
    assert ddcq["Status"].str.contains("No reference genes selected").all()


def test_back_calculated_concentration_range_catches_an_inverted_low_quality_curve():
    # This deliberately poor curve matches a reviewer repro: the sample Cq is
    # inside the observed standard Cq span but its fitted concentration is
    # below the lowest included standard, so it must not appear in range.
    raw = pd.DataFrame([
        {"Plate": "P1", "Well": "A1", "Gene": "G", "Type": "Standard", "Label": "S1", "Cq": 25.0},
        {"Plate": "P1", "Well": "A2", "Gene": "G", "Type": "Standard", "Label": "S2", "Cq": 21.0},
        {"Plate": "P1", "Well": "A3", "Gene": "G", "Type": "Standard", "Label": "S3", "Cq": 24.0},
        {"Plate": "P1", "Well": "B1", "Gene": "G", "Type": "Sample", "Label": "mouse1", "Cq": 25.0},
    ])
    wells = prepare_wells(raw)
    standards = pd.DataFrame([
        {"Plate": "P1", "Gene": "G", "Label": "S1", "Concentration": 1.0, "Include": True, "Reason": ""},
        {"Plate": "P1", "Gene": "G", "Label": "S2", "Concentration": 10.0, "Include": True, "Reason": ""},
        {"Plate": "P1", "Gene": "G", "Label": "S3", "Concentration": 100.0, "Include": True, "Reason": ""},
    ])
    curves, _ = fit_curves(wells, standards)
    curve = curves.iloc[0]
    assert np.isclose(curve["Slope"], -0.5)
    assert curve["R2"] < 0.1
    summary = summarise_samples(wells)
    quantity = quantify(summary.loc[summary.Gene.eq("G") & summary.SampleId.str.startswith("auto:")], curves)
    row = quantity.iloc[0]
    assert row["ObservedRangeStatus"] == "InRange"
    assert row["ConcentrationRangeStatus"] == "BelowRange"
    assert row["RangeStatus"] == "BelowRange"
    assert not bool(row["Eligible"])
    assert pd.isna(row["Quantity"])


def test_excluded_standard_level_requires_a_specific_review_reason():
    raw = pd.DataFrame([
        {"Plate": "P1", "Well": "A1", "Gene": "G", "Type": "Standard", "Label": label, "Cq": cq}
        for label, cq in [("S1", 25.0), ("S2", 22.0), ("S3", 19.0)]
    ])
    wells = prepare_wells(raw)
    standards = make_standard_map(wells, top=16, factor=4)
    standards.loc[1, "Include"] = False
    standards.loc[1, "Reason"] = ""
    curves, points = fit_curves(wells, standards)
    assert "missing a review reason" in points.loc[points.Label.eq("S2"), "Issue"].iloc[0]
    assert not bool(curves.iloc[0]["Eligible"])

    included_map = make_standard_map(wells, top=16, factor=4)
    excluded_well = wells.copy()
    excluded_well.loc[excluded_well.Label.eq("S2"), "Keep"] = False
    excluded_well.loc[excluded_well.Label.eq("S2"), "Reason"] = ""
    curves, points = fit_curves(excluded_well, included_map)
    assert "Excluded standard well is missing a review reason" in points.loc[points.Label.eq("S2"), "Issue"].iloc[0]
    assert not bool(curves.iloc[0]["Eligible"])

def test_standard_missing_companion_is_advisory_but_invalid_text_blocks():
    for companion, eligible in [('Undetermined',True),('',True),('invalid',False)]:
        rows=[{'Plate':'P','Well':str(i),'Gene':'G','Label':label,'Type':'Standard','Cq':cq}
              for i,(label,cq) in enumerate([('S1',20),('S1',companion),('S2',22),('S3',24)])]
        wells=prepare_wells(pd.DataFrame(rows));standards=make_standard_map(wells,top=16,factor=4)
        curves,_=fit_curves(wells,standards)
        assert bool(curves.iloc[0].Eligible) is eligible


def test_curve_efficiency_overflow_is_an_invalid_curve_not_an_export_error():
    wells=prepare_wells(pd.DataFrame([{'Plate':'P','Well':str(i),'Gene':'G','Label':f'S{i}','Type':'Standard','Cq':20+i*.001} for i in range(3)]))
    curves,_=fit_curves(wells,make_standard_map(wells,top=1000,factor=10))
    assert not curves.iloc[0].Eligible
    assert 'efficiency' in curves.iloc[0].Status.lower()


def test_independent_reference_sets_keep_biological_samples_and_missing_results_separate():
    from qpcr_session import analyze_reference_sets, rows_for_reference_set
    from qpcr_storage import new_run

    raw = _base_rows()
    # Both samples share their instrument label, and have different well counts.
    is_sample = raw.Type.eq("Sample")
    raw.loc[is_sample, "Label"] = "Repeated label"
    raw = raw.loc[~(raw.SampleId.eq("B") & raw.Gene.eq("RefB"))]
    wells = prepare_wells(raw)
    standards = make_standard_map(wells, top=16, factor=4)
    settings = new_run()["settings"]
    settings["calibrators"] = ["A"]
    set_a = {"id": "a", "name": "18S only", "references": ["RefA"], "efficiency_confirmed": True}
    set_b = {"id": "b", "name": "CYC1 only", "references": ["RefB"], "efficiency_confirmed": True}
    settings["reference_sets"] = [set_a, set_b]
    settings["target_assays"] = ["Target"]

    output = analyze_reference_sets(wells, standards, settings)
    target = output["ddcq_results"].loc[lambda rows: rows.Gene.eq("Target")]
    assert sorted(target.SampleId) == ["A", "A", "B", "B"]
    assert target.Label.unique().tolist() == ["Repeated label"]
    first = rows_for_reference_set(target, "a")
    assert list(first["Eligible"]) == [True, True]
    np.testing.assert_allclose(first.Value, [1, 2])
    second = rows_for_reference_set(target, "b")
    assert list(second["Eligible"]) == [True, False]
    assert "Missing reference" in second.loc[second.SampleId.eq("B"), "Status"].iloc[0]
    assert second.ReferenceAssays.unique().tolist() == ["RefB"]
    assert target.groupby(["ReferenceSetId", "Endpoint", "Gene", "SampleId"]).size().eq(1).all()
    issues_a = rows_for_reference_set(output["issues"], "a")
    assert not issues_a.Gene.eq("RefB").any()
    assert set(output["standard_results"].Gene) == {"Target"}

    set_b["references"] = ["RefA", "RefB"]
    recalculated = analyze_reference_sets(wells, standards, settings)
    first_again = rows_for_reference_set(recalculated["ddcq_results"], "a")
    pd.testing.assert_frame_equal(
        first.reset_index(drop=True),
        first_again.loc[first_again.Gene.eq("Target")].reset_index(drop=True),
    )


def test_phase_one_results_open_as_original_set_without_recalculation():
    from qpcr_session import rows_for_reference_set
    from qpcr_storage import decode_package, encode_package, new_run, reference_sets

    run = new_run()
    run["settings"]["references"] = ["Original reference"]
    run["results"]["ddcq_results"] = [{"SampleId": "S1", "Value": 2.5, "Eligible": True}]
    package = encode_package(run)
    reopened = decode_package(package)
    assert reopened == run
    assert reference_sets(reopened["settings"])[0]["references"] == ["Original reference"]
    table = pd.DataFrame(reopened["results"]["ddcq_results"])
    assert rows_for_reference_set(table, "original").Value.tolist() == [2.5]
    assert rows_for_reference_set(table, "unknown").empty
