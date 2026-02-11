import numpy as np
import pandas as pd

from qpcr_core import coerce_columns, fit_standard_curve, quantify_samples, mark_outliers, suggest_standard_curve_exclusions


def test_coerce_columns_normalizes_and_casts():
    raw = pd.DataFrame(
        {
            "Plate": ["Plate 1"],
            "Well": ["A1"],
            "Gene": ["gapdh"],
            "Type": ["Sample"],
            "Sample ID": ["Mouse-1"],
            "Replicate": ["1"],
            "Ct": ["20.5"],
        }
    )

    out = coerce_columns(raw)

    assert "Cq" in out.columns
    assert out.loc[0, "Cq"] == 20.5
    assert out.loc[0, "Replicate"] == 1
    assert out.loc[0, "Label"] == "Mouse-1"


def test_fit_standard_curve_basic():
    std_df = pd.DataFrame(
        {
            "Gene": ["G1", "G1", "G1", "G1"],
            "Label": ["Std1", "Std1", "Std2", "Std2"],
            "Cq": [15.0, 15.2, 18.0, 18.1],
        }
    )
    mapping = pd.DataFrame({"Label": ["Std1", "Std2"], "Concentration": [100.0, 10.0]})

    curves, points = fit_standard_curve(std_df, mapping)

    assert curves.shape[0] == 1
    assert points.shape[0] == 2
    assert curves.loc[0, "slope"] < 0
    assert curves.loc[0, "R2"] > 0.9


def test_quantify_samples_uses_curve():
    samples = pd.DataFrame({"Gene": ["G1"], "Cq": [20.0]})
    curves = pd.DataFrame({"Gene": ["G1"], "slope": [-3.3], "intercept": [30.0]})

    out = quantify_samples(samples, curves)

    expected = 10 ** ((20.0 - 30.0) / -3.3)
    assert np.isclose(out.loc[0, "Quantity"], expected, rtol=1e-6)


def test_mark_outliers_flags_large_deltas():
    df = pd.DataFrame(
        {
            "Gene": ["G1", "G1", "G1"],
            "Label": ["L1", "L1", "L1"],
            "Cq": [20.0, 20.1, 28.0],
        }
    )

    out = mark_outliers(df, threshold=2.0)

    assert out["Outlier"].tolist() == [False, False, True]


def test_coerce_columns_normalizes_type_and_infers():
    raw = pd.DataFrame(
        {
            "Plate": ["Plate 1"] * 4,
            "Well": ["A1", "A2", "A3", "A4"],
            "Gene": ["g1"] * 4,
            "Type": ["Std", "Standard curve", "Unknown", ""],
            "Label": ["Std1", "Std2", "SampleA", "Std3"],
            "Replicate": [1, 1, 1, 1],
            "Cq": [18.0, 19.0, 20.0, 21.0],
        }
    )

    out = coerce_columns(raw)

    assert out["Type"].tolist() == ["Standard", "Standard", "Sample", "Standard"]


def test_suggest_standard_curve_exclusions_drops_single_outlier_replicate():
    mapping = pd.DataFrame(
        {
            "Label": [f"Std{i}" for i in range(1, 7)],
            "Concentration": [1000.0, 250.0, 62.5, 15.625, 3.90625, 0.9765625],
        }
    )

    slope = -3.3
    intercept = 35.0
    rows = []
    wells = ["A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2", "E1", "E2", "F1", "F2"]
    for i, label in enumerate(mapping["Label"].tolist(), start=0):
        conc = float(mapping.loc[mapping["Label"] == label, "Concentration"].iloc[0])
        x = np.log10(conc)
        base_cq = slope * x + intercept
        # Two replicates with small noise.
        cq1 = base_cq + 0.03
        cq2 = base_cq - 0.02
        # Make one replicate at Std3 a strong outlier.
        if label == "Std3":
            cq2 = base_cq + 4.0
        rows.append({"Plate": "P1", "Well": wells[2 * i], "Gene": "G1", "Label": label, "Cq": cq1})
        rows.append({"Plate": "P1", "Well": wells[2 * i + 1], "Gene": "G1", "Label": label, "Cq": cq2})

    std_df = pd.DataFrame(rows)
    excl, summary = suggest_standard_curve_exclusions(
        std_df,
        mapping,
        outlier_threshold=0.75,
        min_levels_to_keep=6,
        max_levels_to_drop=1,
        r2_target=0.98,
        eff_min=90.0,
        eff_max=110.0,
    )

    assert excl.shape[0] == 1
    assert excl.iloc[0]["Label"] == "Std3"
    assert excl.iloc[0]["reason"] == "replicate_outlier"
    assert summary.shape[0] == 1
    assert summary.loc[0, "suggested_R2"] > summary.loc[0, "baseline_R2"]


def test_suggest_standard_curve_exclusions_can_drop_whole_level_when_both_reps_bad():
    mapping = pd.DataFrame(
        {
            "Label": [f"Std{i}" for i in range(1, 7)],
            "Concentration": [1000.0, 250.0, 62.5, 15.625, 3.90625, 0.9765625],
        }
    )

    slope = -3.3
    intercept = 35.0
    rows = []
    wells = ["A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2", "E1", "E2", "F1", "F2"]
    for i, label in enumerate(mapping["Label"].tolist(), start=0):
        conc = float(mapping.loc[mapping["Label"] == label, "Concentration"].iloc[0])
        x = np.log10(conc)
        base_cq = slope * x + intercept
        cq1 = base_cq + 0.02
        cq2 = base_cq - 0.02
        # Std5: both replicates are wrong but consistent, so replicate-drop isn't allowed by ΔCq.
        if label == "Std5":
            cq1 = base_cq + 3.5
            cq2 = base_cq + 3.45
        rows.append({"Plate": "P1", "Well": wells[2 * i], "Gene": "G1", "Label": label, "Cq": cq1})
        rows.append({"Plate": "P1", "Well": wells[2 * i + 1], "Gene": "G1", "Label": label, "Cq": cq2})

    std_df = pd.DataFrame(rows)
    excl, summary = suggest_standard_curve_exclusions(
        std_df,
        mapping,
        outlier_threshold=0.75,
        min_levels_to_keep=5,
        max_levels_to_drop=1,
        r2_target=0.98,
        eff_min=90.0,
        eff_max=110.0,
    )

    # Entire level should be dropped (both wells for Std5).
    assert excl.shape[0] == 2
    assert set(excl["Label"].tolist()) == {"Std5"}
    assert set(excl["reason"].tolist()) == {"drop_level_fit"}
    assert summary.shape[0] == 1
    assert summary.loc[0, "dropped_levels"] == 1


def test_suggest_standard_curve_exclusions_returns_empty_when_curve_already_good():
    mapping = pd.DataFrame(
        {
            "Label": [f"Std{i}" for i in range(1, 7)],
            "Concentration": [1000.0, 250.0, 62.5, 15.625, 3.90625, 0.9765625],
        }
    )

    slope = -3.3
    intercept = 35.0
    rows = []
    wells = ["A1", "A2", "B1", "B2", "C1", "C2", "D1", "D2", "E1", "E2", "F1", "F2"]
    for i, label in enumerate(mapping["Label"].tolist(), start=0):
        conc = float(mapping.loc[mapping["Label"] == label, "Concentration"].iloc[0])
        x = np.log10(conc)
        base_cq = slope * x + intercept
        rows.append({"Plate": "P1", "Well": wells[2 * i], "Gene": "G1", "Label": label, "Cq": base_cq + 0.02})
        rows.append({"Plate": "P1", "Well": wells[2 * i + 1], "Gene": "G1", "Label": label, "Cq": base_cq - 0.02})

    std_df = pd.DataFrame(rows)
    excl, summary = suggest_standard_curve_exclusions(
        std_df,
        mapping,
        outlier_threshold=0.75,
        min_levels_to_keep=6,
        max_levels_to_drop=1,
        r2_target=0.98,
        eff_min=90.0,
        eff_max=110.0,
    )

    assert summary.shape[0] == 1
    assert excl.empty


def test_suggest_standard_curve_exclusions_best_r2_mode_can_use_non_outlier_single_rep():
    mapping = pd.DataFrame(
        {
            "Label": [f"Std{i}" for i in range(1, 7)],
            "Concentration": [1000.0, 250.0, 62.5, 15.625, 3.90625, 0.9765625],
        }
    )

    # Realistic profile where replicate disagreement is modest for several levels.
    # Conservative mode keeps all reps; aggressive mode (best_r2 + any_single_replicate)
    # should exclude some wells to improve linearity.
    std_df = pd.DataFrame(
        [
            {"Plate": "P1", "Well": "J1", "Gene": "Ccl2", "Label": "Std1", "Cq": 28.09552373},
            {"Plate": "P1", "Well": "J2", "Gene": "Ccl2", "Label": "Std1", "Cq": 28.87100386},
            {"Plate": "P1", "Well": "J3", "Gene": "Ccl2", "Label": "Std2", "Cq": 30.36116807},
            {"Plate": "P1", "Well": "J4", "Gene": "Ccl2", "Label": "Std2", "Cq": 30.24272422},
            {"Plate": "P1", "Well": "J5", "Gene": "Ccl2", "Label": "Std3", "Cq": 32.89609107},
            {"Plate": "P1", "Well": "J6", "Gene": "Ccl2", "Label": "Std3", "Cq": 32.45149972},
            {"Plate": "P1", "Well": "J7", "Gene": "Ccl2", "Label": "Std4", "Cq": 33.06982266},
            {"Plate": "P1", "Well": "J8", "Gene": "Ccl2", "Label": "Std4", "Cq": 33.78966724},
            {"Plate": "P1", "Well": "J9", "Gene": "Ccl2", "Label": "Std5", "Cq": 33.77495710},
            {"Plate": "P1", "Well": "J10", "Gene": "Ccl2", "Label": "Std5", "Cq": 34.95271561},
            {"Plate": "P1", "Well": "J11", "Gene": "Ccl2", "Label": "Std6", "Cq": 34.98309588},
            {"Plate": "P1", "Well": "J12", "Gene": "Ccl2", "Label": "Std6", "Cq": 34.82608403},
        ]
    )

    excl_default, summary_default = suggest_standard_curve_exclusions(
        std_df,
        mapping,
        outlier_threshold=0.75,
        min_levels_to_keep=6,
        max_levels_to_drop=1,
        r2_target=0.98,
        eff_min=90.0,
        eff_max=110.0,
    )
    assert excl_default.empty
    assert summary_default.shape[0] == 1

    excl_best_r2, summary_best_r2 = suggest_standard_curve_exclusions(
        std_df,
        mapping,
        outlier_threshold=0.75,
        min_levels_to_keep=6,
        max_levels_to_drop=0,
        r2_target=0.98,
        eff_min=90.0,
        eff_max=110.0,
        optimize_for="best_r2",
        replicate_drop_mode="any_single_replicate",
    )
    assert summary_best_r2.shape[0] == 1
    assert summary_best_r2.loc[0, "suggested_R2"] > summary_default.loc[0, "suggested_R2"]
    assert summary_best_r2.loc[0, "dropped_levels"] == 0
    assert summary_best_r2.loc[0, "dropped_replicates"] >= 1
    assert not excl_best_r2.empty
