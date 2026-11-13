import json

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from qpcr_statistics import run_statistics


def _independent_results():
    rows = []
    values = {
        "A": [1.1, 1.8, 2.4, 3.0],
        "B": [2.5, 3.7, 4.2, 5.4],
        "C": [0.4, 1.2, 2.0, 2.9],
    }
    for group, observations in values.items():
        for index, value in enumerate(observations, start=1):
            rows.append({
                "SampleId": "{}-{}".format(group, index), "Group": group,
                "PairId": None, "Gene": "Target", "Endpoint": "Target/reference ratio",
                "Value": value, "Eligible": True, "Status": "Ready",
            })
    # Excluded observations remain visible in the input but do not affect n.
    rows.append({
        "SampleId": "A-excluded", "Group": "A", "PairId": None,
        "Gene": "Target", "Endpoint": "Target/reference ratio",
        "Value": np.nan, "Eligible": False, "Status": "Out of range",
    })
    return pd.DataFrame(rows)


def test_welch_matches_scipy_and_holm_adjusts_the_selected_family():
    results = _independent_results()
    omnibus, comparisons = run_statistics(
        results, "welch", [("A", "B"), ("A", "C"), ("B", "C")], group_order=["A", "B", "C"]
    )
    a = results.loc[results.Group.eq("A") & results.Eligible, "Value"].to_numpy()
    b = results.loc[results.Group.eq("B") & results.Eligible, "Value"].to_numpy()
    expected = stats.ttest_ind(a, b, equal_var=False)
    first = comparisons.loc[(comparisons.Group1 == "A") & (comparisons.Group2 == "B")].iloc[0]
    assert np.isclose(first["PValue"], expected.pvalue)
    assert first["Df"] > 0
    assert first["CI95Lower"] < first["EffectMeanDifference"] < first["CI95Upper"]
    assert first["n1"] == 4 and first["n2"] == 4
    assert omnibus["ExcludedN"] == 1
    assert comparisons["PAdjustedHolm"].ge(comparisons["PValue"]).all()
    assert json.dumps(omnibus, allow_nan=False)


def test_one_way_anova_and_planned_contrast_use_same_pooled_residual_error():
    results = _independent_results().loc[lambda frame: frame.Eligible].copy()
    omnibus, comparisons = run_statistics(results, "anova", [("A", "B"), ("A", "C")], group_order=["A", "B", "C"])
    arrays = [results.loc[results.Group.eq(group), "Value"].to_numpy() for group in ["A", "B", "C"]]
    expected = stats.f_oneway(*arrays)
    assert np.isclose(omnibus["Statistic"], expected.statistic)
    assert np.isclose(omnibus["PValue"], expected.pvalue)
    row = comparisons.loc[(comparisons.Group1 == "A") & (comparisons.Group2 == "B")].iloc[0]
    pooled = float(omnibus["MSE"])
    df = int(omnibus["DfError"])
    a = arrays[0]
    b = arrays[1]
    contrast_se = np.sqrt(pooled * (1 / len(a) + 1 / len(b)))
    contrast_t = (a.mean() - b.mean()) / contrast_se
    assert np.isclose(row["Statistic"], contrast_t)
    assert np.isclose(row["PValue"], 2 * stats.t.sf(abs(contrast_t), df))
    assert json.dumps(omnibus, allow_nan=False)


def _paired_results():
    # Intentionally shuffle biological rows and reverse some PairId order.
    rows = [
        ("A", "a-p2", "p2", 5.0), ("B", "b-p1", "p1", 3.0),
        ("A", "a-unmatched", "p4", 10.0), ("B", "b-p3", "p3", 7.0),
        ("A", "a-p1", "p1", 4.0), ("B", "b-p2", "p2", 4.5),
        ("A", "a-no-pair", None, 8.0),
        ("A", "a-p3", "p3", 8.0),
    ]
    return pd.DataFrame([
        {"Group": group, "SampleId": sample, "PairId": pair, "Gene": "Target",
         "Value": value, "Eligible": True, "Status": "Ready"}
        for group, sample, pair, value in rows
    ])


def test_paired_comparison_matches_by_pair_id_and_reports_unmatched_and_missing_ids():
    results = _paired_results()
    omnibus, comparisons = run_statistics(results, "paired", [("A", "B")], group_order=["A", "B"])
    row = comparisons.iloc[0]
    assert row["N"] == 3
    assert set(row["PairIds"]) == {"p1", "p2", "p3"}
    assert row["UnmatchedPairIds"] == ["p4"]
    assert row["MissingPairIdSampleIds"] == ["a-no-pair"]
    paired_a = np.array([4.0, 5.0, 8.0])
    paired_b = np.array([3.0, 4.5, 7.0])
    expected = stats.ttest_rel(paired_a, paired_b)
    assert np.isclose(row["PValue"], expected.pvalue)
    assert np.isclose(row["EffectMeanDifference"], np.mean(paired_a - paired_b))
    assert omnibus["Test"] == "Paired two-sample t-test"
    assert json.dumps(omnibus, allow_nan=False)
    shuffled = results.sample(frac=1, random_state=7)
    _, shuffled_comparisons = run_statistics(shuffled, "paired", [("A", "B")], group_order=["A", "B"])
    assert np.isclose(shuffled_comparisons.iloc[0]["Statistic"], row["Statistic"])


def test_duplicate_samples_pair_ids_and_multiple_endpoints_are_rejected():
    results = _independent_results().loc[lambda frame: frame.Eligible].copy()
    duplicate = pd.concat([results, results.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="one row per biological SampleId"):
        run_statistics(duplicate, "welch", [("A", "B")], group_order=["A", "B", "C"])

    paired = _paired_results()
    duplicate_pair = pd.concat([paired, paired.loc[paired.SampleId.eq("a-p1")].assign(SampleId="a-p1-duplicate")], ignore_index=True)
    with pytest.raises(ValueError, match="one biological sample per group"):
        run_statistics(duplicate_pair, "paired", [("A", "B")], group_order=["A", "B"])

    mixed_endpoint = pd.concat([results, results.assign(Endpoint="Calibrator-relative fold change", SampleId=results.SampleId + "-fold")], ignore_index=True)
    with pytest.raises(ValueError, match="Select one endpoint"):
        run_statistics(mixed_endpoint, "anova", [("A", "B")], group_order=["A", "B", "C"])

    mixed_gene = results.assign(Gene="Other")
    with pytest.raises(ValueError, match="Select one gene"):
        run_statistics(pd.concat([results, mixed_gene], ignore_index=True), "anova", [("A", "B")], group_order=["A", "B", "C"])


def test_missing_comparison_groups_and_insufficient_sample_sizes_are_rejected():
    results = _independent_results().loc[lambda frame: frame.Eligible].copy()
    small = results.loc[~results.Group.eq("C") | ~results.SampleId.isin(["C-2", "C-3", "C-4"])].copy()
    with pytest.raises(ValueError, match="at least two eligible biological samples"):
        run_statistics(small, "welch", [("A", "C")], group_order=["A", "B", "C"])
    with pytest.raises(ValueError, match="outside group_order"):
        run_statistics(results, "welch", [("A", "missing")], group_order=["A", "B", "C"])


def test_zero_variability_is_an_explicit_error_instead_of_a_plausible_p_value():
    rows = []
    for group, values in {"A": [1.0, 1.0, 1.0], "B": [2.0, 2.0, 2.0]}.items():
        for index, value in enumerate(values):
            rows.append({"SampleId": "{}{}".format(group, index), "Group": group, "Value": value, "Eligible": True})
    results = pd.DataFrame(rows)
    with pytest.raises(ValueError, match="zero or near-zero variability"):
        run_statistics(results, "welch", [("A", "B")], group_order=["A", "B"])
