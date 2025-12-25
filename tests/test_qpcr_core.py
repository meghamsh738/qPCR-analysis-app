import numpy as np
import pandas as pd

from qpcr_core import coerce_columns, fit_standard_curve, quantify_samples, mark_outliers


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
