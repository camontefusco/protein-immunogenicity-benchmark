import numpy as np

from pib.metrics import binary_classification_report, expected_calibration_error


def test_expected_calibration_error_is_bounded() -> None:
    y_true = np.array([0, 0, 1, 1])
    y_prob = np.array([0.1, 0.2, 0.8, 0.9])

    ece = expected_calibration_error(y_true, y_prob)

    assert 0.0 <= ece <= 1.0


def test_binary_classification_report_has_core_metrics() -> None:
    y_true = np.array([0, 0, 1, 1])
    y_prob = np.array([0.1, 0.3, 0.7, 0.9])

    report = binary_classification_report(y_true, y_prob)

    assert set(report) == {
        "auroc",
        "average_precision",
        "accuracy",
        "expected_calibration_error",
    }

