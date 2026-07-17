"""Model evaluation metrics."""

from __future__ import annotations

import numpy as np
from sklearn.metrics import accuracy_score, average_precision_score, roc_auc_score


def expected_calibration_error(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    *,
    n_bins: int = 10,
) -> float:
    """Compute binary expected calibration error with uniform probability bins."""
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    ece = 0.0

    for lower, upper in zip(bins[:-1], bins[1:]):
        in_bin = (y_prob > lower) & (y_prob <= upper)
        if not np.any(in_bin):
            continue
        bin_weight = np.mean(in_bin)
        bin_accuracy = np.mean(y_true[in_bin] == (y_prob[in_bin] >= 0.5))
        bin_confidence = np.mean(np.maximum(y_prob[in_bin], 1.0 - y_prob[in_bin]))
        ece += bin_weight * abs(bin_accuracy - bin_confidence)

    return float(ece)


def binary_classification_report(y_true: np.ndarray, y_prob: np.ndarray) -> dict[str, float]:
    """Return core metrics for an immunogenicity classifier."""
    y_pred = (y_prob >= 0.5).astype(int)
    return {
        "auroc": float(roc_auc_score(y_true, y_prob)),
        "average_precision": float(average_precision_score(y_true, y_prob)),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "expected_calibration_error": expected_calibration_error(y_true, y_prob),
    }

