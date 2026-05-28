def compute_metrics(
    predicted: set[tuple[int, int]],
    reference: set[tuple[int, int]],
) -> dict:
    tp = len(predicted & reference)
    fp = len(predicted - reference)
    fn = len(reference - predicted)

    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    ppv         = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1          = (2 * tp) / (2 * tp + fp + fn) if (2 * tp + fp + fn) > 0 else 0.0
    # simplified for RNA: TN is undefined, so MCC = TP / sqrt((TP+FP)(TP+FN))
    mcc_denom = ((tp + fp) * (tp + fn)) ** 0.5
    mcc = tp / mcc_denom if mcc_denom > 0 else 0.0

    return {
        "tp":          tp,
        "fp":          fp,
        "fn":          fn,
        "sensitivity": sensitivity,
        "ppv":         ppv,
        "f1":          f1,
        "mcc":         mcc,
    }
