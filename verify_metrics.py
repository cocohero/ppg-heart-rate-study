"""Validate published aggregate arithmetic; does NOT reproduce the private algorithm."""
import csv
import json
import math
from pathlib import Path
import statistics


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, label):
    require(math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-10),
            f"{label}: {actual} != {expected}")


def main():
    root = Path(__file__).resolve().parent
    meta = json.loads((root / "metrics.json").read_text(encoding="utf-8"))
    with (root / "subject_metrics.csv").open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    require(len(rows) == 36, "Expected 3 variants x 12 subjects")
    expected_subjects = {f"DATA_{i:02d}_TYPE{'01' if i == 1 else '02'}" for i in range(1, 13)}
    require(set(r["variant"] for r in rows) == set(meta["summary"]), "Variant mismatch")
    for variant, summary in meta["summary"].items():
        selected = [r for r in rows if r["variant"] == variant]
        require(len(selected) == 12, "Expected 12 rows per variant")
        require({r["subject"] for r in selected} == expected_subjects, "Subject mismatch")
        totals, counts, maes, within, successes = [], [], [], [], []
        for row in selected:
            total, valid, missing = (int(row[k]) for k in ("total_windows", "valid_windows", "missing_windows"))
            mae, rate = float(row["mae_bpm"]), float(row["within_5bpm"])
            require(0 < valid <= total and total - valid == missing, "Invalid counts")
            require(math.isfinite(mae) and mae >= 0 and 0 <= rate <= 1, "Invalid metric")
            hit_count = round(rate * valid)
            close(hit_count / valid, rate, "Integer threshold success count")
            totals.append(total)
            counts.append(valid)
            maes.append(mae)
            within.append(rate)
            successes.append(hit_count)
        total, valid, hits = sum(totals), sum(counts), sum(successes)
        require(summary["subjects"] == 12, "Subject count mismatch")
        require((summary["total_windows"], summary["valid_windows"], summary["missing_windows"]) ==
                (total, valid, total - valid), "Summary counts mismatch")
        close(statistics.mean(maes), summary["macro_mae_bpm"], "Macro MAE")
        close(sum(m * n for m, n in zip(maes, counts)) / valid,
              summary["pooled_mae_bpm"], "Pooled MAE")
        close(statistics.mean(within), summary["macro_within_5bpm"], "Macro threshold rate")
        close(hits / valid, summary["pooled_within_5bpm"], "Pooled threshold rate")
        close(hits / total, summary["within_5bpm_all_windows_missing_counted_as_failure"], "All-window rate")
        close(valid / total, summary["output_availability"], "Output availability")
        print(f"OK {variant}: MAE={statistics.mean(maes):.9f}, outputs={valid}/{total}")
    diagnostic = meta["threshold_only_loso"]
    folds = diagnostic["folds"]
    require(len(folds) == 12 and {f["held_out_subject"] for f in folds} == expected_subjects,
            "LOSO fold mismatch")
    close(statistics.mean(f["held_out_mae_bpm"] for f in folds), diagnostic["macro_mae_bpm"],
          "Threshold-only LOSO mean")
    print("OK aggregate arithmetic only; private algorithm NOT independently reproduced.")


if __name__ == "__main__":
    main()