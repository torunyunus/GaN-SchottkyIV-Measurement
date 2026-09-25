"""Exploratory external validation on the published graphene/p-Ge Schottky data.

Source: Liu (2024), https://doi.org/10.5281/zenodo.11481314,
Data.rar / fig e3 / fig e3b / fig e3b-<T>K.xml.

Usage: python external_ge_analysis.py PATH_TO_FIG_E3B_DIRECTORY OUT_DIRECTORY
The original Excel 2003 XML files must first be extracted from Data.rar.
"""

import csv
import itertools
import json
import re
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
from scipy.interpolate import PchipInterpolator


NS = {"s": "urn:schemas-microsoft-com:office:spreadsheet"}
EPSILONS = (0.05, 0.10, 0.15)


def load_curves(directory):
    curves = []
    for path in directory.glob("*fig e3b-?*K.xml"):
        temperature = int(re.search(r"-(\d+)K\.xml$", path.name).group(1))
        rows = []
        for row in ET.parse(path).findall(".//s:Row", NS):
            cells = [c.findtext("s:Data", namespaces=NS) for c in row.findall("s:Cell", NS)]
            if cells and cells[0] == "DataValue":
                rows.append((float(cells[1]), float(cells[2])))  # Vc [V], Ic [A]
        curves.append((temperature, np.asarray(rows, dtype=float)))
    curves.sort(key=lambda x: x[0])
    if len(curves) != 6:
        raise ValueError(f"Expected six temperatures in Extended Data Fig. 3b; found {len(curves)}")
    voltage = curves[0][1][:, 0]
    if any(len(z) != len(voltage) or not np.array_equal(z[:, 0], voltage) for _, z in curves):
        raise ValueError("The six voltage grids are not identical")
    currents = np.stack([z[:, 1] for _, z in curves])
    if not np.isfinite(currents).all() or np.any(currents == 0):
        raise ValueError("Non-finite or exactly zero currents found; log policy needs review")
    return np.array([t for t, _ in curves]), voltage, np.log10(np.abs(currents))


def predict(temperatures, log_current, indices, target, method):
    chosen = np.asarray(indices, dtype=int)
    t = temperatures[chosen]
    y = log_current[chosen]
    if method == "PCHIP":
        return PchipInterpolator(t, y, axis=0)(target)
    lower = np.flatnonzero(t < target)[-1]
    upper = np.flatnonzero(t > target)[0]
    weight = (target - t[lower]) / (t[upper] - t[lower])
    return (1 - weight) * y[lower] + weight * y[upper]


def mae(true, predicted, mask):
    return float(np.mean(np.abs(true[mask] - predicted[mask])))


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def analyze(source, output):
    output.mkdir(parents=True, exist_ok=True)
    temperatures, voltage, log_current = load_curves(source)
    # All published Vc samples are used in the primary analysis. Sensitivity
    # excludes the near-zero region, mirroring the manuscript's voltage audit.
    masks = {"primary_all_bias": np.ones(voltage.size, bool),
             "sensitivity_absV_ge_0p05": np.abs(voltage) >= 0.05,
             "sensitivity_forward_V_0p05_to_1": (voltage >= 0.05) & (voltage <= 1.0),
             "sensitivity_reverse_V_minus1_to_minus0p05": (voltage >= -1.0) & (voltage <= -0.05)}
    folds = []
    for target_index in range(1, len(temperatures) - 1):
        observed = [i for i in range(len(temperatures)) if i != target_index]
        for method in ("Linear", "PCHIP"):
            pred = predict(temperatures, log_current, observed, temperatures[target_index], method)
            folds.append({"method": method, "held_out_K": int(temperatures[target_index]),
                          **{name: mae(log_current[target_index], pred, mask)
                             for name, mask in masks.items()}})
    write_csv(output / "ge_loto.csv", folds, list(folds[0]))

    schedules = []
    for k in (3, 4, 5):
        for inner in itertools.combinations(range(1, len(temperatures) - 1), k - 2):
            selected = (0, *inner, len(temperatures) - 1)
            omitted = [i for i in range(1, len(temperatures) - 1) if i not in inner]
            for method in ("Linear", "PCHIP"):
                errors = [mae(log_current[i], predict(temperatures, log_current, selected,
                                                     temperatures[i], method), masks["primary_all_bias"])
                          for i in omitted]
                schedules.append({"method": method, "k": k,
                                  "measured_K": ";".join(map(str, temperatures[list(selected)])),
                                  "omitted_K": ";".join(map(str, temperatures[omitted])),
                                  "mean_curve_MAE_decade": float(np.mean(errors)),
                                  "worst_curve_MAE_decade": float(max(errors)),
                                  **{f"coverage_{eps:.2f}": float(np.mean(np.array(errors) <= eps))
                                     for eps in EPSILONS}})
    write_csv(output / "ge_schedules.csv", schedules, list(schedules[0]))

    summary = []
    for method in ("Linear", "PCHIP"):
        for k in (3, 4, 5):
            subset = [r for r in schedules if r["method"] == method and r["k"] == k]
            for eps in EPSILONS:
                summary.append({"method": method, "k": k, "threshold_decade": eps,
                                "n_schedules": len(subset),
                                "Rmean": float(np.mean([r["mean_curve_MAE_decade"] <= eps for r in subset])),
                                "Rall": float(np.mean([r["worst_curve_MAE_decade"] <= eps for r in subset])),
                                "median_coverage": float(np.median([r[f"coverage_{eps:.2f}"] for r in subset])),
                                "median_mean_MAE": float(np.median([r["mean_curve_MAE_decade"] for r in subset]))})
    write_csv(output / "ge_summary.csv", summary, list(summary[0]))
    provenance = {"source_doi": "10.5281/zenodo.11481314",
                  "source_figure": "Extended Data Fig. 3b of doi:10.1038/s41586-024-07785-3",
                  "source_channels": "Vc (V), Ic (A)", "response": "log10(abs(Ic/A))",
                  "temperatures_K": temperatures.tolist(), "bias_min_V": float(voltage.min()),
                  "bias_max_V": float(voltage.max()), "points_per_curve": int(voltage.size),
                  "primary_bias_mask": "all published voltage samples",
                  "temperature_endpoints_always_observed": True,
                  "model_selection": "none; frozen Linear and PCHIP interpolators",
                  "exploratory_status": "Source was identified and feasibility screened before formal analysis"}
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print("LOTO:")
    for r in folds:
        print(r)
    print("Schedule summary:")
    for r in summary:
        print(r)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    analyze(Path(sys.argv[1]), Path(sys.argv[2]))
