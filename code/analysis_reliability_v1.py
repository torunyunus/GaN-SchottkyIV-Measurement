#!/usr/bin/env python3
"""Reproducibility analysis for sparse-temperature GaN Schottky I-V reconstruction.

Inputs
------
Schottky_GaN_clean_master_5189.csv

Primary target
--------------
y = log10(I/A)

This script implements:
1) leave-one-temperature-out (LOTO) Linear and PCHIP reconstruction,
2) curve-level MAE under nearest-bracket/common-voltage support,
3) exhaustive endpoint-retaining PCHIP schedule enumeration for k=5,...,18,
4) finite-schedule R_mean, R_all, and temperature-level coverage.

The implementation intentionally mirrors the definitions reported in the manuscript.
"""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

EPS = (0.05, 0.10, 0.15)


def load_curves(csv_path: Path):
    df = pd.read_csv(csv_path)
    curves = {}
    for T, g in df.groupby("T_K", sort=True):
        g = g.sort_values("V_V")
        # Keep finite points and collapse duplicate voltages deterministically by mean.
        gg = g[["V_V", "log10_I_A"]].replace([np.inf, -np.inf], np.nan).dropna()
        gg = gg.groupby("V_V", as_index=False)["log10_I_A"].mean().sort_values("V_V")
        curves[int(T)] = (gg["V_V"].to_numpy(float), gg["log10_I_A"].to_numpy(float))
    temps = np.array(sorted(curves), dtype=int)
    return df, curves, temps


def interp_curve_at_v(curve, vq):
    v, y = curve
    return np.interp(vq, v, y)


def target_mask_and_v(curves, target_T: int, measured_temps):
    measured = np.array(sorted(measured_temps), dtype=int)
    lows = measured[measured < target_T]
    highs = measured[measured > target_T]
    if len(lows) == 0 or len(highs) == 0:
        return None
    TL, TU = int(lows.max()), int(highs.min())
    vt, yt = curves[target_T]
    vL, _ = curves[TL]
    vU, _ = curves[TU]
    lo = max(vL.min(), vU.min())
    hi = min(vL.max(), vU.max())
    mask = (vt >= lo - 1e-12) & (vt <= hi + 1e-12)
    return TL, TU, vt[mask], yt[mask]


def predict_linear(curves, target_T: int, measured_temps):
    sup = target_mask_and_v(curves, target_T, measured_temps)
    if sup is None:
        return None
    TL, TU, vt, yt = sup
    yL = interp_curve_at_v(curves[TL], vt)
    yU = interp_curve_at_v(curves[TU], vt)
    w = (target_T - TL) / (TU - TL)
    pred = yL + w * (yU - yL)
    return vt, yt, pred


def predict_pchip(curves, target_T: int, measured_temps):
    sup = target_mask_and_v(curves, target_T, measured_temps)
    if sup is None:
        return None
    TL, TU, vt, yt = sup
    measured = np.array(sorted(measured_temps), dtype=int)
    pred = np.full(vt.shape, np.nan, dtype=float)
    for j, vq in enumerate(vt):
        Ts = []
        ys = []
        for T in measured:
            v, y = curves[int(T)]
            if vq >= v.min() - 1e-12 and vq <= v.max() + 1e-12:
                Ts.append(float(T))
                ys.append(float(np.interp(vq, v, y)))
        Ts = np.asarray(Ts)
        ys = np.asarray(ys)
        # Need bracketing support and at least two T values. PCHIP with 2 points is linear.
        if len(Ts) >= 2 and Ts.min() <= target_T <= Ts.max():
            pred[j] = float(PchipInterpolator(Ts, ys, extrapolate=False)(target_T))
    ok = np.isfinite(pred) & np.isfinite(yt)
    return vt[ok], yt[ok], pred[ok]


def curve_mae(result):
    if result is None or len(result[1]) == 0:
        return np.nan
    _, yt, pred = result
    return float(np.mean(np.abs(pred - yt)))


def loto(curves, temps):
    rows = []
    interior = temps[1:-1]
    for T in interior:
        measured = [int(x) for x in temps if x != T]
        rows.append({
            "T_K": int(T),
            "Linear_MAE": curve_mae(predict_linear(curves, int(T), measured)),
            "PCHIP_MAE": curve_mae(predict_pchip(curves, int(T), measured)),
        })
    out = pd.DataFrame(rows)
    return out


def schedule_metrics(curves, temps, measured, eps_values=EPS):
    measured = tuple(sorted(int(x) for x in measured))
    omitted = [int(T) for T in temps[1:-1] if int(T) not in measured]
    errs = []
    for T in omitted:
        e = curve_mae(predict_pchip(curves, T, measured))
        if np.isfinite(e):
            errs.append(e)
        else:
            return None
    errs = np.asarray(errs, float)
    row = {
        "k": len(measured),
        "measured_T_K": ";".join(map(str, measured)),
        "n_omitted": len(errs),
        "E_mean": float(np.mean(errs)) if len(errs) else np.nan,
        "E_max_curve": float(np.max(errs)) if len(errs) else np.nan,
    }
    for eps in eps_values:
        tag = f"{eps:.2f}"
        row[f"coverage_eps_{tag}"] = float(np.mean(errs <= eps)) if len(errs) else np.nan
        row[f"pass_mean_eps_{tag}"] = int(row["E_mean"] <= eps) if len(errs) else 0
        row[f"pass_all_eps_{tag}"] = int(row["E_max_curve"] <= eps) if len(errs) else 0
    return row


def exhaustive(curves, temps, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    interior = [int(x) for x in temps[1:-1]]
    lo, hi = int(temps[0]), int(temps[-1])
    summary = []
    for k in range(5, 19):
        rows = []
        choose = k - 2
        for inner in itertools.combinations(interior, choose):
            measured = (lo, *inner, hi)
            r = schedule_metrics(curves, temps, measured)
            if r is not None:
                rows.append(r)
        d = pd.DataFrame(rows)
        d.to_csv(outdir / f"PCHIP_exhaustive_k{k}.csv", index=False)
        s = {"k": k, "N_schedules": len(d)}
        for eps in EPS:
            tag = f"{eps:.2f}"
            s[f"R_mean_eps_{tag}"] = d[f"pass_mean_eps_{tag}"].mean()
            s[f"R_all_eps_{tag}"] = d[f"pass_all_eps_{tag}"].mean()
            s[f"median_coverage_eps_{tag}"] = d[f"coverage_eps_{tag}"].median()
        summary.append(s)
        print(f"k={k}: {len(d):,} schedules")
    sm = pd.DataFrame(summary)
    sm.to_csv(outdir / "PCHIP_exhaustive_summary.csv", index=False)
    return sm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, default=Path("Schottky_GaN_clean_master_5189.csv"))
    ap.add_argument("--mode", choices=["loto", "exhaustive"], default="loto")
    ap.add_argument("--outdir", type=Path, default=Path("results_repro"))
    args = ap.parse_args()

    _, curves, temps = load_curves(args.csv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.mode == "loto":
        out = loto(curves, temps)
        out.to_csv(args.outdir / "LOTO_Linear_PCHIP.csv", index=False)
        print(out.to_string(index=False))
        print("Equal-T Linear MAE:", out["Linear_MAE"].mean())
        print("Equal-T PCHIP MAE:", out["PCHIP_MAE"].mean())
    else:
        out = exhaustive(curves, temps, args.outdir)
        print(out.to_string(index=False))


if __name__ == "__main__":
    main()
