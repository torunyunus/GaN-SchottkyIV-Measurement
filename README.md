# GaN Schottky sparse-temperature reproducibility package

This package accompanies the manuscript **“Reliability-Aware Measurement Planning for Sparse-Temperature I-V Characterization of an Au/Ni/n-GaN Schottky Diode.”**

## Scope
The package contains the archived legacy analysis matrix, the audited canonical 5,189-point GaN dataset, a human-readable canonical CSV, audit/reliability summary outputs, and a Python reference implementation for the principal Linear/PCHIP whole-temperature reconstruction analysis.

## Files

### data/
- `legacy_data.mat` — archived 5,192-row legacy machine-learning matrix.
- `Schottky_GaN_clean_5189.mat` — audited/corrected canonical MAT dataset.
- `Schottky_GaN_clean_master_5189.csv` — audited canonical dataset in tabular form. Key columns: `T_K`, `V_V`, `I_A`, `log10_I_A`.

### results/
- `GaN_data_audit_summary.csv` — legacy-to-canonical reconciliation summary.
- `GaN_reliability_revision_summary.csv` — exact finite-schedule PCHIP R_mean/R_all/coverage summary used in the revised manuscript.
- `LOTO_Linear_PCHIP.csv` — per-temperature Linear and PCHIP LOTO MAEs reproduced from the canonical CSV.

### code/
- `analysis_reliability_v1.py` — reference implementation of LOTO and exhaustive endpoint-retaining PCHIP schedule analysis.

## Reproduce the LOTO results
Python 3 with `numpy`, `pandas`, and `scipy` is required.

```bash
python code/analysis_reliability_v1.py \
  --csv data/Schottky_GaN_clean_master_5189.csv \
  --mode loto \
  --outdir results_recomputed
```

Expected Equal-T MAEs (rounding to manuscript precision):
- Linear: 0.0586 decade
- PCHIP: 0.0567 decade

## Recompute exhaustive PCHIP schedule results

```bash
python code/analysis_reliability_v1.py \
  --csv data/Schottky_GaN_clean_master_5189.csv \
  --mode exhaustive \
  --outdir exhaustive_recomputed
```

This enumerates all endpoint-retaining schedules for k=5,...,18 (130,917 schedules total). Runtime depends on the machine.

## Provenance note
The experimental I-V series originated from the Au/Ni/n-GaN device study by Doğan and Elagöz (Physica E 63, 186–192, 2014; DOI: 10.1016/j.physe.2014.04.019) and was later used in the Torun and Doğan machine-learning study (Superlattices and Microstructures 160, 107062, 2021; DOI: 10.1016/j.spmi.2021.107062).

The legacy matrix contained a missing first raw point at 40 K and a duplicated 340 K curve corresponding to the 320 K data. The canonical 5,189-point reconstruction used in the present analysis corrects those legacy issues without adding new measurements.

## Important interpretation
`R_mean` is the fraction of equally weighted admissible schedules whose *mean curve-level MAE across omitted temperatures* is within the tolerance. `R_all` requires every omitted temperature curve's MAE to satisfy the same tolerance. Neither quantity is a prospective probability for an unseen future device.
