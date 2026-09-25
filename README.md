# GaN Schottky sparse-temperature reproducibility package

This repository accompanies the manuscript **“Reliability-Aware Measurement Planning for Sparse-Temperature I-V Characterization of an Au/Ni/n-GaN Schottky Diode.”**

## Public reproducibility contents

### data/
- `Schottky_GaN_clean_master_5189.csv` — audited canonical analysis dataset containing 5,189 experimental I-V points across 19 temperatures from 40 to 400 K. Key columns include `T_K`, `V_V`, `I_A`, and `log10_I_A`.

### results/
- `GaN_data_audit_summary.csv` — legacy-to-canonical reconciliation summary.
- `GaN_reliability_revision_summary.csv` — exact finite-schedule PCHIP `R_mean` / `R_all` / coverage summary.
- `LOTO_Linear_PCHIP.csv` — per-temperature Linear and PCHIP leave-one-temperature-out MAEs reproduced from the canonical CSV.

### code/
- `analysis_reliability_v1.py` — reference implementation of LOTO and exhaustive endpoint-retaining PCHIP schedule analysis.
- `requirements.txt` — minimal Python dependencies.

## Exploratory external application

`EXPLORATORY_EXTERNAL_CASE.md` documents the separate graphene/p-Ge Schottky source, the reproduction command and limitations. The raw third-party data remain in the original authors' archive; this repository includes only analysis code and derived error summaries.

## Reproduce the LOTO results

Python 3 is required.

```bash
pip install -r code/requirements.txt

python code/analysis_reliability_v1.py \
  --csv data/Schottky_GaN_clean_master_5189.csv \
  --mode loto \
  --outdir results_recomputed
```

Expected Equal-T MAEs (rounding to manuscript precision):

- Linear: **0.0586 decade**
- PCHIP: **0.0567 decade**

## Recompute exhaustive PCHIP schedule results

```bash
python code/analysis_reliability_v1.py \
  --csv data/Schottky_GaN_clean_master_5189.csv \
  --mode exhaustive \
  --outdir exhaustive_recomputed
```

This enumerates all endpoint-retaining schedules for `k=5,...,18`, totaling **130,917 schedules**.

## Data provenance

The experimental I-V series originated from the Au/Ni/n-GaN device study by **Doğan and Elagöz**, *Physica E* 63, 186–192 (2014), DOI: 10.1016/j.physe.2014.04.019, and was later used in the machine-learning study by **Torun and Doğan**, *Superlattices and Microstructures* 160, 107062 (2021), DOI: 10.1016/j.spmi.2021.107062.

A provenance audit of the archived legacy analysis matrix identified a missing first raw point at 40 K and a legacy duplication in which the 340 K curve repeated the 320 K curve. The canonical 5,189-point reconstruction used here resolves those legacy analysis-file issues without adding new experimental measurements. See `results/GaN_data_audit_summary.csv`.

The canonical CSV is the public dataset required to reproduce the analyses in this repository. Binary archival MAT files are not required to run the supplied reference analysis and are therefore not included in the public repository.

## Interpretation of reliability quantities

`R_mean` is the fraction of equally weighted admissible temperature schedules whose **mean curve-level MAE across omitted temperatures** is within the stated tolerance.

`R_all` is more conservative: every omitted temperature curve's MAE must satisfy the tolerance.

Neither quantity should be interpreted as a prospective probability for an unseen future device.

## Citation

Please cite the source experimental/device publications above and the associated measurement-planning manuscript when using this repository.
