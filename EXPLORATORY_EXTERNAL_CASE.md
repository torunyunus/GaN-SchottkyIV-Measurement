# Exploratory external Schottky I–V case

This folder addition accompanies Section 3.6 and Supporting Information Section S12 of the revised GaN measurement-planning manuscript. It applies the same complete-temperature holdout and finite-schedule definitions to an independent published graphene/p-Ge Schottky series. The source was located and feasibility screened before the formal calculation, so this is a post hoc exploratory application, not a blind or preregistered external validation.

## Original data

Chi Liu et al., *A hot-emitter transistor based on stimulated emission of heated carriers*, Nature 632, 782–787 (2024), https://doi.org/10.1038/s41586-024-07785-3. The original source data are in the authors' Zenodo record, https://doi.org/10.5281/zenodo.11481314, `Data.rar` / `fig e3` / `fig e3b` / six `fig e3b-<T>K.xml` files at 224, 232, 241, 251, 261 and 273 K. Download and extract them from that record. The original files are not redistributed in this repository.

## Reproduce the external analysis

Install the existing requirements in `code/requirements.txt`, then run:

```bash
python code/external_ge_analysis.py PATH_TO_DIRECTORY_WITH_SIX_XML_FILES external_recomputed
```

The script evaluates `log10(abs(Ic/A))` using the original voltage grid (1,001 points per curve, −3 to +3 V); excludes the two endpoint temperatures from whole-temperature LOTO; enumerates every endpoint-retaining schedule for k=3, 4, 5; and uses the same 0.05, 0.10, 0.15-decade thresholds as the GaN analysis. Restricted voltage windows are sensitivity checks. The checked-in derived outputs are:

- `results/external_ge_loto.csv`: fold-specific primary and sensitivity errors.
- `results/external_ge_schedules.csv`: all 14 schedules per interpolator with mean and worst-curve MAE and temperature coverage.
- `results/external_ge_summary.csv`: exact Rmean and Rall schedule fractions.
- `results/external_ge_provenance.json`: transformation and source identifiers.

The external 49 K span and 4–6 schedules per budget are much smaller than the 360 K GaN span and its exhaustive schedule populations. Numerical GaN budgets cannot be transferred from these results. The release DOI `10.5281/zenodo.22968314` identifies the earlier GaN repository release; these external files were proposed after that release and are not included in its archived snapshot.
