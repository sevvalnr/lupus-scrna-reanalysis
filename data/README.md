# Data

No data files are committed to this repository. All inputs are public and can be
re-created with the code in this repository.

## Source
- **Study:** Perez, R. K. et al. *Single-cell RNA-seq reveals cell type–specific molecular
  and genetic associations to lupus.* Science 376, eabf1970 (2022).
  https://doi.org/10.1126/science.abf1970
- **Access:** CZ CELLxGENE Discover Census, release `2025-11-08`
  (raw UMI counts, `X["raw"]`). Author metadata (`Processing_Cohort`, `disease_state`)
  from the CELLxGENE h5ad of the same dataset, matched via `observation_joinid`.
- **License:** CC BY 4.0. Genotype data (dbGaP phs002812) are controlled-access and were
  **not** used.

## Why no data is shared
Donor-level count matrices can carry re-identification risk, so neither the raw data nor
the derived donor-level pseudobulk / feature tables are redistributed here.

## How to rebuild
1. Run `scripts/build_pseudobulk_census.py` (Google Colab recommended; ~12 GB RAM).
2. Place the outputs in this folder.
3. Run the notebooks in `notebooks/` in order.