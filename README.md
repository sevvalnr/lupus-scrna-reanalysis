# Lupus scRNA-seq Reanalysis

**Donor-level, confounder-aware reanalysis of a 1.26M-cell lupus PBMC single-cell atlas.**

Published single-cell datasets are often reused as-is. This project asks a simple question first: *can the data be trusted for the comparison we want to make?* It then answers three biological questions about systemic lupus erythematosus (SLE) with methods that treat the **donor**, not the cell, as the unit of analysis.

`Python` · `Scanpy` · `CELLxGENE Census` · `scCODA` · `PyDESeq2` · `decoupler` · `scikit-learn`

---

## Highlights

- **Caught a data-format problem before analysis.** The ready-made `.h5ad` contained scaled values, not raw counts. Switched to raw UMI counts from a version-pinned CELLxGENE Census release and verified count integrity.
- **Audited the study design for confounding.** Processing cohort was strongly tied to disease status (one cohort contained only controls), with additional age, ancestry and sex imbalances. Defined a balanced analysis cohort of **225 donors (129 SLE / 96 controls)** with documented inclusion criteria.
- **Built a memory-efficient, resumable pipeline.** Donor × cell-type pseudobulk over ~900k cells via batched Census queries, running within free-tier limits (8 GB RAM), with built-in integrity checks.
- **Showed that a naive classifier would learn batch, not biology.** A model using only technical features (cohort + sequencing depth) reached **AUROC ≈ 0.82**. Biological features (cell-type composition + interferon score) reached **AUROC ≈ 0.85 on a held-out, unseen processing cohort** and generalized better than raw gene expression.

## Key findings

| Question | Method | Result |
|---|---|---|
| Do cell-type proportions change in SLE? | Compositional analysis (scCODA), checked against three reference cell types | Classical monocytes ↑; CD4 T, B and NK cells ↓, consistent across references |
| Which genes change, and in which cells? | Pseudobulk DE (PyDESeq2), GSEA (Reactome), PROGENy pathway activity | Strong interferon-stimulated gene (ISG) signature in all 7 major cell types; JAK-STAT activity up. Robust to cohort-restricted and low-depth-excluded reruns |
| Can donors be classified, and does it generalize? | Logistic regression, nested CV, donor-level splits, leave-cohort-out evaluation | Batch-only baseline ≈ 0.82; biological model ≈ 0.85 on an unseen cohort |
| How easily can leakage inflate results? | Feature selection on pure noise, inside vs. outside the CV folds | AUROC **1.000** (leaky) vs. **0.506** (correct) |


## Methods in brief

1. **Data** – Raw UMI counts from CELLxGENE Census (release `2025-11-08`), matched to author metadata via `observation_joinid`.
2. **Cohort definition** – Confounder audit (processing cohort, age, ancestry, sex) → balanced 225-donor cohort; all exclusions are pre-specified in `src/lupus_sc/config.py`.
3. **Pseudobulk** – Counts aggregated per donor × cell type, so statistics are computed across donors rather than across non-independent cells.
4. **Analyses** – Compositional analysis, differential expression with pathway follow-up, and machine learning with strict donor-level separation and a batch-only baseline.

The workflow follows the methodology of ***Single-cell best practices*** (Heumos, Schaar et al.), particularly the chapters on quality control, normalization, differential gene expression, compositional analysis and gene set enrichment. Decisions that were fixed before seeing results are recorded in `src/lupus_sc/config.py`.

## Repository structure

```
├── data/README.md        # how to obtain the data (no data files are committed)
├── scripts/              # Census download + pseudobulk construction
├── src/lupus_sc/         # reusable functions and pre-specified constants
├── notebooks/            # analysis notebooks 01–04
├── results/              # gene-, pathway- and model-level summary tables
└── figures/              # final figures
```

## Reproducing the analysis

```bash
git clone https://github.com/sevvalnr/lupus-scrna-reanalysis.git
cd lupus-scrna-reanalysis
pip install -r requirements.txt   # exact versions: requirements-lock.txt
pip install -e .

python scripts/build_pseudobulk_census.py   # builds donor-level inputs in data/
```

Then run the notebooks in `notebooks/` in order. A Colab-specific environment is provided in `requirements-colab.txt`.

## Data and privacy

All inputs are public (CC BY 4.0). **No donor-level data is redistributed**: count matrices can carry re-identification risk, so `data/` is excluded from version control and only summary-level results are shared. Controlled-access genotype data (dbGaP) was not used. See [`data/README.md`](data/README.md).

## References

- Perez, R. K. et al. Single-cell RNA-seq reveals cell type–specific molecular and genetic associations to lupus. *Science* **376**, eabf1970 (2022). https://doi.org/10.1126/science.abf1970
- Heumos, L., Schaar, A. C., Lance, C. et al. Best practices for single-cell analysis across modalities. *Nature Reviews Genetics* **24**, 550–572 (2023). https://doi.org/10.1038/s41576-023-00586-w — online book: https://www.sc-best-practices.org
- CZ CELLxGENE Discover Census. https://chanzuckerberg.github.io/cellxgene-census/
