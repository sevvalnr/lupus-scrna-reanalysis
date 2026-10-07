"""Cell-type composition: scCODA (Bayesian compositional model) and CLR.

Cell-type proportions are compositional: if one cell type expands, the others
must shrink. scCODA models the counts jointly and reports changes relative to
a reference cell type, so results are checked across several references
(sc-best-practices compositional analysis chapter).
"""
import re

import anndata as ad
import numpy as np
import pandas as pd
import pertpy as pt

from lupus_sc import config

COMPOSITION_DIR = config.RESULTS_DIR / "composition"
DISEASE_EFFECT = f"effect_df_disease[T.{config.DISEASE_LABEL}]"


def build_cell_table(cells, donors):
    """One row per cell with its donor, cell type and donor covariates.

    cells:  DataFrame with 'donor_id' and 'cell_type' (selected_cells.parquet)
    donors: DataFrame indexed by donor_id with disease, main_cohort,
            ethnicity and age (main_cohort_donors.csv)
    """
    d = donors[["disease", "main_cohort", "ethnicity", "age"]].copy()
    d["main_cohort"] = "cohort" + d["main_cohort"].astype(int).astype(str)
    d["age_z"] = (d["age"] - d["age"].mean()) / d["age"].std()
    obs = cells[["donor_id", "cell_type"]].astype(str).join(d, on="donor_id")
    assert obs["disease"].notna().all(), "cells from donors missing in donor table"
    obs.index = obs.index.astype(str)
    return ad.AnnData(obs=obs)   # labels only, no expression needed


def run_sccoda(cell_table, reference="automatic", num_samples=10_000,
               num_warmup=1_000, seed=config.SEED):
    """Fit scCODA with one reference cell type. Returns (model, mudata)."""
    model = pt.tl.Sccoda()
    data = model.load(
        cell_table, type="cell_level", generate_sample_level=True,
        cell_type_identifier="cell_type", sample_identifier="donor_id",
        covariate_obs=["disease", "main_cohort", "ethnicity", "age_z"],
    )
    data = model.prepare(data, modality_key="coda", formula=config.SCCODA_FORMULA,
                         reference_cell_type=reference)
    model.run_nuts(data, modality_key="coda", rng_key=seed,
                   num_samples=num_samples, num_warmup=num_warmup)
    model.set_fdr(data, config.SCCODA_FDR)
    return model, data


def effect_tables(data):
    """All covariate effect tables of a fitted model, keyed by covariate name."""
    varm = data["coda"].varm
    return {k: varm[k] for k in varm.keys() if k.startswith("effect_df")}


def save_effects(data, tag, out_dir=COMPOSITION_DIR):
    """Write every effect table as ref-<tag>__<covariate>.csv."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for k, df in effect_tables(data).items():
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", k)
        df.to_csv(out_dir / f"ref-{tag}__{safe}.csv")


def clr(counts, pseudocount=0.5):
    """Centered log-ratio transform of a donor x cell type count table.

    Used as ML features; the pseudocount (0.5) is the same as scCODA's.
    """
    log_counts = np.log(counts + pseudocount)
    return log_counts.sub(log_counts.mean(axis=1), axis=0)