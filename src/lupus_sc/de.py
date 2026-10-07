"""Pseudobulk differential expression (SLE vs. control) per cell type.

Workflow follows the sc-best-practices DE chapter: sum counts per donor and
cell type, drop small pseudobulk samples, filter lowly expressed genes, then
fit a negative binomial GLM (DESeq2) with technical and demographic covariates.
"""
import anndata as ad
import decoupler as dc
import numpy as np
import pertpy as pt

from lupus_sc import config


def load_pseudobulk(path=None):
    """Read the donor x cell type pseudobulk object and add shared covariates.

    Adds:
      - main_cohort as a categorical label ("cohort1" ... "cohort4")
      - age_z: age standardized over donors (one value per donor, not per sample)
      - total_counts: library size of each pseudobulk sample
    """
    pb = ad.read_h5ad(path or config.DATA_DIR / "pseudobulk_lupus.h5ad")
    pb.obs["main_cohort"] = "cohort" + pb.obs["main_cohort"].astype(int).astype(str)
    donor_age = pb.obs.drop_duplicates("donor_id")["age"]
    pb.obs["age_z"] = (pb.obs["age"] - donor_age.mean()) / donor_age.std()
    pb.obs["total_counts"] = np.asarray(pb.X.sum(axis=1)).ravel()
    return pb


def prepare_cell_type(pb, cell_type, exclude=()):
    """Subset one cell type, apply sample and gene filters.

    The gene filter depends on which samples are kept, so it is recomputed
    every time the donor set changes (e.g. in robustness analyses).
    """
    keep = (pb.obs["cell_type"] == cell_type) & ~pb.obs["donor_id"].isin(list(exclude))
    x = pb[keep].copy()
    x = x[(x.obs["n_cells"] >= config.MIN_CELLS_PER_PSEUDOBULK)
          & (x.obs["total_counts"] >= config.MIN_COUNTS_PER_PSEUDOBULK)].copy()
    if not isinstance(x.X, np.ndarray):
        x.X = x.X.toarray()
    dc.pp.filter_by_expr(adata=x, group="disease", **config.FILTER_BY_EXPR)
    return x


def add_depth_covariate(x):
    """Add sequencing depth (log counts per cell), z-scored within the cell type."""
    log_cpc = np.log(x.obs["total_counts"] / x.obs["n_cells"])
    x.obs["log_cpc"] = log_cpc
    x.obs["log_cpc_z"] = (log_cpc - log_cpc.mean()) / log_cpc.std()
    return x


def run_de(x, design=config.DE_DESIGN):
    """Fit DESeq2 and test SLE vs. control. Returns one row per gene,
    sorted by adjusted p-value, with gene symbols in the 'gene' column."""
    model = pt.tl.PyDESeq2(adata=x, design=design)
    model.fit()
    res = model.test_contrasts(model.contrast(
        column="disease",
        baseline=config.CONTROL_LABEL,
        group_to_compare=config.DISEASE_LABEL,
    ))
    res["gene"] = res["variable"].map(x.var["feature_name"])
    return res.sort_values("adj_p_value").reset_index(drop=True)


def summarize_de(res, x, cell_type):
    """One summary row for a cell type: sample counts and significant genes."""
    sig = res[res["adj_p_value"] < config.DE_FDR]
    return {
        "cell_type": cell_type,
        "n_samples": x.n_obs,
        "n_SLE": int((x.obs["disease"] == config.DISEASE_LABEL).sum()),
        "n_genes_tested": x.n_vars,
        "n_sig": len(sig),
        "n_up": int((sig["log_fc"] > 0).sum()),
        "n_down": int((sig["log_fc"] < 0).sum()),
        "n_sig_abs_lfc_gt1": int((sig["log_fc"].abs() > config.STRONG_LFC).sum()),
    }


def cell_type_tag(cell_type):
    """File-name friendly cell type name, e.g. 'classical_monocyte'."""
    return cell_type.replace(",", "").replace(" ", "_")