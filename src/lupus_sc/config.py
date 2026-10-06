"""Pre-specified constants for the whole analysis.

Every analysis decision that was fixed *before* looking at results lives here,
so notebooks and scripts never hard-code them. See the Decision Log in the
README for the reasoning behind each value.
"""
from pathlib import Path

import numpy as np

# --- Paths -------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"

# --- Data source -------------------------------------------------------------
CENSUS_VERSION = "2025-11-08"          # pinned CELLxGENE Census release
DATASET_ID = "218acb0f-9f2f-4f76-b90b-15a4b7c7f629"   # Perez et al. 2022 lupus dataset in Census
ORGANISM = "Homo sapiens"

# --- Analysis cohort -----------------------------------------------------------
SEX = "female"                                  # only 2 male controls exist
ETHNICITIES = ("Asian", "European American")    # other groups too small
DISEASE_STATES = ("managed", "na")              # managed SLE + healthy controls
DISEASE_LABEL = "systemic lupus erythematosus"
CONTROL_LABEL = "normal"

# Donor-level QC exclusions, with the reason recorded next to each donor
EXCLUDED_DONORS = {
    "HC-566": "456 cells in total (next lowest ~1,487); counts per cell ~3x below "
              "median in every cell type; extreme outlier in pseudobulk PCA",
}

# --- Cell types ----------------------------------------------------------------
MAJOR_CELL_TYPES = (
    "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell",
    "classical monocyte",
    "non-classical monocyte",
    "B cell",
    "natural killer cell",
    "conventional dendritic cell",
)

# --- Pseudobulk / differential expression -------------------------------------
MIN_CELLS_PER_PSEUDOBULK = 10       # sc-best-practices DE chapter
MIN_COUNTS_PER_PSEUDOBULK = 1_000   # sc-best-practices DE chapter
FILTER_BY_EXPR = dict(min_count=10, min_total_count=15, large_n=10, min_prop=0.7)
DE_DESIGN = "~ main_cohort + ethnicity + age_z + log_cpc_z + disease"
DE_FDR = 0.05
STRONG_LFC = 1.0                    # |log2FC| threshold for "large" effects

# Pre-specified interferon-stimulated genes (hypothesis fixed before DE)
ISG = ("IFI27", "IFI44L", "IFI44", "ISG15", "MX1",
       "RSAD2", "IFIT1", "IFIT3", "OAS1", "SIGLEC1")

# --- Compositional analysis (scCODA) -------------------------------------------
SCCODA_FORMULA = "disease + main_cohort + ethnicity + age_z"
SCCODA_FDR = 0.05
SCCODA_REFERENCES = ("automatic", "natural killer cell", "classical monocyte")

# --- Gene set enrichment --------------------------------------------------------
GENESET_MIN_SIZE = 15   # exclusive bounds, as in the sc-best-practices GSEA chapter
GENESET_MAX_SIZE = 500

# --- Machine learning -----------------------------------------------------------
SEED = 42
CV_SPLITS = 5
CV_REPEATS = 10
LOGREG_CS = np.logspace(-4, 2, 13)  # same grid as the drafts (exact values)
TOP_K_GENES = 500                   # in-fold SelectKBest for the expression model
EXPR_MIN_CPM = 1.0                  # label-free gene filter: CPM >= 1 ...
EXPR_MIN_DONOR_FRACTION = 0.2       # ... in >= 20% of donors