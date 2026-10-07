"""Pathway-level interpretation of DE results.

Two complementary methods, both run on the DESeq2 Wald statistic of each
cell type separately (sc-best-practices GSEA chapter):
  - GSEA on REACTOME gene sets (MSigDB collection), permutation-based
  - PROGENy pathway activities with a univariate linear model (ULM)

Gene set resources are downloaded once (via OmniPath) and cached as CSV in
results/gsea/, so later runs use exactly the same gene sets.
"""
import decoupler as dc
import pandas as pd

from lupus_sc import config

GSEA_PERMUTATIONS = 1000   # decoupler default; smallest reportable p is ~1/1000
GSEA_DIR = config.RESULTS_DIR / "gsea"


def load_reactome(cache=GSEA_DIR / "reactome_genesets_used.csv"):
    """REACTOME gene sets with size strictly between GENESET_MIN/MAX_SIZE.

    Reads the cached file if it exists; otherwise downloads MSigDB, filters,
    and writes the cache.
    """
    if cache.exists():
        return pd.read_csv(cache)
    msigdb = dc.op.resource("MSigDB", organism="human")
    net = msigdb.query("collection == 'reactome_pathways'")
    net = net[~net.duplicated(("geneset", "genesymbol"))]
    net = net.rename(columns={"geneset": "source", "genesymbol": "target"})[["source", "target"]]
    size = net["source"].map(net.groupby("source").size())
    net = net[(size > config.GENESET_MIN_SIZE) & (size < config.GENESET_MAX_SIZE)]
    cache.parent.mkdir(parents=True, exist_ok=True)
    net.to_csv(cache, index=False)
    return net


def load_progeny(cache=GSEA_DIR / "progeny_used.csv", top=100):
    """PROGENy footprint genes (top responsive genes per pathway), cached."""
    if cache.exists():
        return pd.read_csv(cache)
    net = dc.op.progeny(organism="human", top=top)
    cache.parent.mkdir(parents=True, exist_ok=True)
    net.to_csv(cache, index=False)
    return net


def gene_stats(res, name):
    """Turn one DE result table into a 1 x genes matrix of Wald statistics."""
    d = (res.dropna(subset=["stat", "gene"])
            .drop_duplicates("gene")
            .set_index("gene")[["stat"]].T)
    d.index = [name]
    return d


def run_enrichment(de_results, reactome, progeny):
    """Run GSEA and PROGENy for each cell type.

    de_results: dict {cell_type: DE result table from de.run_de}
    Returns a dict of four DataFrames (pathway x cell type):
    gsea_score, gsea_padj, progeny_score, progeny_padj.
    """
    out = {k: {} for k in ("gsea_score", "gsea_padj", "progeny_score", "progeny_padj")}
    for ct, res in de_results.items():
        data = gene_stats(res, ct)
        s, p = dc.mt.gsea(data=data, net=reactome,
                          times=GSEA_PERMUTATIONS, seed=config.SEED)
        out["gsea_score"][ct], out["gsea_padj"][ct] = s.iloc[0], p.iloc[0]
        s, p = dc.mt.ulm(data=data, net=progeny)
        out["progeny_score"][ct], out["progeny_padj"][ct] = s.iloc[0], p.iloc[0]
    return {k: pd.DataFrame(v) for k, v in out.items()}