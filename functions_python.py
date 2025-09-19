import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gc
from CSCORE import CSCORE
from tqdm import tqdm
from scipy.sparse import issparse
from scipy.stats import spearmanr, pearsonr
from typing import List

def compute_expression_margins(
    adata,
    batch_size: int = 1000,
    celltype_col: str = "cell_type_reannotation",
    sample_type_col: str = "sample_type",
    gene_symbol_col: str = "GeneSymbol",
    cancer_label: str = "Cancer cell",
    normal_label: str = "normal",
    epithelial_label: str = "Epithelial cell",
    immune_labels: List[str] = ["T cell", "B cell", "Myeloid cell", "NK", "Plasma cell", "Neutrophil", "Mast cell"],
    eps: float = 1e-6
):
    # cell index sets
    cancer_idx = (adata.obs[celltype_col] == cancer_label).values
    epithelial_idx = ((adata.obs[celltype_col] == epithelial_label) & (adata.obs[sample_type_col] == normal_label)).values
    immune_idx_dict = {ct: (adata.obs[celltype_col] == ct).values for ct in immune_labels}
    gene_names = adata.var[gene_symbol_col].values
    n_genes = adata.n_vars
    results = []
    for start in tqdm(range(0, n_genes, batch_size), desc="Processing gene batches"):
        end = min(start + batch_size, n_genes)
        batch_genes = np.arange(start, end)
        # gene slice
        X_batch = adata[:, batch_genes].X
        if issparse(X_batch):
            X_batch = X_batch.toarray()
        # mean calculation
        mean_expr_cancer = X_batch[cancer_idx, :].mean(axis=0)
        mean_expr_epi = X_batch[epithelial_idx, :].mean(axis=0)
        immune_means = []
        for ct in immune_labels:
            idx = immune_idx_dict[ct]
            mean_expr = X_batch[idx, :].mean(axis=0)
            immune_means.append(mean_expr)
        immune_means = np.vstack(immune_means)  # shape: (#immune_types, #genes)
        max_immune = immune_means.max(axis=0)                         # shape: (n_genes,)
        max_immune_idx = immune_means.argmax(axis=0)                  # shape: (n_genes,)
        max_immune_celltypes = [immune_labels[i] for i in max_immune_idx]
        log2fc_epi = np.log2((mean_expr_cancer + eps) / (mean_expr_epi + eps))
        log2fc_immune = np.log2((mean_expr_cancer + eps) / (max_immune + eps))
        for i, gidx in enumerate(batch_genes):
            results.append({
                "Gene": gene_names[gidx],
                "mean_expr_cancer": mean_expr_cancer[i],
                "mean_expr_normal": mean_expr_epi[i],
                "max_mean_expr_immune": max_immune[i],
                "immune_celltype_max": max_immune_celltypes[i],
                "log2FC_cancer_vs_normal": log2fc_epi[i],
                "log2FC_cancer_vs_maximmune": log2fc_immune[i]
            })
        del X_batch, mean_expr_cancer, mean_expr_epi, immune_means, max_immune, max_immune_idx, max_immune_celltypes
        gc.collect()
    return pd.DataFrame(results)

def cell_specific_gene_expression(
    adata,
    gene: str,
    celltype_col: str = "cell_type_reannotation",
    celltype_label: str = "Cancer cell",
    gene_symbol_col: str = "GeneSymbol"
):
    if gene not in adata.var[gene_symbol_col].values: # Finding gene indices
        raise ValueError(f"Gene '{gene}' not found in adata.var[{gene_symbol_col}]")
    gene_idx = adata.var[gene_symbol_col].values.tolist().index(gene)
    if celltype_col not in adata.obs.columns: # Filtering cells
        raise ValueError(f"'{celltype_col}' not found in adata.obs")
    cell_mask = adata.obs[celltype_col].values == celltype_label
    if cell_mask.sum() == 0: 
        raise ValueError(f"No cells found for label '{celltype_label}' in column '{celltype_col}'")
    X = adata.X[cell_mask, gene_idx] # Extracting expression matrix
    if not isinstance(X, np.ndarray):
        X = X.toarray().flatten()
    else:
        X = X.flatten()
    return X

def compute_positivity(adata, gene1, gene2, gene_symbol_col="GeneSymbol"):
    gene_list = adata.var[gene_symbol_col].tolist()
    try:
        idx1 = gene_list.index(gene1)
        idx2 = gene_list.index(gene2)
    except ValueError as e:
        raise ValueError(f"Gene not found: {e}")
    # Extracting expression values
    x1 = adata.X[:, idx1]
    x2 = adata.X[:, idx2]
    if not isinstance(x1, np.ndarray):
        x1 = x1.toarray().flatten()
        x2 = x2.toarray().flatten()
    else:
        x1 = x1.flatten()
        x2 = x2.flatten()
    # binary positivity
    gene1_pos = x1 > 0
    gene2_pos = x2 > 0
    double_pos = gene1_pos & gene2_pos
    n_total = len(x1)
    n_gene1_pos = gene1_pos.sum()
    n_gene2_pos = gene2_pos.sum()
    n_double_pos = double_pos.sum()
    pct_gene1 = n_gene1_pos / n_total
    pct_gene2 = n_gene2_pos / n_total
    pct_double = n_double_pos / n_total
    cond_prob1 = n_double_pos / n_gene1_pos if n_gene1_pos > 0 else np.nan
    cond_prob2 = n_double_pos / n_gene2_pos if n_gene2_pos > 0 else np.nan
    stats = pd.DataFrame([{
        f"{gene1}_pos_rate": pct_gene1,
        f"{gene2}_pos_rate": pct_gene2,
        f"{gene1}_{gene2}_double_pos_rate": pct_double,
        f"P(double+|{gene1}+)": cond_prob1,
        f"P(double+|{gene2}+)": cond_prob2
    }])
    return stats
