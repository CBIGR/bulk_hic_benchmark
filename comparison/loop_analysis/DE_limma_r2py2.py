import numpy as np
import pandas as pd

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter

pandas2ri.activate()

# Import R packages
limma = importr("limma")
stats = importr("stats")
base = importr("base")

def limma_de(expr_df: pd.DataFrame,
             group: pd.Series,
             contrast_name="KD_vs_CTRL",
             p_adj_threshold=0.05) -> pd.DataFrame:
    """
    Run limma differential expression.
    expr_df: genes x samples matrix (continuous values, often log2 expression)
    group: sample group labels (length == number of columns of expr_df)
    """

    # --- checks ---
    assert list(expr_df.columns) == list(group.index), "group index must match expr_df sample columns"

    # Convert Python data -> R
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_expr = ro.conversion.py2rpy(expr_df)
        r_group = ro.FactorVector(group.astype(str).values)

    ro.globalenv["expr"] = r_expr
    ro.globalenv["group"] = r_group

    # R code
    ro.r("""
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
    fit <- lmFit(expr, design)

    # Contrast: group2 - group1
    # e.g. KD - CTRL
    """)
    # Build contrast string automatically from group values
    levels = list(group.astype(str).unique())

    if len(levels) != 2:
        raise ValueError(f"limma_de expects exactly 2 groups, got {levels}")

    g1, g2 = levels[0], levels[1]
    contrast = f"{g2}-{g1}"   # group2 vs group1

    ro.globalenv["contrast_str"] = ro.StrVector([contrast])

    ro.r("""
    cont.matrix <- makeContrasts(contrasts=contrast_str, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2)
    res <- topTable(fit2, adjust = "BH", number=Inf,sort.by="none")
    """)

    # Extract result
    r_res = ro.globalenv["res"]
    with localconverter(ro.default_converter + pandas2ri.converter):
        res_df = ro.conversion.rpy2py(r_res)

    # Add gene IDs (rownames)
    gene_ids = list(expr_df.index)
    res_df.insert(0, "gene", gene_ids)

    # label DE genes
    res_df["label"] = (res_df["adj.P.Val"] < p_adj_threshold).astype(int)

    return res_df


# ---------------- Example usage ----------------
if __name__ == "__main__":
    expr_df = pd.read_csv("/media/thu/CN2033-DATA/phd/evaluation_dl_hic/loop_caller/knockdown/DataSet_04_049.txt", sep='\t', index_col=0)   # genes x samples

    # expr_df.rename(columns={expr_df.columns[0]: "id"}, inplace=True)
    # Detect knockoff and control columns
    knockoff_cols = [col for col in expr_df.columns if col.startswith("k")]
    control_cols = [col for col in expr_df.columns if col.startswith("c")]

    expr_df = expr_df[control_cols+knockoff_cols]
    group = pd.Series(['CTRL'] * len(control_cols) + ['KD'] * len(knockoff_cols), index=expr_df.columns)
    # group.index = expr_df.drop('id',axis = 1).columns

    res = limma_de(expr_df, group)
