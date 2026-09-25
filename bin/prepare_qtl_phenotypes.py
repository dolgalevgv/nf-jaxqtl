#!/usr/bin/env python

from pathlib import Path
import re
import sys

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd


adata_path, regions_path = [Path(path) for path in sys.argv[1:3]]
group_cols = [col.strip() for col in sys.argv[3].split(",")]
if not all(group_cols):
    raise ValueError("group_col must contain non-empty column names")

adata = ad.read_h5ad(adata_path)
adata.var = adata.var.set_index("gene_id")
regions = pd.read_csv(regions_path, sep="\t")

groups = adata.obs.groupby(group_cols, observed=True, dropna=False, sort=False)
written_names = set()

for values, row_indices in groups.indices.items():
    if not isinstance(values, tuple):
        values = (values,)

    labels = []
    for value in values:
        value = "NA" if pd.isna(value) else str(value)
        labels.append(re.sub(r"[^A-Za-z0-9._=-]+", "_", value))

    name = "__".join(labels)
    if name in written_names:
        raise ValueError(f"Different groups produce the same filename: {name}.h5ad")
    written_names.add(name)
    subset = adata[row_indices, :].copy()

    pdata = dc.pp.pseudobulk(
        subset,
        sample_col='donor',
        groups_col=None,
        layer='counts',
        mode='sum'
    )

    counts = pd.DataFrame(pdata.X.T, index=pdata.var_names, columns=pdata.obs_names).astype(int)
    counts = regions.join(counts, on="gene_id", how="inner")
    offset = pd.DataFrame({"iid": pdata.obs["donor"].to_list(), "offset": np.log(pdata.obs["psbulk_counts"].to_numpy())})

    counts.to_csv(f"{name}.phenotype.bed.gz", sep="\t", index=False)
    offset.to_csv(f"{name}.offset.tsv", sep="\t", index=False)