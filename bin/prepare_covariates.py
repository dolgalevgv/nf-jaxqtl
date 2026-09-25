#!/usr/bin/env python

import sys
from pathlib import Path

import pandas as pd


donors_path, pca_path = [Path(path) for path in sys.argv[1:3]]

donors = pd.read_csv(donors_path, index_col=0)
pca = pd.read_csv(pca_path, sep="\t", index_col=0)

donors = donors.drop(columns="vcf_id").join(pca)
donors = donors.rename(columns={"donor_id": "iid"})

donors.to_csv("covariates.tsv", sep="\t")