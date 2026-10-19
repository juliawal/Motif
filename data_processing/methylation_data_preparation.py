"""Prepare CpG values and convert beta-values to M-values."""

import numpy as np
import pandas as pd


def keep_cg_cpgs(cpg_methylation):
    prefix = cpg_methylation.columns.str.slice(0, 2)
    is_cg_cpg = prefix == "cg"
    return cpg_methylation.loc[:, is_cg_cpg]


def drop_masked_cpgs(cpg_methylation, cpg_annotations_file):
    """Remove CpGs flagged by MASK_general in the Zhou manifest."""
    cpg_mask_annotation = pd.read_csv(cpg_annotations_file, sep="\t", usecols=["probeID", "MASK_general"]).rename(columns={"probeID": "cpg_id"})
    masked_cpg_ids = cpg_mask_annotation.loc[cpg_mask_annotation["MASK_general"].astype(bool), "cpg_id"]
    masked_cpg_ids = set(masked_cpg_ids)
    cpgs_not_masked = ~cpg_methylation.columns.isin(masked_cpg_ids)
    return cpg_methylation.loc[:, cpgs_not_masked]


def drop_incomplete_cpgs(cpg_methylation): 
    return cpg_methylation.loc[:, np.isfinite(cpg_methylation.to_numpy()).all(axis=0)]


def beta_to_m(cpg_methylation, m_value_boundary_shift):
    beta_values = cpg_methylation.to_numpy(copy=True)
    beta_values[beta_values == 0] = m_value_boundary_shift 
    beta_values[beta_values == 1] = 1 - m_value_boundary_shift
    cpg_m_values = np.log2(beta_values / (1 - beta_values))
    return pd.DataFrame(cpg_m_values, index=cpg_methylation.index, columns=cpg_methylation.columns)
