"""Test gene and CpG sex effects and flag pairs where both endpoints are sex-dependent."""

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from data_processing.methylation_data_preparation import beta_to_m
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


def benjamini_hochberg(p_values):
    adjusted_p_values = pd.Series(np.nan, index=p_values.index)
    testable = p_values.notna()
    if testable.any():
        adjusted_p_values[testable] = multipletests(p_values[testable].to_numpy(), method="fdr_bh")[1]
    return adjusted_p_values


def deseq2_genes_sex_effect_test(counts, metadata, fdr_alpha, n_cpus, genes_in_pairs):
    """Fit DESeq2 on all genes and test male/female effects in genes expression."""
    inference = DefaultInference(n_cpus=n_cpus)
    dds = DeseqDataSet(counts=counts, metadata=metadata[["sex"]], design="~sex", inference=inference, quiet=True)
    dds.deseq2()
    gene_statistics = DeseqStats(dds, contrast=["sex", "male", "female"], inference=inference, independent_filter=False, quiet=True) 
    gene_statistics.summary()
    gene_results = gene_statistics.results_df
    gene_results = gene_results.loc[gene_results.index.intersection(genes_in_pairs)]
    gene_results["padj"] = benjamini_hochberg(gene_results["pvalue"])
    gene_results["contrast"] = "male_minus_female"
    gene_results["n_male"] = (metadata["sex"] == "male").sum()
    gene_results["n_female"] = (metadata["sex"] == "female").sum()
    gene_results = gene_results.rename_axis("gene_id").reset_index()
    gene_results["sex_dependent"] = gene_results["padj"] <= fdr_alpha
    return gene_results


def cpgs_sex_effect_test(cpg_methylation, metadata, m_value_boundary_shift, fdr_alpha):
    """Compare male and female CpGs M-values."""
    samples_male = metadata.index[metadata["sex"] == "male"]
    samples_female = metadata.index[metadata["sex"] == "female"]
    cpg_m_values = beta_to_m(cpg_methylation, m_value_boundary_shift)
    m_male = cpg_m_values[samples_male].to_numpy()
    m_female = cpg_m_values[samples_female].to_numpy()
    welch_test = ttest_ind(m_male, m_female, axis=1, equal_var=False)
    cpg_results = pd.DataFrame({"cpg_id": cpg_methylation.index, "n_male": len(samples_male), "n_female": len(samples_female), "p_value": welch_test.pvalue,})
    cpg_results["padj"] = benjamini_hochberg(cpg_results["p_value"])
    cpg_results["sex_dependent"] = cpg_results["padj"] <= fdr_alpha
    return cpg_results


def sex_filter(pairs, counts, cpg_methylation, metadata, fdr_alpha, n_cpus, m_value_boundary_shift):
    pairs["both_sex_dependent"] = False
    sexes = metadata["sex"].unique()
    if "male" not in sexes or "female" not in sexes: # skip sex filter if only one sex is present
        return pairs, None, None
    genes_in_pairs = pairs["MF_ensembl"].unique()
    gene_results = deseq2_genes_sex_effect_test(counts, metadata, fdr_alpha, n_cpus, genes_in_pairs)
    cpg_results = cpgs_sex_effect_test(cpg_methylation, metadata, m_value_boundary_shift, fdr_alpha)
    gene_adjusted_p_values = gene_results.set_index("gene_id")["padj"]
    cpg_adjusted_p_values = cpg_results.set_index("cpg_id")["padj"]
    flagged_pairs = []
    for gene_id, cpg_id in zip(pairs["MF_ensembl"], pairs["CpG"]):
        gene_adjusted_p_value = gene_adjusted_p_values[gene_id]
        cpg_adjusted_p_value = cpg_adjusted_p_values[cpg_id]
        if gene_adjusted_p_value <= fdr_alpha and cpg_adjusted_p_value <= fdr_alpha:
            flagged_pairs.append(True)
        else:
            flagged_pairs.append(False)
    pairs["both_sex_dependent"] = np.array(flagged_pairs)
    return pairs, gene_results, cpg_results
