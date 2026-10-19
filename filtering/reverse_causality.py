import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def pairs_with_cpgs_near_genes_tss(pairs, cpg_coordinates, gene_tss_annotations, reverse_tss_window):

    cpgs_coordinates = cpg_coordinates[["chrom", "pos"]]
    cpgs_coordinates = cpgs_coordinates.rename(columns={"chrom": "cpg_chrom"})

    annotated_pairs = pairs.join(cpgs_coordinates, on="CpG") 

    genes_tss = gene_tss_annotations[["chrom", "tss", "strand"]]
    genes_tss = genes_tss.rename(columns={"chrom": "gene_chrom"}) 

    pairs_with_all_tss = annotated_pairs.reset_index(names="pair_index")
    pairs_with_all_tss = pairs_with_all_tss.join(genes_tss, on="MF_ensembl") # add all transcript starts belonging to each gene.

    # Choose nearest TSS for each pair and add TSS to annotated_pairs
    distance = pairs_with_all_tss["pos"] - pairs_with_all_tss["tss"]
    for row_position, strand in enumerate(pairs_with_all_tss["strand"]):
        if strand == "-":
            distance.iloc[row_position] *= -1 # reverse the distance sign for each gene on the minus strand
    pairs_on_the_same_chromosome = pairs_with_all_tss["cpg_chrom"] == pairs_with_all_tss["gene_chrom"]
    pairs_with_all_tss["distance_to_tss"] = distance.where(pairs_on_the_same_chromosome)
    pairs_with_all_tss["absolute_distance"] = pairs_with_all_tss["distance_to_tss"].abs()
    nearest_tss = (pairs_with_all_tss.sort_values(["pair_index", "absolute_distance"]).drop_duplicates("pair_index"))
    for column in ("gene_chrom", "tss", "strand", "distance_to_tss"):
        annotated_pairs[column] = nearest_tss[column].to_numpy()

    annotated_pairs["near_genes_tss"] = (annotated_pairs["distance_to_tss"].abs() <= reverse_tss_window)
    pairs_near_tss = annotated_pairs.loc[annotated_pairs["near_genes_tss"], ["MF_ensembl", "CpG"]]
    return pairs_near_tss


def compute_spearman_correlations(pairs, gene_expression, cpg_methylation):
    n_samples = gene_expression.shape[1]
    correlations = []
    for gene_id, cpg_id in pairs[["MF_ensembl", "CpG"]].itertuples(index=False, name=None):
        expression_values = gene_expression.loc[gene_id].to_numpy()
        methylation_values = cpg_methylation.loc[cpg_id].to_numpy()
        rho = np.nan
        if n_samples < 3:
            correlation_status = "too_few_samples"
        elif np.ptp(expression_values) == 0:
            correlation_status = "constant_gene_expression"
        elif np.ptp(methylation_values) == 0:
            correlation_status = "constant_cpg_methylation"
        else:
            rho = spearmanr(expression_values, methylation_values).statistic
            correlation_status = "tested"
        correlations.append((rho, n_samples, correlation_status))
    return pd.DataFrame(correlations, index=pairs.index, columns=["correlation_rho", "correlation_n", "correlation_status"])


def reverse_causality_filter(pairs, gene_expression, cpg_methylation_blocks, cpg_coordinates, gene_tss_annotations, reverse_tss_window):
    pairs_near_tss = pairs_with_cpgs_near_genes_tss(pairs, cpg_coordinates, gene_tss_annotations, reverse_tss_window)
    pairs = pairs.copy()  # create result table with all original pairs
    pairs["near_genes_tss"] = pairs.index.isin(pairs_near_tss.index)
    
    # Calculate correlation for the pairs in each methylation block.
    correlations_by_pair = {}
    for cpg_methylation_block in cpg_methylation_blocks:
        cpgs_in_block = pairs_near_tss["CpG"].isin(cpg_methylation_block.index)
        pairs_in_block = pairs_near_tss.loc[cpgs_in_block]
        if pairs_in_block.empty:
            continue
        correlations = compute_spearman_correlations(pairs_in_block, gene_expression, cpg_methylation_block)
        for pair_index, rho, n_samples, correlation_status in correlations.itertuples(name=None):
            correlations_by_pair[pair_index] = (rho, n_samples, correlation_status)

    # Give every pair its correlation result.
    correlation_values = []
    sample_counts = []
    correlation_statuses = []
    for pair_index, near_tss in pairs["near_genes_tss"].items():
        if not near_tss:
            rho = np.nan
            n_samples = 0
            correlation_status = "not_evaluated_outside_tss_window"
        elif pair_index in correlations_by_pair:
            rho, n_samples, correlation_status = correlations_by_pair[pair_index]
        else:
            rho = np.nan
            n_samples = 0
            correlation_status = "not_in_matrix"

        correlation_values.append(rho)
        sample_counts.append(n_samples)
        correlation_statuses.append(correlation_status)

    pairs["correlation_rho"] = np.array(correlation_values, dtype=float)
    pairs["correlation_n"] = np.array(sample_counts, dtype=int)
    pairs["correlation_status"] = np.array(correlation_statuses, dtype=object)

    # Flag pairs with CpGs near its genes tss that are negatively correlated.
    negative_correlation = pairs["correlation_rho"] < 0
    pairs["reverse_causality_flag"] = pairs["near_genes_tss"] & negative_correlation
    return pairs
