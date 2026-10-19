import numpy as np


def calculate_degree_matched_p_value(proximity_scores_by_gene, mf_entrez_ids, sampled_gene_sets):
    mfs_scores = []
    for mf_entrez_id in mf_entrez_ids:
        mf_score = proximity_scores_by_gene.get(mf_entrez_id, 0.0) # use 0.0 if there is no entry in proximity_scores_by_gene for that MF
        mfs_scores.append(mf_score)
    mfs_score = np.mean(mfs_scores)
    n_sampled_gene_set_scores_greater_equal_mf_score = 0
    for sampled_gene_set in sampled_gene_sets:
        genes_scores = []
        for gene in sampled_gene_set: 
            gene_score = proximity_scores_by_gene.get(gene, 0.0)
            genes_scores.append(gene_score)
        sampled_gene_set_score = np.mean(genes_scores)
        if sampled_gene_set_score >= mfs_score:
            n_sampled_gene_set_scores_greater_equal_mf_score += 1
    n_random_sets = len(sampled_gene_sets)
    p_value = (n_sampled_gene_set_scores_greater_equal_mf_score + 1) / (n_random_sets + 1)
    return p_value
