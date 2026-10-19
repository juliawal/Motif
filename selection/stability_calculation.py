import numpy as np


def kuncheva_stability(mf_sets_per_split, n_genes_in_universe, n_mfs):
    expected_mf_overlap = n_mfs * n_mfs / n_genes_in_universe
    stability_values_per_split = []
    for mf_sets in mf_sets_per_split:
        for mfs in mf_sets:
            if len(mfs) < n_mfs:
                mean_kuncheva_stability = np.nan
                enough_mfs = False
                return mean_kuncheva_stability, enough_mfs
        observed_mf_overlap = len(mf_sets[0] & mf_sets[1])
        stability_values_per_split.append((observed_mf_overlap - expected_mf_overlap) / (n_mfs - expected_mf_overlap))
    mean_kuncheva_stability = np.mean(stability_values_per_split)
    enough_mfs = True
    return mean_kuncheva_stability, enough_mfs
