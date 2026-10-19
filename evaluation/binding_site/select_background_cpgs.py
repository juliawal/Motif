from utils.selection_utils import random_generator_for_equal_scores
import numpy as np
import pandas as pd


def build_cpg_matching_categories(cpg_annotations, gc_fraction_by_cpg, n_gc_bins=10):
    gc_fraction_by_cpg = gc_fraction_by_cpg.reindex(cpg_annotations.index)
    cpg_has_gc_fraction = gc_fraction_by_cpg.notna()
    cpg_matching_categories = cpg_annotations.loc[cpg_has_gc_fraction].copy()
    gc_fraction_by_cpg = gc_fraction_by_cpg.loc[cpg_has_gc_fraction]
    if gc_fraction_by_cpg.nunique() <= 1: # all CpGs have same gc_fraction
        cpg_matching_categories["gc_bin"] = 0
    else:
        cpg_matching_categories["gc_bin"] = pd.qcut(gc_fraction_by_cpg, q = n_gc_bins, labels = False, duplicates="drop").astype(int) # pd.qcut divides values into q bins with each bin containing aproximately the same number of observations (<< quantiles)
    cpg_matching_categories["matching_group"] = cpg_matching_categories["island_context"] + "," + cpg_matching_categories["gene_region"] + "," + cpg_matching_categories["gc_bin"].astype(str)
    return cpg_matching_categories


def draw_random_matched_background_cpgs(target_cpgs, cpg_matching_categories, mf_symbol, bg_to_target_ratio = 10, bg_sampling_seed = 42):
    """Draw bg_to_target_ratio (or all) control CpGs per MF from the same matching_group."""
    background_random_generator = random_generator_for_equal_scores([f"{bg_sampling_seed}:{mf_symbol}"])

    target_cpg_matching_group = cpg_matching_categories.loc[sorted(target_cpgs), "matching_group"]
    n_target_cpgs_by_matching_group = target_cpg_matching_group.value_counts().sort_index()

    cpg_is_target = cpg_matching_categories.index.isin(target_cpgs)
    background_cpgs_for_sampling = cpg_matching_categories.loc[~cpg_is_target] # exclude MF's targets
    background_cpgs_by_matching_group = {}
    for matching_group, background_cpgs_in_matching_group in background_cpgs_for_sampling.groupby("matching_group"):
        background_cpgs_by_matching_group[matching_group] = background_cpgs_in_matching_group.index.sort_values().to_numpy()
    background_cpgs = []
    for matching_group, n_target_cpgs in n_target_cpgs_by_matching_group.items():
        n_background_cpgs_needed = n_target_cpgs * bg_to_target_ratio 
        background_cpg_candidates = background_cpgs_by_matching_group.get(matching_group)
        if background_cpg_candidates is None:
            continue 
        if len(background_cpg_candidates) < n_background_cpgs_needed:
            sampled_background_cpgs = background_cpg_candidates.tolist()  # use all available CpGs if needed
        else:
            sampled_background_cpgs = background_random_generator.choice(background_cpg_candidates, size = n_background_cpgs_needed, replace = False).tolist()
        background_cpgs.extend(sampled_background_cpgs)
    return background_cpgs
