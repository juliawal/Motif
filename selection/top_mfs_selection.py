import pandas as pd
from utils.selection_utils import random_generator_for_equal_scores


def select_top_mfs(aggregated_edges, n_mfs, mf_score_method, run_seeds):
    edge_weights_grouped_by_mfs = aggregated_edges.groupby("MF_ensembl")["weight"]
    if mf_score_method == "sum":
        mf_scores = edge_weights_grouped_by_mfs.sum()
    elif mf_score_method == "max":
        mf_scores = edge_weights_grouped_by_mfs.max()
    elif mf_score_method == "mean":
        mf_scores = edge_weights_grouped_by_mfs.mean()
    elif mf_score_method == "degree":
        mf_scores = edge_weights_grouped_by_mfs.size()
    ranked_mfs = mf_scores.sort_values(ascending=False)
    mf_random_generator = random_generator_for_equal_scores(run_seeds)
    if n_mfs is None or len(ranked_mfs) <= n_mfs:
        final_mfs = ranked_mfs
        mfs_with_same_score_as_unselected_mfs = pd.Series(False, index=final_mfs.index)
    else:
        boundary_score = ranked_mfs.iloc[n_mfs - 1]
        mfs_above_boundary = ranked_mfs[ranked_mfs > boundary_score]
        mfs_on_boundary = ranked_mfs[ranked_mfs == boundary_score]
        slots_for_mfs_on_boundary = n_mfs - len(mfs_above_boundary)
        random_order = mf_random_generator.permutation(len(mfs_on_boundary))
        chosen_positions = random_order[:slots_for_mfs_on_boundary]
        chosen_boundary_mfs = mfs_on_boundary.iloc[chosen_positions] # draw remaining slots from mfs on the boundary
        final_mfs = pd.concat([mfs_above_boundary, chosen_boundary_mfs])
        mfs_with_same_score_as_unselected_mfs = ((final_mfs == boundary_score) & (len(mfs_on_boundary) > slots_for_mfs_on_boundary))
    return pd.DataFrame({"MF_ensembl": final_mfs.index.astype(str), "mfs_scores": final_mfs.to_numpy(), "rank": range(1, len(final_mfs) + 1), "mfs_with_same_score_as_unselected_mfs": mfs_with_same_score_as_unselected_mfs.to_numpy(),})
