import hashlib
import numpy as np
import pandas as pd


def filter_runs_dataframes(edges_per_run_files, pair_ids_after_filters):
    """Read the edge files and keep only the given MF–CpG pairs."""
    list_of_filtered_dataframes = []
    for edges_per_run_file in edges_per_run_files: 
        run_edges = pd.read_csv(edges_per_run_file, sep="\t", usecols=["MF_ensembl", "CpG", "importance"])
        filtered_run_edges = run_edges.rename(columns={"importance": "weight"})
        if pair_ids_after_filters is not None: 
            pairs_to_keep = [(mf, cpg) in pair_ids_after_filters for mf, cpg in zip(filtered_run_edges["MF_ensembl"], filtered_run_edges["CpG"])]
            filtered_run_edges = filtered_run_edges.loc[pairs_to_keep].reset_index(drop=True)
        list_of_filtered_dataframes.append(filtered_run_edges)
    return list_of_filtered_dataframes


def random_generator_for_equal_scores(run_seeds):
    run_seeds = sorted(str(run_seed) for run_seed in run_seeds)
    run_seeds_text = ",".join(run_seeds)
    seed_bytes = hashlib.blake2b(run_seeds_text.encode(), digest_size=8).digest()
    seed = int.from_bytes(seed_bytes)
    return np.random.default_rng(seed) 


def ordered_edges(edges):
    edges = edges[["MF_ensembl", "CpG", "weight"]]
    return edges.sort_values(["weight", "MF_ensembl", "CpG"], ascending=[False, True, True]).reset_index(drop=True)
