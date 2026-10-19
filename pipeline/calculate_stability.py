from argparse import ArgumentParser
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from selection.top_edges_selection import select_top_edges_per_seed
from utils.selection_utils import filter_runs_dataframes, random_generator_for_equal_scores
from selection.run_aggregation import combine_and_aggregate_runs
from selection.top_mfs_selection import select_top_mfs
from selection.stability_calculation import kuncheva_stability


def create_runs_aggregation_methods(aggregation_options):
    runs_aggregation_methods = []
    for method_name in aggregation_options["runs_aggregation_method"]:
        if method_name == "freq":
            for alpha in aggregation_options["alpha"]:
                for beta in aggregation_options["beta"]:
                    runs_aggregation_methods.append({"runs_aggregation_method": method_name, "alpha": alpha, "beta": beta})
        else:
            runs_aggregation_methods.append({"runs_aggregation_method": method_name, "alpha": None, "beta": None})
    return runs_aggregation_methods


def create_configuration_grid(runs_aggregation_methods):
    configurations = []
    for mf_score_method in ("sum", "max", "mean"):
        for runs_aggregation_method in runs_aggregation_methods:
            if runs_aggregation_method["runs_aggregation_method"] == "freq" and runs_aggregation_method["alpha"] == 1 and runs_aggregation_method["beta"] == 0: # frequency with alpha = 1 and beta = 0 is the same as mean
                continue
            configuration = runs_aggregation_method.copy()
            configuration["mf_score_method"] = mf_score_method
            configuration_name = configuration["runs_aggregation_method"]
            if configuration_name == "freq":
                configuration_name += f"_a{configuration['alpha']}_b{configuration['beta']}"
            configuration["configuration_name"] = configuration_name + "_" + mf_score_method
            configurations.append(configuration)
    configurations.append({"runs_aggregation_method": "mean", "alpha": None, "beta": None, "mf_score_method": "degree", "configuration_name": "degree"}) # degree needs only one run aggregation setting
    return configurations



def complementary_splits(n_runs):
    if n_runs % 2:
        raise ValueError("Split-half stability needs an even number of runs.")
    run_splits = []
    for left_run_indices in combinations(range(n_runs), n_runs // 2):
        if 0 not in left_run_indices: # keep each split e.g. left: (0,1), right: (2,3) and left: (2,3), right: (0,1) only once, with run 0 in the left half
            continue
        right_run_indices = []
        for index in range(n_runs):
            if index not in left_run_indices:
                right_run_indices.append(index)
        run_splits.append((left_run_indices, tuple(right_run_indices)))
    return run_splits


def sample_splits(run_splits, n_splits, split_seed):
    split_random_generator = np.random.default_rng(split_seed)
    selected_split_indices = sorted(split_random_generator.choice(len(run_splits), size=n_splits, replace=False))
    return [run_splits[index] for index in selected_split_indices]


def calculate_stability_for_all_settings(list_of_filtered_dataframes, run_seeds, n_genes_in_universe, cohort, n_mfs, edge_cuts, configurations, selected_run_splits):
    stability_rows = []
    for k_edges in edge_cuts:
        cut_dataframes = []
        for dataframe, run_seed in zip(list_of_filtered_dataframes, run_seeds, strict=True):
            edge_random_generator = random_generator_for_equal_scores([run_seed])
            cut_dataframes.append(select_top_edges_per_seed(dataframe, k_edges, edge_random_generator))
        for configuration in configurations:
            runs_aggregation_method = {key: configuration[key] for key in ("runs_aggregation_method", "alpha", "beta")}
            mf_sets_per_split = []
            for left_run_indices, right_run_indices in selected_run_splits:
                mf_sets = []
                for run_indices in (left_run_indices, right_run_indices):
                    half_dataframes = [cut_dataframes[index] for index in run_indices]
                    half_run_seeds = [run_seeds[index] for index in run_indices]
                    aggregated_edges = combine_and_aggregate_runs(half_dataframes, runs_aggregation_method) # E2
                    selected_mfs = select_top_mfs(aggregated_edges, n_mfs, configuration["mf_score_method"], half_run_seeds) # E3
                    mf_sets.append(set(selected_mfs["MF_ensembl"]))
                mf_sets_per_split.append(mf_sets)
            mean_kuncheva_stability, enough_mfs = kuncheva_stability(mf_sets_per_split, n_genes_in_universe, n_mfs)
            stability_rows.append({"cohort": cohort, "configuration_name": configuration["configuration_name"], "runs_aggregation_method": configuration["runs_aggregation_method"], "alpha": configuration["alpha"], "beta": configuration["beta"], "mf_score_method": configuration["mf_score_method"], "k_edges": str(k_edges), "n_mfs": n_mfs, "enough_mfs": enough_mfs, "mean_kuncheva_stability": mean_kuncheva_stability})
    cohort_stability_results = pd.DataFrame(stability_rows)
    return cohort_stability_results


def main():
    parser = ArgumentParser()
    parser.add_argument("--edges-per-run-files", nargs="+", required=True)
    parser.add_argument("--run-seeds", nargs="+", required=True)
    parser.add_argument("--gene-names", required=True)
    parser.add_argument("--pairs-after-filters", required=True, help="Retained pairs after the two filters.")
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-mfs", type=int, required=True, help="Maximal number of selected MFs")
    parser.add_argument("--k-edges", nargs="+", type=int, required=True, help="Edge cut K.")
    parser.add_argument("--n-splits", type=int, required=True, help="Number of used distinct split-half pairs.")
    parser.add_argument("--split-seed", type=int, required=True)
    args = parser.parse_args()

    pairs_after_filters = pd.read_csv(args.pairs_after_filters, sep="\t", usecols=["MF_ensembl", "CpG"])
    pair_ids_after_filters = set(zip(pairs_after_filters["MF_ensembl"], pairs_after_filters["CpG"]))
    list_of_filtered_dataframes = filter_runs_dataframes(args.edges_per_run_files, pair_ids_after_filters=pair_ids_after_filters)
    gene_names = pd.read_csv(args.gene_names, sep="\t", header=None)[0]
    runs_aggregation_methods = create_runs_aggregation_methods({"runs_aggregation_method": ["mean", "borda", "freq"], "alpha": [0.0, 1.0, 2.0], "beta": [0.0, 0.5, 1.0]})
    configurations = create_configuration_grid(runs_aggregation_methods) # runs aggregation methods plus mf score methods and configuration name
    run_splits = complementary_splits(len(list_of_filtered_dataframes))
    selected_run_splits = sample_splits(run_splits, args.n_splits, args.split_seed)
    cohort_stability_results = calculate_stability_for_all_settings(list_of_filtered_dataframes, args.run_seeds, len(gene_names), args.cohort, args.n_mfs, args.k_edges, configurations, selected_run_splits)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    cohort_stability_results.to_csv(args.output, sep="\t", index=False, encoding="utf-8", float_format="%.6f")


if __name__ == "__main__":
    main()
