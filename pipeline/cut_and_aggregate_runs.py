from pathlib import Path
from argparse import ArgumentParser
import pandas as pd
from selection.top_edges_selection import select_top_edges_per_seed
from utils.selection_utils import filter_runs_dataframes, random_generator_for_equal_scores
from selection.run_aggregation import combine_and_aggregate_runs
from pipeline.apply_selection import load_selected_setting


def main():
    parser = ArgumentParser()
    parser.add_argument("--edges-per-run-files", nargs="+", required=True)
    parser.add_argument("--run-seeds", nargs="+", required=True)
    parser.add_argument("--pairs-after-filters", required=True)
    parser.add_argument("--selected-setting-file", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    selected_setting_dict = load_selected_setting(args.selected_setting_file)
    pairs_after_filters = pd.read_csv(args.pairs_after_filters, sep="\t", usecols=["MF_ensembl", "CpG"])
    pair_ids_after_filters = set(zip(pairs_after_filters["MF_ensembl"], pairs_after_filters["CpG"]))
    list_of_filtered_dataframes = filter_runs_dataframes(args.edges_per_run_files, pair_ids_after_filters=pair_ids_after_filters)

    # E1
    k_edges = int(selected_setting_dict["k_edges"])
    cut_dataframes = []
    for dataframe, run_seed in zip(list_of_filtered_dataframes, args.run_seeds, strict=True):
        edge_random_generator = random_generator_for_equal_scores([run_seed])
        cut_dataframes.append(select_top_edges_per_seed(dataframe, k_edges, edge_random_generator))

    # E2
    runs_aggregation_method = {key: selected_setting_dict[key] for key in ("runs_aggregation_method", "alpha", "beta")}
    aggregated_edges = combine_and_aggregate_runs(cut_dataframes, runs_aggregation_method)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    aggregated_edges.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
