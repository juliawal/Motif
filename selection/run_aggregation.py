import pandas as pd
from utils.selection_utils import ordered_edges


def combine_and_aggregate_runs(cut_dataframes, runs_aggregation_method):
    run_weight_tables, weight_columns = [], []
    for run_index, dataframe in enumerate(cut_dataframes):
        run_weight_column = f"weight_{run_index}"
        run_weight_tables.append(dataframe.rename(columns={"weight": run_weight_column}).set_index(["MF_ensembl", "CpG"])[[run_weight_column]])
        weight_columns.append(run_weight_column)
    edges_with_runs_weights = pd.concat(run_weight_tables, axis=1)
    edges_with_runs_weights[weight_columns] = edges_with_runs_weights[weight_columns].fillna(0.0)
    return aggregate_runs(edges_with_runs_weights.reset_index(), weight_columns, runs_aggregation_method)


def mean_edge_weight_over_runs(edges_with_runs_weights, weight_columns):
    edges_with_runs_weights["weight"] = edges_with_runs_weights[weight_columns].mean(axis=1)
    return ordered_edges(edges_with_runs_weights)


def freq_edge_weight_over_runs(edges_with_runs_weights, weight_columns, alpha, beta):
    edge_present_per_run = edges_with_runs_weights[weight_columns] > 0
    n_runs_with_edge = edge_present_per_run.sum(axis=1)
    mean_present_edge_weight = edges_with_runs_weights[weight_columns].sum(axis=1) / n_runs_with_edge
    max_edge_weight = edges_with_runs_weights[weight_columns].max(axis=1)
    edges_with_runs_weights["weight"] = (beta * max_edge_weight + (1 - beta) * mean_present_edge_weight) * (n_runs_with_edge / len(weight_columns)).pow(alpha)
    return ordered_edges(edges_with_runs_weights)


def borda_edge_weight_over_runs(edges_with_runs_weights, weight_columns):
    borda_points = pd.DataFrame(0.0, index=edges_with_runs_weights.index, columns=weight_columns)
    for run_weight_column in weight_columns:
        edge_present_per_run = edges_with_runs_weights[run_weight_column] > 0
        present_weights = edges_with_runs_weights[run_weight_column][edge_present_per_run]
        n_present = len(present_weights)
        if n_present:
            ranks = present_weights.rank(method="average", ascending=False)
            rank_points = (n_present - ranks + 1) / n_present
            borda_points.loc[edge_present_per_run, run_weight_column] = rank_points
    edges_with_runs_weights["weight"] = borda_points.mean(axis=1)
    return ordered_edges(edges_with_runs_weights)


def aggregate_runs(edges_with_runs_weights, weight_columns, runs_aggregation_method):
    edges_with_runs_weights = edges_with_runs_weights.copy()
    aggregated_edges = edges_with_runs_weights
    method_name = runs_aggregation_method["runs_aggregation_method"]
    if method_name == "mean":
        aggregated_edges = mean_edge_weight_over_runs(edges_with_runs_weights, weight_columns)
    if method_name == "borda":
        aggregated_edges = borda_edge_weight_over_runs(edges_with_runs_weights, weight_columns)
    if method_name == "freq":
        aggregated_edges = freq_edge_weight_over_runs(edges_with_runs_weights, weight_columns, runs_aggregation_method["alpha"], runs_aggregation_method["beta"])
    return aggregated_edges
