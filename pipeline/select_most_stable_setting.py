from argparse import ArgumentParser
from pathlib import Path
import pandas as pd


def load_and_concatenate_stability_tables(grid_stability_files):
    cohort_stability_tables = []
    for grid_stability_file in grid_stability_files:
        cohort_stability_results = pd.read_csv(grid_stability_file, sep="\t", dtype={"k_edges": str})
        cohort_stability_tables.append(cohort_stability_results)
    stability_results = pd.concat(cohort_stability_tables, ignore_index=True)
    return stability_results


def prepare_stability_comparison(stability_results):
    stability_results_with_enough_mfs = stability_results
    for _, cohort_results_for_one_setting in stability_results_with_enough_mfs.groupby(["k_edges", "configuration_name"]):
        if not all(cohort_results_for_one_setting["enough_mfs"]):
            stability_results_with_enough_mfs = stability_results_with_enough_mfs.drop(index=cohort_results_for_one_setting.index)
    stability_results_with_enough_mfs = stability_results_with_enough_mfs.reset_index(drop=True)
    stability_matrix = stability_results_with_enough_mfs.copy()
    stability_matrix["setting_id"] = stability_matrix["k_edges"] + "__" + stability_matrix["configuration_name"]
    stability_matrix = stability_matrix.pivot(index="cohort", columns="setting_id", values="mean_kuncheva_stability") # create matrix with stability values as values
    return stability_results_with_enough_mfs, stability_matrix


def rank_stability_settings_by_median_kuncheva_stability(stability_results_with_enough_mfs, stability_matrix):
    ranking_rows = []
    for setting_id in stability_matrix.columns:
        k_edges, configuration_name = setting_id.split("__", 1)
        cohort_results_for_one_setting = stability_results_with_enough_mfs[(stability_results_with_enough_mfs["k_edges"] == k_edges) & (stability_results_with_enough_mfs["configuration_name"] == configuration_name)]
        first_cohort_result = cohort_results_for_one_setting.iloc[0]
        runs_aggregation_method = first_cohort_result["runs_aggregation_method"]
        ranking_rows.append({"setting_id": setting_id, "k_edges": k_edges, "configuration_name": configuration_name, "runs_aggregation_method": runs_aggregation_method, "alpha": first_cohort_result["alpha"], "beta": first_cohort_result["beta"], "mf_score_method": first_cohort_result["mf_score_method"], "kuncheva_median": stability_matrix[setting_id].median()})
    ranked_settings_with_stability = pd.DataFrame(ranking_rows).sort_values(["kuncheva_median", "setting_id"], ascending=[False, True]).reset_index(drop=True)
    ranked_settings_with_stability.insert(0, "selection_rank", range(1, len(ranked_settings_with_stability) + 1))
    return ranked_settings_with_stability


def select_best_settings(stability_results_with_enough_mfs, ranked_settings_with_stability, n_selected_top_settings):
    best_settings = ranked_settings_with_stability.head(n_selected_top_settings).copy()
    best_settings["n_mfs"] = stability_results_with_enough_mfs["n_mfs"].iloc[0]
    cohorts_results_table_for_best_settings = []
    for _, setting in best_settings.iterrows():
        matching_setting = ((stability_results_with_enough_mfs["k_edges"] == setting["k_edges"]) & (stability_results_with_enough_mfs["configuration_name"] == setting["configuration_name"]))
        cohort_results_for_one_setting = stability_results_with_enough_mfs.loc[matching_setting].sort_values("cohort")
        cohort_results_for_one_setting = cohort_results_for_one_setting[["cohort", "k_edges", "configuration_name", "n_mfs", "mean_kuncheva_stability"]]
        cohort_results_for_one_setting["setting_id"] = setting["setting_id"]
        cohort_results_for_one_setting["selection_rank"] = setting["selection_rank"]
        cohorts_results_table_for_best_settings.append(cohort_results_for_one_setting)
    cohorts_results_for_best_settings = pd.concat(cohorts_results_table_for_best_settings, ignore_index=True) # cohorts stability values for selected best settings
    return best_settings, cohorts_results_for_best_settings


def main():
    parser = ArgumentParser()
    parser.add_argument("--grid-stability-files", nargs="+", required=True, help="grid_stability.tsv file per cohort.")
    parser.add_argument("--n-selected-top-settings", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    stability_results = load_and_concatenate_stability_tables(args.grid_stability_files)
    stability_results_with_enough_mfs, stability_matrix = prepare_stability_comparison(stability_results)
    ranked_settings_with_stability = rank_stability_settings_by_median_kuncheva_stability(stability_results_with_enough_mfs, stability_matrix)
    best_settings, cohorts_results_for_best_settings = select_best_settings(stability_results_with_enough_mfs, ranked_settings_with_stability, args.n_selected_top_settings)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for output_name, output_table in (("stability_results_with_enough_mfs", stability_results_with_enough_mfs), ("ranking", ranked_settings_with_stability), ("selection", best_settings), ("selected_stability", cohorts_results_for_best_settings)):
        output_table.to_csv(output_dir / f"{output_name}.tsv", sep="\t", index=False)
    stability_matrix.to_csv(output_dir / "kuncheva_by_cohort.tsv", sep="\t")


if __name__ == "__main__":
    main()
