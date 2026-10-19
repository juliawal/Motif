from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
from pipeline.apply_selection import load_selected_setting
from evaluation.summary import evaluation_filenames, selection_columns, summarize_evaluation_results_per_cohort, summarize_cohort_testability, select_enrichment_top_terms


def collect_results_by_cohort(selected_setting_folder, cohorts, selected_setting_info, min_overlap, n_enrichment_top_terms):
    evaluation_results_tables = []
    enrichment_top_terms_tables = []
    for cohort in cohorts:
        cohort_folder = selected_setting_folder / cohort
        selected_mfs = pd.read_csv(cohort_folder / "selected_mfs" / "mfs.tsv", sep="\t", usecols=["MF_ensembl"])
        mf_cpg_pairs = pd.read_csv(cohort_folder / "evaluation_input" / "mf_cpg_pairs.tsv", sep="\t", usecols=["MF_ensembl", "MF"])
        evaluation_folder = cohort_folder / "evaluation"
        evaluation_results_by_stage = {}
        for stage, evaluation_filename in evaluation_filenames.items():
            evaluation_results_by_stage[stage] = pd.read_csv(evaluation_folder / evaluation_filename, sep="\t")
        evaluation_results_for_cohort = summarize_evaluation_results_per_cohort(cohort, selected_mfs, mf_cpg_pairs, evaluation_results_by_stage, selected_setting_info, min_overlap)
        evaluation_results_tables.append(evaluation_results_for_cohort)
        enrichment_top_terms_for_cohort = select_enrichment_top_terms(evaluation_results_by_stage["gene_set_enrichment"], cohort, selected_setting_info, min_overlap, n_enrichment_top_terms)
        enrichment_top_terms_tables.append(enrichment_top_terms_for_cohort)
    evaluation_results = pd.concat(evaluation_results_tables, ignore_index=True)
    enrichment_top_terms = pd.concat(enrichment_top_terms_tables, ignore_index=True)
    return evaluation_results, enrichment_top_terms
 


def load_selected_setting_for_summary(selected_setting_file):
    selected_setting_dict = load_selected_setting(selected_setting_file)
    selected_setting_info = {column: selected_setting_dict[column] for column in selection_columns}
    return selected_setting_info


def main():
    parser = ArgumentParser()
    parser.add_argument("--selected-setting-folder", type=Path, required=True)
    parser.add_argument("--cohorts", nargs="+", required=True)
    parser.add_argument("--min-overlap", type=int, required=True)
    parser.add_argument("--n-enrichment-top-terms", type=int, required=True)
    parser.add_argument("--selected-setting-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    selected_setting_info = load_selected_setting_for_summary(args.selected_setting_file)
    evaluation_results_by_cohort, enrichment_top_terms = collect_results_by_cohort(args.selected_setting_folder, args.cohorts, selected_setting_info, args.min_overlap, args.n_enrichment_top_terms)
    evaluation_summary = summarize_cohort_testability(evaluation_results_by_cohort)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    evaluation_results_by_cohort.to_csv(args.output_dir / "evaluation_by_cohort.tsv", sep="\t", index=False)
    evaluation_summary.to_csv(args.output_dir / "evaluation_summary.tsv", sep="\t", index=False)
    enrichment_top_terms.to_csv(args.output_dir / "enrichment_top_terms.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
