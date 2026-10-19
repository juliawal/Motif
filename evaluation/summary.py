import numpy as np
import pandas as pd


evaluation_filenames = {"binding_sites": "binding_sites.tsv", "ppi_proximity": "ppi_proximity.tsv", "gene_set_enrichment": "gene_set_enrichment.tsv", "mediation": "mediation.tsv"}
selection_columns = ("selection_rank",  # Setting rank by median Kuncheva stability across cohorts.
    "setting_id", "n_mfs")

cohort_result_columns = ("cohort",
    "evaluation_stage",  # Binding sites, PPI proximity, enrichment or mediation.
    "n_selected_mfs",  # Actual number of selected MFs from E3.
    "n_mf_symbols",  # Number of distinct MF gene symbols in F0 MF-CpG pairs.
    "mf_list_is_full",  # Retained MF count equals the configured n_mfs.
    "cohort_testable",  # Evaluation possible for this cohort; does not imply a significant result.

    "binding_site_fimo_n_mf",  # Number of MFs for which the difference (= target_hit_rate - bg_hit_rate) could be calculated.
    "binding_site_fimo_n_not_measured_mfs",  # Number of MFs in the FIMO input for which the hit-rate difference could not be calculated.
    "binding_site_fimo_n_nonzero_differences",  # Number of MFs with a hit-rate difference != 0.
    "binding_site_fimo_n_positive_differences",  # Number of MFs with a higher target than background hit rate.
    "binding_site_fimo_median_difference",  # Median of differences.
    "binding_site_fimo_p_value_sign_test",  # Sign-test p-value: are positive differences more common than negative differences?

    "binding_site_chipseq_n_mf",  # Number of MFs for which the difference (= target_hit_rate - bg_hit_rate) could be calculated.
    "binding_site_chipseq_n_not_measured_mfs",  # Number of MFs in the ChIP-seq input for which the hit-rate difference could not be calculated.
    "binding_site_chipseq_n_nonzero_differences",  # Number of MFs with a hit-rate difference != 0.
    "binding_site_chipseq_n_positive_differences",  # Number of MFs with a higher target than background hit rate.
    "binding_site_chipseq_median_difference",  # Median of valid differences.
    "binding_site_chipseq_p_value_sign_test",  # Sign-test p-value: are positive differences more common than negative differences?

    "ppi_n_mfs_nominated",  # Number of distinct F0 MF Ensembl IDs before mapping to PPI.
    "ppi_n_mfs_in_network",  # Number of distinct mapped MF Entrez nodes present in the PPI.
    "ppi_n_random_sets", "ppi_pagerank_p_value", "ppi_harmonic_p_value",

    "enrichment_min_overlap",  # Minimum distinct MF symbols in a term for reporting.
    "enrichment_n_tested_terms",  # Number of all tested terms across libraries.
    "enrichment_n_terms_with_min_mf_overlap",  # Number of terms with min_overlap, also if they are not significant.

    "mediation_n_samples",  # Number of shared samples with complete age and sex metadata.
    "mediation_n_input_pairs",  # Number of MF-CpG pairs before filtering for usable CpG coordinates.
    "mediation_n_candidate_triplets",  # Number of selected MF-CpG-G' combinations; also untestable triplets.
    "mediation_n_tested_triplets",  # Number of triplets with both paths testable.
    "mediation_n_supported_triplets",  # Number of triplets with BY adjusted signifikant p-value.
    "mediation_fdr_alpha", "mediation_max_cpg_tss_distance_bp") + selection_columns

summary_columns = ("stage", "n_testable_cohorts", "n_not_testable_cohorts",
    "n_total_cohorts",  # Number of all included cohorts, regardless of testability.
    "n_testable_cohorts_with_full_mf_list", "n_cohorts_with_full_mf_list", "n_testable_cohorts_without_full_mf_list", "n_cohorts_without_full_mf_list") + selection_columns

enrichment_top_term_columns = ("cohort",) + selection_columns + (
    "term_rank",  # Rank by adjusted p-value.
    "library_name", "term_name", "p_value",
    "adj_p_value",  # BH adjusted p-value across all tested terms and libraries.
    "observed_to_expected_term_overlap_ratio",  # Observed MF overlap divided by its random expectation.
    "n_mf_symbols_in_term", "mf_symbols_in_term")


def collect_binding_results(binding_site_evaluation_summary):
    binding_summary_by_track = binding_site_evaluation_summary.set_index("track")
    results = {}
    columns = ("n_mf", "n_not_measured_mfs", "n_nonzero_differences", "n_positive_differences", "median_difference", "p_value_sign_test")
    for track, column_prefix in (("FIMO", "binding_site_fimo"), ("ChIP_seq", "binding_site_chipseq")):
        for column in columns:
            results[f"{column_prefix}_{column}"] = binding_summary_by_track.at[track, column]
    binding_site_is_testable = np.isfinite(binding_site_evaluation_summary["p_value_sign_test"]).any()
    return binding_site_is_testable, results


def collect_proximity_results(ppi_network_proximity_results):
    results = {}
    for column in ("n_mfs_nominated", "n_mfs_in_network", "n_random_sets"):
        results[f"ppi_{column}"] = ppi_network_proximity_results[column].iloc[0]
    for measure_result in ppi_network_proximity_results.itertuples(index=False):
        results[f"ppi_{measure_result.measure}_p_value"] = measure_result.p_value
    return True, results


def collect_enrichment_results(gene_enrichment_results, min_overlap):
    term_has_min_mf_overlap = gene_enrichment_results["n_mf_symbols_in_term"] >= min_overlap
    results = {"enrichment_min_overlap": min_overlap, "enrichment_n_tested_terms": len(gene_enrichment_results), "enrichment_n_terms_with_min_mf_overlap": term_has_min_mf_overlap.sum()}
    enrichment_is_testable = not gene_enrichment_results.empty
    return enrichment_is_testable, results


def select_enrichment_top_terms(gene_enrichment_results, cohort, selected_setting_info, min_overlap, n_enrichment_top_terms):
    if n_enrichment_top_terms < 1:
        raise ValueError("n_enrichment_top_terms must be positive")
    term_has_min_mf_overlap = gene_enrichment_results["n_mf_symbols_in_term"] >= min_overlap
    terms_with_min_mf_overlap = gene_enrichment_results.loc[term_has_min_mf_overlap]
    terms_with_min_mf_overlap = terms_with_min_mf_overlap.sort_values(["adj_p_value", "p_value", "library_name", "term_name"])
    enrichment_top_terms = terms_with_min_mf_overlap.head(n_enrichment_top_terms).copy()
    enrichment_top_terms["cohort"] = cohort
    for column, value in selected_setting_info.items():
        enrichment_top_terms[column] = value
    enrichment_top_terms["term_rank"] = range(1, len(enrichment_top_terms) + 1)
    enrichment_top_terms = enrichment_top_terms.loc[:, list(enrichment_top_term_columns)]
    return enrichment_top_terms


def collect_mediation_results(mediation_summary):
    mediation_summary_row = mediation_summary.iloc[0]
    columns = ("n_samples", "n_input_pairs", "n_candidate_triplets", "n_tested_triplets", "n_supported_triplets", "max_cpg_tss_distance_bp", "fdr_alpha")
    results = {f"mediation_{column}": mediation_summary_row[column] for column in columns} 
    mediation_is_testable = mediation_summary_row["testable"]
    return mediation_is_testable, results 


def summarize_cohort_testability(evaluation_results_by_cohort):
    first_cohort_result = evaluation_results_by_cohort.iloc[0] # all cohorts use the same selected setting; >> read the setting from the first row
    selected_setting_info = first_cohort_result[list(selection_columns)].to_dict()
    summary_rows = []
    for stage in evaluation_filenames:
        cohort_is_in_stage = evaluation_results_by_cohort["evaluation_stage"] == stage
        cohort_results = evaluation_results_by_cohort.loc[cohort_is_in_stage]
        mf_list_is_full = cohort_results["mf_list_is_full"]
        cohorts_with_full_mf_list = cohort_results.loc[mf_list_is_full]
        cohorts_without_full_mf_list = cohort_results.loc[~mf_list_is_full]
        cohort_testable = cohort_results["cohort_testable"]
        stage_summary = {"stage": stage, "n_testable_cohorts": cohort_testable.sum(), "n_not_testable_cohorts": (~cohort_testable).sum(), "n_total_cohorts": len(cohort_results), "n_testable_cohorts_with_full_mf_list": cohorts_with_full_mf_list["cohort_testable"].sum(),
            "n_cohorts_with_full_mf_list": len(cohorts_with_full_mf_list), "n_testable_cohorts_without_full_mf_list": cohorts_without_full_mf_list["cohort_testable"].sum(), "n_cohorts_without_full_mf_list": len(cohorts_without_full_mf_list)}
        stage_summary.update(selected_setting_info)
        summary_rows.append(stage_summary)
    evaluation_summary = pd.DataFrame(summary_rows, columns=summary_columns)
    return evaluation_summary


def summarize_evaluation_results_per_cohort(cohort, selected_mfs, mf_cpg_pairs, evaluation_results_by_stage, selected_setting_info, min_overlap):
    n_selected_mfs = len(selected_mfs)
    n_mfs = selected_setting_info["n_mfs"]
    n_mf_symbols = mf_cpg_pairs["MF"].nunique() # different Ensembl IDs can share a symbol

    summary_results_by_stage = {"binding_sites": collect_binding_results(evaluation_results_by_stage["binding_sites"]),
        "ppi_proximity": collect_proximity_results(evaluation_results_by_stage["ppi_proximity"]),
        "gene_set_enrichment": collect_enrichment_results(evaluation_results_by_stage["gene_set_enrichment"], min_overlap),
        "mediation": collect_mediation_results(evaluation_results_by_stage["mediation"])}
    
    cohort_result_rows = []
    for evaluation_stage, (cohort_testable, results) in summary_results_by_stage.items():
        cohort_result_row = {"cohort": cohort, "evaluation_stage": evaluation_stage, "n_selected_mfs": n_selected_mfs, "n_mf_symbols": n_mf_symbols, "mf_list_is_full": n_selected_mfs == n_mfs, "cohort_testable": cohort_testable}
        cohort_result_row.update(selected_setting_info)
        cohort_result_row.update(results)
        cohort_result_rows.append(cohort_result_row)
    cohort_results = pd.DataFrame(cohort_result_rows, columns=cohort_result_columns)
    return cohort_results
