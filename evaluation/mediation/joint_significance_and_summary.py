import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def apply_joint_significance(mediation_triplets_with_results, fdr_alpha):
    triplet_has_fitted_models = mediation_triplets_with_results["triplet_status"] == "modelled" # include only triplets that were testable
    if triplet_has_fitted_models.any():
        modelled_triplets = mediation_triplets_with_results.loc[triplet_has_fitted_models]
        joint_p_values = modelled_triplets[["a_pvalue", "b_pvalue"]].max(axis = 1)
        mediation_supported, joint_adj_p_values, _, _ = multipletests(joint_p_values, alpha = fdr_alpha, method = "fdr_by") 
        mediation_triplets_with_results.loc[triplet_has_fitted_models, "joint_adj_p_value"] = joint_adj_p_values
        mediation_triplets_with_results.loc[triplet_has_fitted_models, "mediation_supported"] = mediation_supported 
        mediation_triplets_with_results.loc[triplet_has_fitted_models, "triplet_status"] = np.where(mediation_supported, "significant", "not_significant")
    return mediation_triplets_with_results


def summarize_mediation(mediation_triplets_with_results, cohort, n_mf_cpg_pairs, n_samples, age_and_sex_columns, min_samples, fdr_alpha, max_cpg_tss_distance_bp):
    triplet_was_tested = mediation_triplets_with_results["triplet_status"].isin(("not_significant", "significant"))
    n_tested_triplets = int(triplet_was_tested.sum())
    # A cohort is testable if it has the minimum number of samples and one or more triplets were tested.
    mediation_summary = pd.DataFrame([{"cohort": cohort, "n_samples": n_samples, "n_input_pairs": n_mf_cpg_pairs, "n_candidate_triplets": len(mediation_triplets_with_results), "n_tested_triplets": n_tested_triplets, "n_supported_triplets": mediation_triplets_with_results["mediation_supported"].sum(),
                                       "testable": n_samples >= min_samples and n_tested_triplets > 0, "fdr_alpha": fdr_alpha, "max_cpg_tss_distance_bp": max_cpg_tss_distance_bp, "age_and_sex_columns": ",".join(age_and_sex_columns)}])
    return mediation_summary
