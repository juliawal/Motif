"""Compare target and background hit rates and test each binding track."""

import numpy as np
import pandas as pd
from scipy.stats import binomtest


def summarize_mf_hit_rates(target_cpgs_by_mf, background_cpgs_by_mf, hit_cpgs_by_mf, track):
    hit_rate = []
    for mf_symbol, target_cpgs in target_cpgs_by_mf.items():
        background_cpgs = background_cpgs_by_mf[mf_symbol]
        hit_cpgs = hit_cpgs_by_mf[mf_symbol]
        n_target_cpgs = len(target_cpgs)
        n_background_cpgs = len(background_cpgs)
        if not background_cpgs:
            hit_rate.append({"MF": mf_symbol, "track": track, "status": "no_background", "n_target_cpgs": n_target_cpgs, "n_background_cpgs": 0})
            continue
        n_target_hits = len(target_cpgs & hit_cpgs)
        n_background_hits = len(background_cpgs & hit_cpgs)
        hit_rate.append({"MF": mf_symbol, "track": track, "status": "tested", "n_target_hits": n_target_hits, "n_background_hits": n_background_hits, "n_target_cpgs": n_target_cpgs, "n_background_cpgs": n_background_cpgs, "target_hit_rate": n_target_hits / n_target_cpgs, "bg_hit_rate": n_background_hits / n_background_cpgs})
    return pd.DataFrame(hit_rate)


def calculate_hit_rate_sign_test(mf_hit_rates):
    """Calculate hit rate differences of the measured MFs."""
    hit_rate_differences = np.array([])
    if not mf_hit_rates.empty:
        measured_mfs = mf_hit_rates["status"] == "tested"
        measured_mf_hit_rates = mf_hit_rates.loc[measured_mfs]
        if not measured_mf_hit_rates.empty:
            hit_rate_differences = (measured_mf_hit_rates["target_hit_rate"] - measured_mf_hit_rates["bg_hit_rate"]).to_numpy()
    nonzero_differences = hit_rate_differences[hit_rate_differences != 0] # don't use zero differences
    n_measured_mfs = len(hit_rate_differences)
    n_nonzero_differences = len(nonzero_differences)
    n_positive_differences = (hit_rate_differences > 0).sum()
    sign_test_p_value = np.nan
    median_hit_rate_difference = np.nan
    if n_measured_mfs:
        median_hit_rate_difference = np.median(hit_rate_differences)
        sign_test_p_value = 1.0
        if n_nonzero_differences:
            sign_test = binomtest(n_positive_differences, n = n_nonzero_differences, p = 0.5, alternative = "greater")
            sign_test_p_value = sign_test.pvalue
    return {"n_mf": n_measured_mfs, "n_not_measured_mfs": len(mf_hit_rates) - n_measured_mfs, "n_nonzero_differences": n_nonzero_differences, "n_positive_differences": n_positive_differences, "median_difference": median_hit_rate_difference, "p_value_sign_test": sign_test_p_value}


def create_binding_site_evaluation_results(binding_site_evaluation_results, n_original_target_cpgs_by_mf, n_target_cpgs_excluded_by_mf):
    if binding_site_evaluation_results.empty:
        return pd.DataFrame(columns=["MF", "track", "status", "n_target_hits", "n_background_hits", "n_target_cpgs", "n_background_cpgs", "target_hit_rate", "bg_hit_rate", "jaspar_ids", "chipseq_peak_file", "n_target_input_cpgs", "n_target_cpgs_excluded"]) 
    binding_site_evaluation_results["n_target_input_cpgs"] = binding_site_evaluation_results["MF"].map(n_original_target_cpgs_by_mf)
    binding_site_evaluation_results["n_target_cpgs_excluded"] = binding_site_evaluation_results["MF"].map(n_target_cpgs_excluded_by_mf)
    return binding_site_evaluation_results
