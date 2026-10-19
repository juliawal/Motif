import numpy as np
import pandas as pd
import statsmodels.api as sm


def fit_mediation_models(mf_expression, cpg_m_values, g_prime_expression, age_and_sex_metadata):
    # Standardize all input matrices.
    standardized_values = [] # expression or methylation
    for values in (mf_expression, cpg_m_values, g_prime_expression):
        standard_deviation = values.std()
        standardized_values.append((values - values.mean()) / standard_deviation)
    mf_expression, cpg_m_values, g_prime_expression = standardized_values
    standardized_age_and_sex_values = []
    for age_or_sex_column in age_and_sex_metadata:
        age_or_sex_values = age_and_sex_metadata[age_or_sex_column].to_numpy()
        standard_deviation = age_or_sex_values.std() 
        standardized_age_or_sex_values = (age_or_sex_values - age_or_sex_values.mean()) / standard_deviation
        standardized_age_and_sex_values.append(standardized_age_or_sex_values)

    # X = standardized MF expression; M = standardized CpG methylation; Y = standardized G' expression; C = standardized age and, when both sexes are present, the standardized sex values; 1 = constant column to estimate a baseline value

    # Prepare inputs for first regression model: [1, age, sex (if included), X].
    cpg_methylation_model_inputs = sm.add_constant(np.column_stack([*standardized_age_and_sex_values, mf_expression]), has_constant="add")

    # cpg_methylation_model_inputs could be for example:
    #             Constant    Age    Sex    MF expression
    # Sample 1        1      -1.30    -1         -0.45
    # Sample 2        1      -0.45     1          1.30
    # Sample 3        1       0.45    -1         -1.30
    # Sample 4        1       1.30     1          0.45
    
    # Prepare inputs for second regression model: [1, age, sex (if included), X, M].
    g_prime_expression_model_inputs = sm.add_constant(np.column_stack([*standardized_age_and_sex_values, mf_expression, cpg_m_values]), has_constant="add")

    # First regression model: M ~ 1 + X + C.
    # M = alpha_1 + a * X + gamma_age * age + gamma_sex * sex + error_1
    cpg_model = sm.OLS(cpg_m_values, cpg_methylation_model_inputs).fit().get_robustcov_results(cov_type="HC3", use_t=True)

    # Second regression model: Y ~ 1 + X + M + C.
    # Y = alpha_2 + c_prime * X + b * M + delta_age * age + delta_sex * sex + error_2
    cpg_adjusted_gene_model = sm.OLS(g_prime_expression, g_prime_expression_model_inputs).fit().get_robustcov_results(cov_type="HC3", use_t=True)

    # HC3 standard errors allow unequal error variances across samples; two-sided t-tests use the residual degrees of freedom (<< Hayes and Cai 2007)

    # column order: constant column, age, optional sex, MF, CpG (second model only).
    a_estimate_index = len(standardized_age_and_sex_values) + 1
    b_estimate_index = a_estimate_index + 1
 
    mediation_model_results = {}
    mediation_model_results["a_estimate"] = cpg_model.params[a_estimate_index]
    mediation_model_results["a_pvalue"] = cpg_model.pvalues[a_estimate_index]
    mediation_model_results["b_estimate"] = cpg_adjusted_gene_model.params[b_estimate_index]
    mediation_model_results["b_pvalue"] = cpg_adjusted_gene_model.pvalues[b_estimate_index] 
    mediation_model_results["alpha_1_estimate"] = cpg_model.params[0]
    mediation_model_results["alpha_2_estimate"] = cpg_adjusted_gene_model.params[0] 
    mediation_model_results["c_prime_estimate"] = cpg_adjusted_gene_model.params[a_estimate_index]
    mediation_model_results["gamma_age_estimate"] = cpg_model.params[1]
    mediation_model_results["delta_age_estimate"] = cpg_adjusted_gene_model.params[1]
    if len(age_and_sex_metadata.columns) == 2:
        mediation_model_results["gamma_sex_estimate"] = cpg_model.params[2]
        mediation_model_results["delta_sex_estimate"] = cpg_adjusted_gene_model.params[2] 
    return mediation_model_results


def fit_mediation_triplets(triplets, gene_expression, cpg_m_values, age_and_sex_metadata):
    triplet_rows = []
    for triplet_number, triplet in enumerate(triplets.itertuples(index=False), start=1):
        mf_expression = gene_expression.loc[triplet.MF_ensembl].to_numpy()
        cpg_m_values_for_triplet = cpg_m_values.loc[triplet.CpG].to_numpy()
        g_prime_expression = gene_expression.loc[triplet.G_prime_ensembl].to_numpy()
        triplet_row = {"triplet_id": f"T{triplet_number}", "MF": triplet.MF, "MF_ensembl": triplet.MF_ensembl, "CpG": triplet.CpG, "G_prime": triplet.G_prime, "G_prime_ensembl": triplet.G_prime_ensembl, "G_prime_distance_bp": triplet.G_prime_distance_bp, "n_samples": len(age_and_sex_metadata)}
        for column in ("a_estimate", "a_pvalue", "b_estimate", "b_pvalue", "alpha_1_estimate", "gamma_age_estimate", "gamma_sex_estimate", "alpha_2_estimate", "delta_age_estimate", "delta_sex_estimate", "c_prime_estimate"):
            triplet_row[column] = pd.NA

        # Are all G' expression values finite?
        if not np.isfinite(g_prime_expression).all():
            triplet_row["triplet_status"] = "g_prime_expression_nonfinite"
            triplet_rows.append(triplet_row)
            continue

        # Are there constant G' expressions?
        if g_prime_expression.std() == 0:
            triplet_row["triplet_status"] = "g_prime_expression_constant"
            triplet_rows.append(triplet_row)
            continue

        # Do the model inputs have full rank?
        age_and_sex_values = age_and_sex_metadata.to_numpy(dtype=float) # convert age and True/False sex values; dtype=float converts True into 1.0 and False into 0.0
        g_prime_expression_model_inputs = np.column_stack([np.ones(len(mf_expression)), age_and_sex_values, mf_expression, cpg_m_values_for_triplet])
        if np.linalg.matrix_rank(g_prime_expression_model_inputs) != g_prime_expression_model_inputs.shape[1]:
            triplet_row["triplet_status"] = "model_inputs_not_full_rank"
            triplet_rows.append(triplet_row)
            continue

        mediation_model_results = fit_mediation_models(mf_expression, cpg_m_values_for_triplet, g_prime_expression, age_and_sex_metadata)

        # Are the calculated a/b estimates and their p-values finite?
        a_and_b_results = [mediation_model_results["a_estimate"], mediation_model_results["a_pvalue"], mediation_model_results["b_estimate"], mediation_model_results["b_pvalue"]]
        if not np.isfinite(a_and_b_results).all():
            triplet_row["triplet_status"] = "model_results_nonfinite"
            triplet_rows.append(triplet_row)
            continue

        triplet_row.update(mediation_model_results)
        triplet_row["triplet_status"] = "modelled"
        triplet_rows.append(triplet_row)

    mediation_triplets = pd.DataFrame.from_records(triplet_rows, columns=["triplet_id", "MF", "MF_ensembl", "CpG", "G_prime", "G_prime_ensembl", "G_prime_distance_bp", "n_samples", "triplet_status", "a_estimate", "a_pvalue", "b_estimate", "b_pvalue", "alpha_1_estimate", "gamma_age_estimate", "gamma_sex_estimate", "alpha_2_estimate", "delta_age_estimate", "delta_sex_estimate", "c_prime_estimate"])
    mediation_triplets["joint_adj_p_value"] = pd.NA
    mediation_triplets["mediation_supported"] = False
    return mediation_triplets
