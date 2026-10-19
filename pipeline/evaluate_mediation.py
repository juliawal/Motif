"""Read one cohort, run the F4 mediation screen and save its results."""

from pathlib import Path
from argparse import ArgumentParser
import pandas as pd
from evaluation.mediation.mediation_triplet_selection import select_pairs_with_cpg_coordinates, select_triplets
from evaluation.mediation.age_and_sex_metadata_preparation import prepare_age_and_sex_metadata
from evaluation.mediation.mediation_regression import fit_mediation_triplets
from evaluation.mediation.joint_significance_and_summary import apply_joint_significance, summarize_mediation
from data_processing.expression_data_preparation import strip_ensembl_version
from data_processing.methylation_data_preparation import beta_to_m
from utils.methylation_utils import read_cpg_methylation_blocks


def load_cpg_m_values(cpg_methylation_file, cpgs_from_mf_cpg_pairs_with_coordinates, prepared_samples_with_age_and_sex, m_value_boundary_shift):
    cpg_methylation_blocks = list(read_cpg_methylation_blocks(cpg_methylation_file, cpgs_from_mf_cpg_pairs_with_coordinates))
    cpg_methylation = pd.concat(cpg_methylation_blocks)
    cpg_methylation = cpg_methylation.loc[:, prepared_samples_with_age_and_sex]
    cpg_methylation = cpg_methylation.loc[cpgs_from_mf_cpg_pairs_with_coordinates]
    cpg_m_values = beta_to_m(cpg_methylation, m_value_boundary_shift)
    return cpg_m_values


def main():
    parser = ArgumentParser()
    parser.add_argument("--mf-cpg-pairs", required=True)
    parser.add_argument("--prepared-samples", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--gene-expression-file", required=True)
    parser.add_argument("--cpg-methylation-file", required=True)
    parser.add_argument("--cpg-coordinates-file", required=True)
    parser.add_argument("--gene-annotation-file", required=True)
    parser.add_argument("--cohort", required=True)
    parser.add_argument("--output-triplets", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--max-cpg-tss-distance-bp", type=int, required=True)
    parser.add_argument("--min-samples", type=int, required=True)
    parser.add_argument("--fdr-alpha", type=float, required=True)
    parser.add_argument("--m-value-boundary-shift", type=float, required=True)
    args = parser.parse_args()

    # Establish one complete age/sex sample set before selecting or fitting triplets.
    mf_cpg_pairs_with_weight = pd.read_csv(args.mf_cpg_pairs, sep="\t")
    cpg_coordinates = pd.read_csv(args.cpg_coordinates_file, sep="\t", index_col=0).rename_axis("cpg_id")
    pairs_with_cpg_coordinates = select_pairs_with_cpg_coordinates(mf_cpg_pairs_with_weight, cpg_coordinates)
    samples = pd.read_csv(args.prepared_samples, sep="\t")
    sample_ids = pd.Index(samples["sample_id"].tolist())
    metadata = pd.read_csv(args.metadata, sep="\t", index_col=0)
    samples_with_age_and_sex, age_and_sex_metadata, age_and_sex_columns = prepare_age_and_sex_metadata(sample_ids, metadata)
    mediation_triplets_with_results = pd.DataFrame(columns=("triplet_id", "MF", "MF_ensembl", "CpG", "G_prime", "G_prime_ensembl", "G_prime_distance_bp", "n_samples", "triplet_status", "a_estimate", "a_pvalue", "b_estimate", "b_pvalue", "alpha_1_estimate", "gamma_age_estimate", "gamma_sex_estimate", "alpha_2_estimate", "delta_age_estimate", "delta_sex_estimate", "c_prime_estimate", "joint_adj_p_value", "mediation_supported"))

    if len(samples_with_age_and_sex) >= args.min_samples and not pairs_with_cpg_coordinates.empty:
        gene_expression = pd.read_csv(args.gene_expression_file, sep="\t", index_col=0)
        gene_expression.index = strip_ensembl_version(gene_expression.index)
        gene_expression = gene_expression.loc[:, samples_with_age_and_sex]
        cpgs_from_mf_cpg_pairs_with_coordinates = pd.Index(pairs_with_cpg_coordinates["CpG"].unique())
        cpg_m_values = load_cpg_m_values(args.cpg_methylation_file, cpgs_from_mf_cpg_pairs_with_coordinates, samples_with_age_and_sex, args.m_value_boundary_shift)
        gene_annotation = pd.read_csv(args.gene_annotation_file, sep="\t").rename(columns={"selected_gene_symbol": "symbol", "gene_tss": "tss"})
        mediation_triplets = select_triplets(pairs_with_cpg_coordinates, gene_annotation, gene_expression.index, max_cpg_tss_distance_bp = args.max_cpg_tss_distance_bp)
        mediation_triplets_with_results = fit_mediation_triplets(mediation_triplets, gene_expression, cpg_m_values, age_and_sex_metadata)
        mediation_triplets_with_results = apply_joint_significance(mediation_triplets_with_results, fdr_alpha=args.fdr_alpha)

    mediation_summary = summarize_mediation(mediation_triplets_with_results, args.cohort, len(mf_cpg_pairs_with_weight), len(samples_with_age_and_sex), age_and_sex_columns, args.min_samples, args.fdr_alpha, args.max_cpg_tss_distance_bp)
    Path(args.output_triplets).parent.mkdir(parents=True, exist_ok=True)
    mediation_triplets_with_results.to_csv(args.output_triplets, sep="\t", index=False) 
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    mediation_summary.to_csv(args.summary, sep="\t", index=False)


if __name__ == "__main__":
    main()
