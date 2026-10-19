import numpy as np
import pandas as pd
from .expression_data_preparation import drop_incomplete_genes, intersect_gene_frames, round_counts, select_hvgs, strip_ensembl_version, use_only_coding_genes
from .cpg_filtering_and_aggregation import aggregate_cpgs_per_gene, build_cpg_gene_mapping, filter_cpgs
from .methylation_data_preparation import drop_incomplete_cpgs, drop_masked_cpgs, keep_cg_cpgs
from .sample_preparation import drop_incomplete_samples, read_shared_samples, select_tumour_samples


def read_feature_matrix(path):
    feature_matrix = pd.read_csv(path, sep="\t", index_col=0).transpose()
    feature_matrix.index.name = None 
    feature_matrix.columns.name =  None
    feature_matrix.columns = strip_ensembl_version(feature_matrix.columns) 
    return feature_matrix


def load_raw_counts(count_path):
    counts = read_feature_matrix(count_path)
    values = counts.to_numpy()
    return pd.DataFrame(round_counts(values).astype(int), index=counts.index.copy(), columns=counts.columns.copy()) 


def load_input_files(cpg_methylation_file, gene_expression_file, raw_counts_file, cpg_annotations_file):
    cpg_methylation = read_feature_matrix(cpg_methylation_file)
    gene_expression = read_feature_matrix(gene_expression_file)
    cpg_methylation = keep_cg_cpgs(cpg_methylation)
    cpg_methylation = drop_masked_cpgs(cpg_methylation, cpg_annotations_file)
    raw_counts =  load_raw_counts(raw_counts_file) 
    return cpg_methylation, gene_expression, raw_counts


def combine_expression_and_methylation(gene_expression, cpg_methylation):
    if not gene_expression.index.equals(cpg_methylation.index):
        raise ValueError("Expression and methylation frames have different sample indices.")
    feature_matrix = pd.concat([gene_expression, cpg_methylation], axis=1).astype(np.float32) # convert to float32 for less storage
    gene_names = list(gene_expression.columns)
    cpg_names = list(cpg_methylation.columns)
    return feature_matrix, gene_names, cpg_names

def align_matrices(gene_expression, cpg_methylation, samples):
    return gene_expression.loc[samples], cpg_methylation.loc[samples]


def prepare_data(cpg_methylation_file, gene_expression_file, raw_counts_file, cpg_annotations_file, should_filter_cpgs, should_aggregate_cpgs, cpg_filter_genes_file, cpg_probe_coordinates_file, gene_regions_gtf_file, cpg_keep_promoters, cpg_keep_bodies, only_coding_genes, promoter_length, promoter_overlap, n_hvgs, shared_samples_file, primary_tumour_type, sample_max_missing_share, coding_gene_list_file=None, exclude_patients_of_type=None):
    cpg_methylation, gene_expression, raw_counts = load_input_files(cpg_methylation_file, gene_expression_file, raw_counts_file, cpg_annotations_file)
  
    samples = pd.Index(read_shared_samples(shared_samples_file)) 
    samples = select_tumour_samples(samples, primary_tumour_type, exclude_patients_of_type)

    gene_expression, cpg_methylation = align_matrices(gene_expression, cpg_methylation, samples)
    gene_expression, cpg_methylation = drop_incomplete_samples(gene_expression, cpg_methylation, sample_max_missing_share)
    gene_expression = drop_incomplete_genes(gene_expression)
    cpg_methylation = drop_incomplete_cpgs(cpg_methylation) 
    gene_expression = use_only_coding_genes(gene_expression, only_coding_genes, coding_gene_list_file)
    gene_expression = intersect_gene_frames(gene_expression, raw_counts)
    gene_expression = select_hvgs(gene_expression, n_hvgs) 

    if should_aggregate_cpgs or should_filter_cpgs:
        cpg_gene_mapping = build_cpg_gene_mapping(cpg_methylation, should_filter_cpgs, cpg_filter_genes_file, cpg_probe_coordinates_file, gene_regions_gtf_file, keep_promoters=cpg_keep_promoters, keep_bodies=cpg_keep_bodies, promoter_length=promoter_length, promoter_overlap=promoter_overlap)
        if should_filter_cpgs:
            cpg_methylation = filter_cpgs(cpg_methylation, cpg_gene_mapping)
        if should_aggregate_cpgs:
            cpg_methylation = aggregate_cpgs_per_gene(cpg_methylation, cpg_gene_mapping)
    else:
        cpg_gene_mapping = pd.DataFrame(columns=["aggregated_cpg", "ensembl_gene_id", "cpg_id", "region_type", "n_cpgs"])
    
    feature_matrix, gene_names, cpg_names = combine_expression_and_methylation(gene_expression, cpg_methylation)
    counts = raw_counts.loc[feature_matrix.index]
    return feature_matrix, gene_names, cpg_names, cpg_gene_mapping, counts
