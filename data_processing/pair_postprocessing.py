"""Prepare the selected MF-CpG pairs for evaluation."""

from pathlib import Path
import pandas as pd
from data_processing.expression_data_preparation import get_gene_information


def prepare_evaluation_pairs(aggregated_edges_file, pairs_path, selected_mfs_file, gene_annotation_file, mapping_output): 
    Path(pairs_path).parent.mkdir(parents=True, exist_ok=True)
    Path(mapping_output).parent.mkdir(parents=True, exist_ok=True)

    selected_mf_ids = pd.read_csv(selected_mfs_file, sep="\t")["MF_ensembl"]
    aggregated_edges = pd.read_csv(aggregated_edges_file, sep="\t") 
    edges_to_keep = aggregated_edges["MF_ensembl"].isin(selected_mf_ids)
    edges_for_selected_mfs = aggregated_edges.loc[edges_to_keep].copy()
    edge_counts = edges_for_selected_mfs["MF_ensembl"].value_counts()

    gene_mapping = get_gene_information(gene_annotation_file, edges_for_selected_mfs["MF_ensembl"].unique())
    gene_mapping["sex_chromosome"] = gene_mapping["chromosome"].isin(["X", "Y"])
    gene_mapping["n_edges"] = gene_mapping["ensembl_id"].map(edge_counts)
    gene_mapping["mapped"] = gene_mapping["symbol"].notna()
    gene_mapping["entrez_available"]  = gene_mapping["entrez_id"].notna() 

    # Attach symbols and chromosomes to edges_for_selected_mfs.
    genes_by_ensembl = gene_mapping.set_index("ensembl_id")
    gene_chromosomes = genes_by_ensembl["chromosome"]
    edges_for_selected_mfs["MF"] = edges_for_selected_mfs["MF_ensembl"].map(genes_by_ensembl["symbol"])
    edges_for_selected_mfs["MF_chromosome"] = edges_for_selected_mfs["MF_ensembl"].map(gene_chromosomes) 
    mf_cpg_pairs_with_weight = edges_for_selected_mfs.dropna(subset=["MF"]) 
    mf_cpg_pairs_with_weight = mf_cpg_pairs_with_weight[["MF", "CpG", "MF_ensembl", "weight", "MF_chromosome"]]

    mf_cpg_pairs_with_weight.to_csv(pairs_path, sep="\t", index=False, encoding="utf-8")
    gene_mapping[["ensembl_id", "symbol", "entrez_id", "chromosome", "sex_chromosome", "n_edges", "mapped", "entrez_available"]].to_csv(mapping_output, sep="\t", index=False, encoding="utf-8")
    
    return mf_cpg_pairs_with_weight
