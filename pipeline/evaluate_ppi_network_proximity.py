from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
from data_processing.expression_data_preparation import load_ensembl_to_entrez_mapping
from evaluation.network_proximity.gene_sets_preparation_based_on_ppi import get_background_genes, get_genes_degrees_in_ppi, get_mf_entrez_ids_in_ppi, build_degree_bins, get_degree_matched_genes, sample_degree_matched_genes
from evaluation.network_proximity.harmonic_centrality import harmonic_centrality_from_enzymes
from evaluation.network_proximity.hippie import load_hippie_network
from evaluation.network_proximity.pagerank import personalized_pagerank
from evaluation.network_proximity.p_value_calculation import calculate_degree_matched_p_value
import random

dnmt_tet_genes = {"DNMT1": 1786, "DNMT3A": 1788, "DNMT3B": 1789, "TET1": 80312, "TET2": 54790, "TET3": 200424}


def main():
    parser = ArgumentParser()
    parser.add_argument("--mf-cpg-pairs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--hippie-network-file", type=Path, required=True)
    parser.add_argument("--hippie-release", required=True)
    parser.add_argument("--gene-annotation-file", type=Path, required=True)
    parser.add_argument("--gene-ids", type=Path, required=True)
    parser.add_argument("--degree-bin-min-size", type=int, required=True)
    parser.add_argument("--n-random-sets", type=int, required=True)
    parser.add_argument("--pagerank-alpha", type=float, required=True)
    parser.add_argument("--random-seed", type=int, required=True)
    args = parser.parse_args()

    mf_ensembl_ids = pd.read_csv(args.mf_cpg_pairs, sep="\t", usecols=["MF_ensembl"])["MF_ensembl"].drop_duplicates().tolist()
    gene_ids = args.gene_ids.read_text().splitlines() # from Stage A
    hippie_network = load_hippie_network(args.hippie_network_file)
    ensembl_to_entrez_dict = load_ensembl_to_entrez_mapping(args.gene_annotation_file)
    proximity_background_genes = get_background_genes(hippie_network, gene_ids, ensembl_to_entrez_dict)
    genes_degrees_in_ppi = get_genes_degrees_in_ppi(hippie_network, proximity_background_genes)
    mf_entrez_ids = get_mf_entrez_ids_in_ppi(hippie_network, mf_ensembl_ids, ensembl_to_entrez_dict)

    ppi_network_proximity_results = pd.DataFrame(columns = ["measure", "n_random_sets", "p_value", "n_mfs_nominated", "n_mfs_in_network"])
    if mf_entrez_ids:
        degree_bins = build_degree_bins(genes_degrees_in_ppi, args.degree_bin_min_size)
        degree_matched_genes_by_mf = get_degree_matched_genes(mf_entrez_ids, genes_degrees_in_ppi, degree_bins, dnmt_tet_genes)
        proximity_random_generator = random.Random(args.random_seed)
        sampled_gene_sets = [] # draw all random gene sets once and use them for both measures
        for _ in range(args.n_random_sets):
            sampled_gene_set = sample_degree_matched_genes(degree_matched_genes_by_mf, mf_entrez_ids, proximity_random_generator)
            sampled_gene_sets.append(sampled_gene_set)

        for row, measure in enumerate(("pagerank", "harmonic")):
            if measure == "pagerank":
                proximity_scores_by_gene = personalized_pagerank(hippie_network, dnmt_tet_genes, args.pagerank_alpha)
            if measure == "harmonic":
                proximity_scores_by_gene = harmonic_centrality_from_enzymes(hippie_network, dnmt_tet_genes)
            p_value = calculate_degree_matched_p_value(proximity_scores_by_gene, mf_entrez_ids, sampled_gene_sets)
            ppi_network_proximity_results.loc[row] = [measure, args.n_random_sets, p_value, len(mf_ensembl_ids), len(mf_entrez_ids)]

    dnmt_tet_gene_symbols = ",".join(dnmt_tet_genes)
    dnmt_tet_entrez_ids = ",".join(str(entrez_id) for entrez_id in dnmt_tet_genes.values())

    ppi_network_proximity_results["network_release"] = args.hippie_release
    ppi_network_proximity_results["network_file"] = args.hippie_network_file.name
    ppi_network_proximity_results["mapping_file"] = args.gene_annotation_file.name
    ppi_network_proximity_results["background_gene_ids_file"] = args.gene_ids.name
    ppi_network_proximity_results["n_graph_nodes"] = hippie_network.number_of_nodes()
    ppi_network_proximity_results["n_graph_edges"] = hippie_network.number_of_edges()
    ppi_network_proximity_results["degree_bin_min_size"] = args.degree_bin_min_size
    ppi_network_proximity_results["random_seed"] = args.random_seed
    ppi_network_proximity_results["pagerank_alpha"] = args.pagerank_alpha
    ppi_network_proximity_results["dnmt_tet_gene_symbols"] = dnmt_tet_gene_symbols
    ppi_network_proximity_results["dnmt_tet_entrez_ids"] = dnmt_tet_entrez_ids

    args.output.parent.mkdir(parents=True, exist_ok=True) 
    ppi_network_proximity_results.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
