from collections import defaultdict
import networkx as nx


def harmonic_centrality_from_enzymes(hippie_network, dnmt_tet_genes):
    harmonic_centrality_scores_to_dnmts_and_tets_by_gene = defaultdict(float)
    for enzyme_entrez_id in dnmt_tet_genes.values(): 
        shortest_paths_to_enzyme = nx.single_source_shortest_path_length(hippie_network, enzyme_entrez_id)
        for gene_entrez_id, distance_to_enzyme in shortest_paths_to_enzyme.items():
            if distance_to_enzyme > 0: # exclude enzyme itself with distance = 0
                harmonic_centrality_scores_to_dnmts_and_tets_by_gene[gene_entrez_id] += 1.0 / distance_to_enzyme
    return harmonic_centrality_scores_to_dnmts_and_tets_by_gene
