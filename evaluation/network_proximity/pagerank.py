import networkx as nx


def personalized_pagerank(hippie_network, dnmt_tet_genes, pagerank_alpha):
    personalization_by_enzyme = {}
    for enzyme_entrez_id in dnmt_tet_genes.values():
        personalization_by_enzyme[enzyme_entrez_id] = 1.0 / len(dnmt_tet_genes)
    pagerank_scores_by_gene = nx.pagerank(hippie_network, alpha=pagerank_alpha, personalization=personalization_by_enzyme)
    return pagerank_scores_by_gene