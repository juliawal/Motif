import networkx as nx
import pandas as pd


def load_hippie_network(hippie_network_file):
    ppi = pd.read_csv(hippie_network_file, sep="\t", header=None, usecols=[1, 3], names=["gene_entrez_id_1", "gene_entrez_id_2"])
    hippie_network = nx.Graph() # create unweighted networkx graph
    for gene_entrez_id_1, gene_entrez_id_2 in zip(ppi["gene_entrez_id_1"], ppi["gene_entrez_id_2"]):
        if gene_entrez_id_1 != gene_entrez_id_2: # dont't add self loops
            hippie_network.add_edge(gene_entrez_id_1, gene_entrez_id_2)
    return hippie_network
