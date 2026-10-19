def get_background_genes(hippie_network, gene_ids, ensembl_to_entrez_dict):
    background_genes_in_hippie_ppi = set()
    for gene_ensembl_id in gene_ids:
        gene_entrez_id = ensembl_to_entrez_dict.get(gene_ensembl_id)
        if gene_entrez_id in hippie_network:
            background_genes_in_hippie_ppi.add(gene_entrez_id)
    return background_genes_in_hippie_ppi

    
def get_genes_degrees_in_ppi(hippie_network, background_gene_ids):
    genes_degrees_in_ppi = {}
    for gene_entrez_id, gene_degree in hippie_network.degree():
        if gene_entrez_id in background_gene_ids:
            genes_degrees_in_ppi[gene_entrez_id] = gene_degree
    return genes_degrees_in_ppi


def get_mf_entrez_ids_in_ppi(hippie_network, mf_ensembl_ids, ensembl_to_entrez_dict):
    mapped_mf_entrez_ids = []
    for mf_ensembl_id in mf_ensembl_ids:
        mf_entrez_id = ensembl_to_entrez_dict.get(mf_ensembl_id)
        if mf_entrez_id in hippie_network:
            mapped_mf_entrez_ids.append(mf_entrez_id)
    mf_entrez_ids = list(dict.fromkeys(mapped_mf_entrez_ids)) # different Ensembl IDs can map to the same Entrez ID; >> remove duplicate Entrez IDs
    return mf_entrez_ids


def build_degree_bins(genes_degrees_in_ppi, degree_bin_min_size):
    if degree_bin_min_size > len(genes_degrees_in_ppi):
        raise ValueError(f"degree_bin_min_size={degree_bin_min_size} exceeds the {len(genes_degrees_in_ppi)} background genes")
    gene_ids_by_degree = {}
    for gene_entrez_id, gene_degree in genes_degrees_in_ppi.items():
        if gene_degree not in gene_ids_by_degree:
            gene_ids_by_degree[gene_degree] = [] 
        gene_ids_by_degree[gene_degree].append(gene_entrez_id)

    degree_bins = {}
    for gene_degree in sorted(gene_ids_by_degree):
        gene_ids_in_bin = []
        for gene_entrez_id in gene_ids_by_degree[gene_degree]:
            gene_ids_in_bin.append(gene_entrez_id) # start with all genes that have exactly this degree
        for degree_difference in range(1, max(gene_ids_by_degree) + 1):
            if len(gene_ids_in_bin) >= degree_bin_min_size:
                break
            for degree in (gene_degree - degree_difference, gene_degree + degree_difference):
                for gene_entrez_id in gene_ids_by_degree.get(degree, []):
                    gene_ids_in_bin.append(gene_entrez_id)
        degree_bins[gene_degree] = gene_ids_in_bin
    return degree_bins


def get_degree_matched_genes(mf_entrez_ids, genes_degrees_in_ppi, degree_bins, dnmt_tet_genes):
    enzyme_entrez_ids = set(dnmt_tet_genes.values())
    degree_matched_genes_by_mf = {}
    for mf_entrez_id in mf_entrez_ids:
        mf_degree = genes_degrees_in_ppi[mf_entrez_id]
        gene_ids_in_bin = degree_bins[mf_degree]
        # Keep all genes from this bin, except the MF itself and the enzymes.
        same_degree_genes = []
        for gene_entrez_id in gene_ids_in_bin:
            if gene_entrez_id == mf_entrez_id or gene_entrez_id in enzyme_entrez_ids:
                continue
            same_degree_genes.append(gene_entrez_id)
        degree_matched_genes_by_mf[mf_entrez_id] = same_degree_genes
    return degree_matched_genes_by_mf

 
def sample_degree_matched_genes(degree_matched_genes_by_mf, mfs, proximity_random_generator):
    sampled_genes = set()
    for mf in mfs:
        degree_matched_genes = degree_matched_genes_by_mf[mf]
        sampled_gene_selected = False
        while not sampled_gene_selected:
            sampled_gene = proximity_random_generator.choice(degree_matched_genes)
            if sampled_gene not in sampled_genes:
                sampled_genes.add(sampled_gene) 
                sampled_gene_selected = True
    return list(sampled_genes)
