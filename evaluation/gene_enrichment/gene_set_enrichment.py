import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

def calculate_gene_set_enrichment(mfs_as_gene_symbols, background_genes_symbols, terms_by_library):
    enrichment_columns = ["library_name", "term_name", "n_background_genes", "n_term_genes", "n_mf_symbols", "n_mf_symbols_in_term", "p_value", "adj_p_value", "observed_to_expected_term_overlap_ratio", "mf_symbols_in_term", "n_tested_terms"]
    background_genes_symbols = set(background_genes_symbols)
    n_background_genes = len(background_genes_symbols)
    mfs_as_gene_symbols = set(mfs_as_gene_symbols)
    n_mf_symbols = len(mfs_as_gene_symbols)
    if n_mf_symbols == 0:
        return pd.DataFrame(columns = enrichment_columns)
    enrichment_rows = []
    for library_name, genes_by_term in terms_by_library.items():
        for term_name, gene_symbols in genes_by_term.items():
            term_symbols_in_background = gene_symbols & background_genes_symbols
            n_term_genes = len(term_symbols_in_background)
            if n_term_genes == 0:
                continue
            mf_symbols_in_term = mfs_as_gene_symbols & term_symbols_in_background
            n_mf_symbols_in_term = len(mf_symbols_in_term)
            p_value = hypergeom.sf(n_mf_symbols_in_term - 1, n_background_genes, n_term_genes, n_mf_symbols) # .sf() gives the probability of more than x hits; >> - 1 to get the probability of at least n_mf_symbols_in_term hits
            observed_to_expected_term_overlap_ratio = (n_mf_symbols_in_term / n_mf_symbols) / (n_term_genes / n_background_genes)
            enrichment_rows.append({"library_name": library_name, "term_name": term_name, "n_background_genes": n_background_genes, "n_term_genes": n_term_genes, "n_mf_symbols": n_mf_symbols, "n_mf_symbols_in_term": n_mf_symbols_in_term, "p_value": p_value, "observed_to_expected_term_overlap_ratio": observed_to_expected_term_overlap_ratio, "mf_symbols_in_term": ";".join(mf_symbols_in_term)})
    if not enrichment_rows:
        return pd.DataFrame(columns = enrichment_columns)
    gene_enrichment_results = pd.DataFrame(enrichment_rows, columns = enrichment_columns)
    gene_enrichment_results["adj_p_value"] = multipletests(gene_enrichment_results["p_value"], method="fdr_bh")[1]
    gene_enrichment_results["n_tested_terms"] = len(gene_enrichment_results)
    return gene_enrichment_results
