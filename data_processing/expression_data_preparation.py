import re
import numpy as np
import pandas as pd


def strip_ensembl_version(ids):
    ensembl_version = re.compile(r"\.\d+(?=_PAR_Y$|$)")
    return pd.Index(ids).str.replace(ensembl_version, "", regex=True)


def get_gene_information(gene_annotation_file, ensembl_ids): # for F0 and F3
    gene_information = pd.read_csv(gene_annotation_file, sep="\t", dtype=str)
    gene_information = gene_information.rename(columns={"selected_gene_symbol": "symbol"})
    gene_information = gene_information[["ensembl_id", "symbol", "entrez_id", "chromosome"]]
    gene_information["chromosome"] = gene_information["chromosome"].str.removeprefix("chr")
    gene_ids =  strip_ensembl_version([str(i) for i in ensembl_ids])
    gene_ids = list(dict.fromkeys(gene_ids))
    gene_information = gene_information.set_index("ensembl_id")
    gene_information = gene_information.reindex(gene_ids).reset_index()
    gene_information = gene_information.rename(columns={"index": "ensembl_id"})
    return gene_information


def load_symbol_aliases(gene_annotation_file): # for F1
    gene_information = pd.read_csv(gene_annotation_file, sep="\t", dtype=str).rename(columns={"selected_gene_symbol": "symbol", "previous_gene_symbols": "prev_symbols", "gene_symbol_aliases": "alias_symbols"})
    aliases = {}
    for symbol, prev, alias in gene_information[["symbol", "prev_symbols", "alias_symbols"]].itertuples(index=False):
        if not isinstance(symbol, str):
            continue
        names = []
        for field in (prev, alias):
            if isinstance(field, str) and field:
                names.extend(part for part in field.split("|") if part)
        seen = dict.fromkeys(n for n in names if n != symbol)
        if seen:
            aliases[symbol] = list(seen)
    return aliases


def load_ensembl_to_entrez_mapping(gene_annotation_file):
    ensembl_to_entrez_mapping = pd.read_csv(gene_annotation_file, sep="\t")
    ensembl_to_entrez_mapping = ensembl_to_entrez_mapping.dropna(subset=["entrez_id"])
    ensembl_to_entrez_mapping["entrez_id"] = ensembl_to_entrez_mapping["entrez_id"].astype(int)
    ensembl_to_entrez_dict = dict(zip(ensembl_to_entrez_mapping["ensembl_id"], ensembl_to_entrez_mapping["entrez_id"]))
    return ensembl_to_entrez_dict


def use_only_coding_genes(gene_expression, only_coding_genes, coding_gene_list_file):
    if not only_coding_genes:
        return gene_expression
    with open(coding_gene_list_file) as file:
        gene_ids = [line.strip() for line in file if line.strip()]
    coding_gene_ids = list(dict.fromkeys(gene_ids))
    coding_genes = gene_expression.loc[:, gene_expression.columns.isin(coding_gene_ids)] 
    return coding_genes

def select_hvgs(gene_expression, n_hvgs):
    """Select top n_hvgs genes by their median absolute deviation (MAD)."""
    values = gene_expression.to_numpy()
    gene_medians = np.median(values, axis=0)
    mads = np.median(np.abs(values - gene_medians), axis=0)
    genes_ranked_by_mad = pd.DataFrame({"original_column_number": np.arange(gene_expression.shape[1]), "gene_id": gene_expression.columns.astype(str), "mad": mads})
    genes_ranked_by_mad = genes_ranked_by_mad.sort_values(["mad", "gene_id"], ascending=[False, True])
    retained = genes_ranked_by_mad.iloc[:n_hvgs]
    selected_gene_column_indices = np.sort(retained["original_column_number"].to_numpy())
    gene_expression =  gene_expression.iloc[:, selected_gene_column_indices]
    return gene_expression

def round_counts(values):
    rounded = np.rint(values)
    return rounded


def drop_incomplete_genes(gene_expression):
    return gene_expression.loc[:, np.isfinite(gene_expression.to_numpy()).all(axis=0)]
 

def intersect_gene_frames(gene_expression, raw_counts):
    genes_also_in_raw_counts = gene_expression.columns.isin(raw_counts.columns)
    common_genes = gene_expression.columns[genes_also_in_raw_counts]
    return gene_expression.loc[:, common_genes]
