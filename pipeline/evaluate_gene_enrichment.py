from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
from evaluation.gene_enrichment.gene_set_enrichment import calculate_gene_set_enrichment
from data_processing.expression_data_preparation import get_gene_information


def read_gene_sets(gene_set_dir, gene_set_libraries):
    terms_by_library = {}
    for library_name in gene_set_libraries:
        gene_set_file = gene_set_dir / f"{library_name}.gmt"
        genes_by_term = {}
        with gene_set_file.open() as gene_set_library_file:
            for line in gene_set_library_file:
                term_name, _, *gene_symbols = line.rstrip("\n").split("\t") 
                genes_by_term[term_name] = {gene_symbol for gene_symbol in gene_symbols if gene_symbol}
        terms_by_library[library_name] = genes_by_term
    return terms_by_library

 
def main():
    parser = ArgumentParser()
    parser.add_argument("--mf-cpg-pairs", type=Path, required=True)
    parser.add_argument("--gene-annotation-file", type=Path, required=True)
    parser.add_argument("--gene-set-dir", type=Path, required=True)
    parser.add_argument("--gene-set-libraries", nargs="+", required=True)
    parser.add_argument("--background-gene-ids", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-overlap", type=int, required=True)
    args = parser.parse_args()

    mf_cpg_pairs_with_weight = pd.read_csv(args.mf_cpg_pairs, sep="\t")
    mfs_as_gene_symbols = mf_cpg_pairs_with_weight["MF"]
    background_gene_ids = args.background_gene_ids.read_text().splitlines() # ensembl IDs of all genes after Stage A
    gene_information = get_gene_information(args.gene_annotation_file, background_gene_ids)
    background_genes_symbols = gene_information["symbol"].dropna() 
    terms_by_library = read_gene_sets(args.gene_set_dir, args.gene_set_libraries)
    gene_enrichment_results = calculate_gene_set_enrichment(mfs_as_gene_symbols, background_genes_symbols, terms_by_library) 
    gene_enrichment_results["term_has_min_mf_overlap"] = gene_enrichment_results["n_mf_symbols_in_term"] >= args.min_overlap
    args.output.parent.mkdir(parents=True, exist_ok=True)
    gene_enrichment_results.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
