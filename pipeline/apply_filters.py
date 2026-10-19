from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
from filtering.reverse_causality import reverse_causality_filter
from filtering.sex import sex_filter
from data_processing.expression_data_preparation import strip_ensembl_version
from data_processing.sample_preparation import read_shared_samples
from utils.methylation_utils import read_cpg_methylation_blocks


def gene_cpg_pairs_of_all_runs(edges_per_run_files):
    run_tables = []
    for edges_per_run_file in edges_per_run_files:
        run_pairs = pd.read_csv(edges_per_run_file, sep="\t", usecols=["MF_ensembl", "CpG"], dtype={"MF_ensembl": str, "CpG": str})
        run_tables.append(run_pairs)
    pairs = pd.concat(run_tables, ignore_index=True).drop_duplicates(ignore_index=True) # every (MF_ensembl, CpG) pair any run produced
    return pairs


def apply_sex_filter(args):
    pairs = gene_cpg_pairs_of_all_runs(args.edges_per_run_files)
    samples = read_shared_samples(args.samples)
    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    metadata = pd.read_csv(args.metadata, sep="\t", index_col=0)
    metadata = metadata.loc[samples]
    cpg_methylation_blocks = list(read_cpg_methylation_blocks(args.cpg_methylation, pairs["CpG"]))
    cpg_methylation = pd.concat(cpg_methylation_blocks)
    cpg_methylation = cpg_methylation.loc[:, samples]
    pairs_with_sex_flags, gene_results, cpg_results = sex_filter(pairs, counts, cpg_methylation, metadata, args.fdr_alpha, args.n_cpus, args.m_value_boundary_shift)
    pairs_after_sex_filter = pairs_with_sex_flags.loc[~pairs_with_sex_flags["both_sex_dependent"]]
    output = Path(args.output) 
    output.parent.mkdir(parents=True, exist_ok=True)
    pairs_after_sex_filter.to_csv(output, sep="\t", index=False)
    gene_columns = ["gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "contrast", "n_male", "n_female", "sex_dependent"]
    cpg_columns = ["cpg_id", "n_male", "n_female", "p_value", "padj", "sex_dependent"]
    if gene_results is None:
        gene_results = pd.DataFrame(columns=gene_columns)
    if cpg_results is None: 
        cpg_results = pd.DataFrame(columns=cpg_columns)
    gene_results.reindex(columns=gene_columns).to_csv(output.with_name("gene_sex_test.tsv"), sep="\t", index=False)
    cpg_results.reindex(columns=cpg_columns).to_csv(output.with_name("cpg_sex_test.tsv"), sep="\t", index=False)


def apply_reverse_causality_filter(args):
    pairs_after_sex_filter = pd.read_csv(args.pairs_after_sex_filter, sep="\t")
    samples = read_shared_samples(args.samples)
    gene_expression = pd.read_csv(args.gene_expression, sep="\t", index_col=0)
    gene_expression.index = strip_ensembl_version(gene_expression.index)
    gene_expression = gene_expression.loc[:, samples]
    cpg_methylation_blocks = (cpg_methylation_block.loc[:, samples] for cpg_methylation_block in read_cpg_methylation_blocks(args.cpg_methylation, pairs_after_sex_filter["CpG"]))
    gene_tss_annotations = pd.read_csv(args.gene_tss_annotations_file, sep="\t", index_col=0).rename(columns={"chromosome": "chrom", "gene_tss": "tss", "gene_strand": "strand"})
    cpg_coordinates = pd.read_csv(args.cpg_coordinates_file, sep="\t", index_col=0).rename_axis("cpg_id")
    pairs_with_reverse_causality_flags = reverse_causality_filter(pairs_after_sex_filter, gene_expression, cpg_methylation_blocks, cpg_coordinates, gene_tss_annotations, reverse_tss_window=args.reverse_tss_window)
    pairs_after_filters = pairs_with_reverse_causality_flags.loc[~pairs_with_reverse_causality_flags["reverse_causality_flag"]]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pairs_after_filters.to_csv(output, sep="\t", index=False)


def main():
    parser = ArgumentParser()
    steps = parser.add_subparsers(dest="step", required=True)

    sex_filter_parser = steps.add_parser("sex-filter")
    sex_filter_parser.add_argument("--cpg-methylation", required=True)
    sex_filter_parser.add_argument("--samples", required=True)
    sex_filter_parser.add_argument("--output", required=True)
    sex_filter_parser.add_argument("--edges-per-run-files", nargs="+", required=True)
    sex_filter_parser.add_argument("--counts", required=True)
    sex_filter_parser.add_argument("--metadata", required=True)
    sex_filter_parser.add_argument("--fdr-alpha", type=float, required=True)
    sex_filter_parser.add_argument("--m-value-boundary-shift", type=float, required=True)
    sex_filter_parser.add_argument("--n-cpus", type=int, required=True)

    reverse_causality_filter_parser = steps.add_parser("reverse-causality-filter")
    reverse_causality_filter_parser.add_argument("--cpg-methylation", required=True)
    reverse_causality_filter_parser.add_argument("--samples", required=True)
    reverse_causality_filter_parser.add_argument("--output", required=True)
    reverse_causality_filter_parser.add_argument("--pairs-after-sex-filter", required=True) 
    reverse_causality_filter_parser.add_argument("--gene-expression", required=True)
    reverse_causality_filter_parser.add_argument("--cpg-coordinates-file", required=True)
    reverse_causality_filter_parser.add_argument("--gene-tss-annotations-file", required=True)
    reverse_causality_filter_parser.add_argument("--reverse-tss-window", type=int, required=True)
    
    args = parser.parse_args()
    if args.step == "sex-filter":
        apply_sex_filter(args) 
    else:
        apply_reverse_causality_filter(args)


if __name__ == "__main__":
    main()
