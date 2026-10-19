from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
from data_processing.data_preparation import prepare_data


def main():
    parser = ArgumentParser()
    parser.add_argument("--cpg-methylation-file", required=True, help="CpG methylation .tsv file (features × samples).")
    parser.add_argument("--gene-expression-file", required=True, help="Gene expression .tsv file (features × samples).")
    parser.add_argument("--raw-counts", required=True)
    parser.add_argument("--cpg-annotations-file", required=True) 
    parser.add_argument("--shared-samples", required=True)
    parser.add_argument("--coding-gene-list-file", default=None, help="Coding gene list.")
    parser.add_argument("--exclude-patients-of-type", default=None, help="For SKCM. The sample type of the other half when a cohort is split. Patients contributing both are removed from both halves. Only SKCM needs it.")
    parser.add_argument("--primary-tumour-type", required=True, help="The primary tumour type of a cohort (01: for solid tumors, 03: for LAML).")
    parser.add_argument("--out-feature-matrix", required=True)
    parser.add_argument("--out-gene-names", required=True)
    parser.add_argument("--out-cpg-names", required=True)
    parser.add_argument("--out-cpg-gene-region-map", required=True)
    parser.add_argument("--out-samples", required=True)
    parser.add_argument("--out-counts", required=True, help="Rounded integer counts.")
    parser.add_argument("--only-coding-genes", required=True)
    parser.add_argument("--n-hvgs", type=int, required=True)
    parser.add_argument("--sample-max-missing-share", type=float, required=True, help="Drop a sample if more than sample-max-missing-share of its CpGs or genes are missing.")
    parser.add_argument("--should-filter-cpgs", required=True)
    parser.add_argument("--should-aggregate-cpgs", required=True)
    parser.add_argument("--cpg-filter-genes-file", default=None, help="Gene list for the CpG filter.")
    parser.add_argument("--cpg-probe-coordinates-file", default=None, help="Hg38 CpG probe coordinates.")
    parser.add_argument("--gene-regions-gtf-file", default=None, help="Hg38 Gencode GTF.")
    parser.add_argument("--cpg-keep-promoters", required=True)
    parser.add_argument("--cpg-keep-bodies", required=True)
    parser.add_argument("--promoter-length", type=int, required=True)
    parser.add_argument("--promoter-overlap", type=int, required=True)
    args = parser.parse_args()

    feature_matrix, gene_names, cpg_names, cpg_gene_mapping, counts = prepare_data(cpg_methylation_file = args.cpg_methylation_file, gene_expression_file = args.gene_expression_file, raw_counts_file = args.raw_counts,
        cpg_annotations_file = args.cpg_annotations_file, only_coding_genes = args.only_coding_genes.lower() == "true",  coding_gene_list_file = args.coding_gene_list_file,
        shared_samples_file = args.shared_samples, n_hvgs = args.n_hvgs, primary_tumour_type = args.primary_tumour_type, sample_max_missing_share = args.sample_max_missing_share,
        exclude_patients_of_type = args.exclude_patients_of_type, should_filter_cpgs = args.should_filter_cpgs.lower() == "true", should_aggregate_cpgs = args.should_aggregate_cpgs.lower() == "true", cpg_filter_genes_file = args.cpg_filter_genes_file,
        cpg_probe_coordinates_file = args.cpg_probe_coordinates_file, gene_regions_gtf_file = args.gene_regions_gtf_file, cpg_keep_promoters = args.cpg_keep_promoters.lower() == "true",
        cpg_keep_bodies = args.cpg_keep_bodies.lower() == "true", promoter_length = args.promoter_length, promoter_overlap = args.promoter_overlap) 

    for output_path in (args.out_feature_matrix, args.out_gene_names, args.out_cpg_names, args.out_cpg_gene_region_map, args.out_samples, args.out_counts):
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    feature_matrix.to_csv(args.out_feature_matrix, sep="\t")
    pd.Series(gene_names).to_csv(args.out_gene_names, index=False, header=False)
    pd.Series(cpg_names).to_csv(args.out_cpg_names, index=False, header=False)
    cpg_gene_mapping.to_csv(args.out_cpg_gene_region_map, sep="\t", index=False)
    pd.Series(feature_matrix.index, name="sample_id").to_csv(args.out_samples, sep="\t", index=False)
    counts.to_csv(args.out_counts, sep="\t")


if __name__ == "__main__":
    main()
