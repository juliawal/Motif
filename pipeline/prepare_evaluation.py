from argparse import ArgumentParser
from pathlib import Path
from data_processing.pair_postprocessing import prepare_evaluation_pairs


def main():
    parser = ArgumentParser()
    parser.add_argument("--aggregated-edges-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--selected-mfs-file", type=Path, required=True)
    parser.add_argument("--gene-annotation-file", type=Path, required=True)
    parser.add_argument("--output-mf-mapping", type=Path, required=True)
    args = parser.parse_args()

    prepare_evaluation_pairs(args.aggregated_edges_file, args.output, args.selected_mfs_file, args.gene_annotation_file, args.output_mf_mapping)


if __name__ == "__main__":
    main()
