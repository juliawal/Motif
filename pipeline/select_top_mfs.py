from pathlib import Path
from argparse import ArgumentParser
import pandas as pd
from selection.top_mfs_selection import select_top_mfs
from pipeline.apply_selection import load_selected_setting


def main():
    parser = ArgumentParser()
    parser.add_argument("--aggregated-edges", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--selected-setting-file", required=True)
    parser.add_argument("--n-mfs", type=int, required=True)
    parser.add_argument("--run-seeds", nargs="+", required=True)
    args = parser.parse_args()

    selected_setting_dict = load_selected_setting(args.selected_setting_file)
    aggregated_edges = pd.read_csv(args.aggregated_edges, sep="\t")

    # E3
    selected_mfs = select_top_mfs(aggregated_edges, args.n_mfs, selected_setting_dict["mf_score_method"], args.run_seeds)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    selected_mfs.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
