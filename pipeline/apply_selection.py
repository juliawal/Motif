from argparse import ArgumentParser
from pathlib import Path
import pandas as pd


def load_selected_setting(selected_setting_file):
    selected_setting_df = pd.read_csv(selected_setting_file, sep="\t", dtype={"k_edges": str})
    selected_setting_dict = selected_setting_df.iloc[0].to_dict()
    return selected_setting_dict


def main():
    parser = ArgumentParser()
    parser.add_argument("--best-settings-file", required=True)
    parser.add_argument("--selection-rank", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    best_settings = pd.read_csv(args.best_settings_file, sep="\t", dtype={"k_edges": str})
    selected_setting_df = best_settings.iloc[[args.selection_rank - 1]]
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    selected_setting_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
