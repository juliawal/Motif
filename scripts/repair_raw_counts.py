from argparse import ArgumentParser
from pathlib import Path
import numpy as np
import pandas as pd


def repair_raw_counts(raw_counts_file):
    """Transform log2(count + 1) values back to integer counts."""
    raw_counts = pd.read_csv(raw_counts_file, sep="\t", index_col=0)
    count_values = raw_counts.to_numpy() 
    reconstructed_counts = np.exp2(count_values) - 1.0  # back transform log2(count + 1)
    integer_counts = np.rint(reconstructed_counts).astype(np.int64) # round to integers
    pd.DataFrame(integer_counts, index=raw_counts.index, columns=raw_counts.columns).to_csv(raw_counts_file, sep="\t")


def main():
    parser = ArgumentParser() 
    parser.add_argument("raw_counts_files", type=Path, nargs="+")
    args = parser.parse_args()
    for raw_counts_file in args.raw_counts_files:
        repair_raw_counts(raw_counts_file)


if __name__ == "__main__":
    main()
