from argparse import ArgumentParser
from pathlib import Path
import pandas as pd


def main():
    parser = ArgumentParser()
    parser.add_argument("--cpg-map", type=Path, default=Path("resources/HM450.hg38.manifest.gencode.v36.probeMap"))
    parser.add_argument("--output-file", type=Path, default=Path("resources/cpg_annot.tsv"))
    args = parser.parse_args()

    cpg_probe_coordinates = pd.read_csv(args.cpg_map, sep="\t")
    cpg_probe_coordinates = cpg_probe_coordinates.rename(columns={"#id": "cpg_id"})
    cpg_probe_coordinates["pos"] = cpg_probe_coordinates["chromStart"] + 1 # convert 0 based CpG start to 1 based CpG position (for C2 and F4)
    cpg_coordinates = cpg_probe_coordinates[["cpg_id", "chrom", "pos"]]
    cpg_coordinates = cpg_coordinates.set_index("cpg_id")
    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    cpg_coordinates.to_csv(args.output_file, sep="\t")


if __name__ == "__main__":
    main()
