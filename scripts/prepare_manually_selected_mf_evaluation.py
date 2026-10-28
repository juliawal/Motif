from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
import yaml


def main():
    parser = ArgumentParser()
    parser.add_argument("--mfs", nargs="+", required=True)
    args = parser.parse_args()

    config = yaml.safe_load(Path("config.yaml").read_text())
    background_gene_ids = set()
    background_cpg_ids = set()
    mf_cpg_pair_tables = []

    for cohort in config["cohorts"]:
        prepared_folder = Path("results") / cohort / "prepared"
        background_gene_ids.update(pd.read_csv(prepared_folder / "regulator_names.tsv", header=None)[0])
        background_cpg_ids.update(pd.read_csv(prepared_folder / "target_names.tsv", header=None)[0])
        for rank in range(1, config["n_selected_stable_settings"] + 1):
            cohort_folder = Path(f"results/selected_settings/rank_{rank}") / cohort
            mf_cpg_pairs = pd.read_csv(cohort_folder / "evaluation_input/mf_cpg_pairs.tsv", sep="\t", dtype=str, usecols=["MF", "CpG", "MF_ensembl", "MF_chromosome"])
            pairs_to_keep = mf_cpg_pairs["MF"].isin(args.mfs) | mf_cpg_pairs["MF_ensembl"].isin(args.mfs)
            mf_cpg_pair_tables.append(mf_cpg_pairs.loc[pairs_to_keep])

    mf_cpg_pairs = pd.concat(mf_cpg_pair_tables).drop_duplicates(["MF_ensembl", "CpG"])
    selected_mfs = mf_cpg_pairs[["MF_ensembl", "MF"]].drop_duplicates("MF_ensembl")

    missing_mfs = set(args.mfs) - set(selected_mfs["MF"]) - set(selected_mfs["MF_ensembl"])
    if missing_mfs:
        parser.error("No MF-CpG pairs saved for: " + ", ".join(sorted(missing_mfs)))

    output_folder = Path("results/manual_mfs")
    for folder in ("selected_mfs", "evaluation_input", "prepared"):
        (output_folder / folder).mkdir(parents=True, exist_ok=True)
    selected_mfs.to_csv(output_folder / "selected_mfs/mfs.tsv", sep="\t", index=False)
    mf_cpg_pairs.to_csv(output_folder / "evaluation_input/mf_cpg_pairs.tsv", sep="\t", index=False)
    pd.Series(sorted(background_gene_ids)).to_csv(output_folder / "prepared/regulator_names.tsv", index=False, header=False)
    pd.Series(sorted(background_cpg_ids)).to_csv(output_folder / "prepared/target_names.tsv", index=False, header=False)
    print(f"{len(selected_mfs)} MFs, {len(mf_cpg_pairs)} MF-CpG pairs")


if __name__ == "__main__":
    main()
