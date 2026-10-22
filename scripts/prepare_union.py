"""Create Union inputs from rank 1–3: python -m scripts.prepare_union."""

from pathlib import Path
import pandas as pd
import yaml


def main():
    config = yaml.safe_load(Path("config.yaml").read_text())
    for cohort in config["cohorts"]:
        selected_mf_tables = []
        mf_cpg_pair_tables = []
        for rank in (1, 2, 3):
            cohort_folder = Path(f"results/selected_settings/rank_{rank}") / cohort
            selected_mfs = pd.read_csv(cohort_folder / "selected_mfs/mfs.tsv", sep="\t", usecols=["MF_ensembl"])
            mf_cpg_pairs = pd.read_csv(cohort_folder / "evaluation_input/mf_cpg_pairs.tsv", sep="\t", usecols=["MF", "CpG", "MF_ensembl", "MF_chromosome"])
            selected_mf_tables.append(selected_mfs)
            mf_cpg_pair_tables.append(mf_cpg_pairs)

        union_mfs = pd.concat(selected_mf_tables).drop_duplicates("MF_ensembl")
        union_mf_cpg_pairs = pd.concat(mf_cpg_pair_tables)
        union_mf_cpg_pairs = union_mf_cpg_pairs.drop_duplicates(["MF_ensembl", "CpG"])

        output_file = Path("results/union") / cohort / "selected_mfs/mfs.tsv"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        union_mfs.to_csv(output_file, sep="\t", index=False)

        output_file = Path("results/union") / cohort / "evaluation_input/mf_cpg_pairs.tsv"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        union_mf_cpg_pairs.to_csv(output_file, sep="\t", index=False)
        print(f"{cohort}: {len(union_mfs)} MFs, {len(union_mf_cpg_pairs)} MF-CpG pairs")


if __name__ == "__main__":
    main()
