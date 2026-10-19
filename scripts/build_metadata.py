from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
import yaml
from data_processing.sample_preparation import read_shared_samples, load_usable_metadata


def main():
    parser = ArgumentParser()
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("config.yaml"))
    parser.add_argument("--data-directory", type=Path, default=Path("data"))
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    metadata = load_usable_metadata(args.metadata)
    for cohort in config["cohorts"]:
        cohort_dir = args.data_directory / cohort
        sample_ids = read_shared_samples(cohort_dir / "shared_samples.tsv")
        cohort_metadata = metadata.loc[sample_ids]

        # Save sex metadata for sex filter.
        sample_metadata = pd.DataFrame({"sample_id": sample_ids})
        sample_metadata["sex"] = cohort_metadata["Gender"].str.lower().to_numpy()
        sample_metadata.to_csv(cohort_dir / "metadata.tsv", sep="\t", index=False)

        # Save age metadata for mediation. 
        sample_metadata["age"] = pd.to_numeric(cohort_metadata["Age at Diagnosis in Years"], errors="coerce").to_numpy()
        sample_metadata.to_csv(cohort_dir / "age_and_sex.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
