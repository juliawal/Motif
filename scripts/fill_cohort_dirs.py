"""Move (and split) PanCan files into cohort directories."""

from argparse import ArgumentParser
import shutil
from pathlib import Path
import yaml


def main():
    parser = ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("config.yaml"))
    parser.add_argument("--split-dir", type=Path, default=Path("splitted_by_project"))
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    split_config = config["split"]

    for cohort in sorted(config["cohorts"]):
        is_split_cohort = cohort in split_config["cohorts"]
        project = cohort
        if is_split_cohort:
            project = split_config["project"]
        cohort_dir = args.data_dir / cohort
        cohort_dir.mkdir(parents=True, exist_ok=True)
        for matrix_name in ("expression", "methylation", "raw_counts", "shared_samples"):
            source_file = args.split_dir / f"{matrix_name}_TCGA-{project}.tsv"
            output_file = cohort_dir / f"{matrix_name}.tsv"
            if is_split_cohort:
                shutil.copy2(source_file, output_file) # copy shared files so that the other SKCM cohort can read them too; both SKCM cohorts use the files from the same project
            else:
                shutil.move(str(source_file), str(output_file))


if __name__ == "__main__":
    main()
