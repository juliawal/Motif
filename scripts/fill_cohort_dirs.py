"""Split PanCan matrices by project and move them into cohort directories."""

from argparse import ArgumentParser
import contextlib
import shutil
from pathlib import Path
import yaml
from data_processing.sample_preparation import read_sample_ids, find_shared_samples, load_usable_metadata, write_shared_samples


def split_data_to_cohorts(input_path, data_type, cohort_to_samples, all_samples_in_file, output_dir):
    if data_type == "methylation":
        id_column = "cg_id"
    else:
        id_column = "gene_id"
    sample_id_to_column_number = {sample: index + 1 for index, sample in enumerate(all_samples_in_file)}

    with contextlib.ExitStack() as stack: # >> close all cohort files in the end, also if an error occurs
        # Open one file per cohort before reading because every input row is written to all of files at the same time.
        cohort_outputs = []
        for cohort, samples in cohort_to_samples.items():
            out_path = output_dir / f"{data_type}_{cohort}.tsv"
            output_file = stack.enter_context(open(out_path, "w", newline=""))
            output_file.write("\t".join([id_column] + samples) + "\n")
            cohort_outputs.append((output_file, [sample_id_to_column_number[sample] for sample in samples]))

        # Read the data matrix one line at a time.
        input_file = stack.enter_context(open(input_path, "rt"))
        next(input_file) # skip the header
        for line in input_file:
            row = line.rstrip("\r\n").split("\t")
            for output_file, columns in cohort_outputs:
                output_file.write(row[0] + "\t" + "\t".join([row[i] for i in columns]) + "\n") # row[0] is feature_id (gene or CpG)


def main():
    parser = ArgumentParser()
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--expression", required=True)
    parser.add_argument("--methylation", required=True)
    parser.add_argument("--raw-counts", required=True)
    parser.add_argument("--config", type=Path, default=Path("config.yaml"))
    parser.add_argument("--split-dir", type=Path, default=Path("splitted_by_project"))
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    split_config = config["split"]

    input_files_by_data_type = {"expression": args.expression, "methylation": args.methylation, "raw_counts": args.raw_counts}
    args.split_dir.mkdir(parents=True, exist_ok=True)

    metadata = load_usable_metadata(args.metadata)
    samples_by_data_source = {"metadata": list(metadata.index)}
    for data_type, path in input_files_by_data_type.items():
        samples_by_data_source[data_type] = read_sample_ids(path)
    samples_present_in_all_data_sources = find_shared_samples(samples_by_data_source, order_from="expression")
    cohort_to_samples_dict = {}
    for sample in samples_present_in_all_data_sources:
        cohort = metadata.loc[sample, "project_id"]
        if cohort not in cohort_to_samples_dict:
            cohort_to_samples_dict[cohort] = []
        cohort_to_samples_dict[cohort].append(sample)

    for data_type, path in input_files_by_data_type.items():
        split_data_to_cohorts(path, data_type, cohort_to_samples_dict, samples_by_data_source[data_type], args.split_dir)

    for cohort, samples in cohort_to_samples_dict.items():
        write_shared_samples(args.split_dir / f"shared_samples_{cohort}.tsv", samples)

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
