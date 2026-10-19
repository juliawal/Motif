import csv
from pathlib import Path
import numpy as np
import pandas as pd


def read_sample_ids(matrix_path):
    path = Path(matrix_path)
    with open(path, "rt", newline="") as file:
        header = next(csv.reader(file, delimiter="\t"))
    return header[1:]


def read_shared_samples(path):
    path = Path(path)
    sample_table = pd.read_csv(path, sep="\t", dtype=str)
    return sample_table["sample_id"].tolist()


def write_shared_samples(path, samples):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.Series(list(samples), name="sample_id").to_csv(path, sep="\t", index=False)
    return path


def load_usable_metadata(metadata_path):
    metadata = pd.read_csv(metadata_path, sep="\t")
    metadata = metadata.rename(columns={"sampleID": "sample"})
    metadata = metadata.dropna(subset=["project_id"])
    metadata = metadata.dropna(subset=["Gender"])
    return metadata.drop_duplicates(subset="sample", keep="first").set_index("sample")


def find_shared_samples(samples_by_data_source, *, order_from):
    sample_sets = {name: set(values) for name, values in samples_by_data_source.items()}
    common_samples = set.intersection(*sample_sets.values())
    return [sample for sample in dict.fromkeys(samples_by_data_source[order_from]) if sample in common_samples]


def get_primary_tumour_type(sample):
    return str(sample).split("-")[3][:2]


def get_patient_id(sample):
    return "-".join(str(sample).split("-")[:3])


def select_tumour_samples(samples, primary_tumour_type, exclude_patients_of_type=None):
    """Select one sample of the given tumour type per patient; if a patient has more than one sample: keep the alphabetically first sample for each patient."""
    samples = list(samples)
    excluded_patients = set()
    if exclude_patients_of_type is not None:
        excluded_patients = {get_patient_id(sample) for sample in samples if get_primary_tumour_type(sample) == exclude_patients_of_type}
    matching_samples = [sample for sample in samples if get_primary_tumour_type(sample) == primary_tumour_type and get_patient_id(sample) not in excluded_patients]
    # Select the alphabetically first sample for each patient.
    chosen_patients = set()
    chosen_samples = set()
    for sample in sorted(matching_samples):
        patient = get_patient_id(sample)
        if patient not in chosen_patients:
            chosen_patients.add(patient)
            chosen_samples.add(sample)
    return pd.Index([sample for sample in samples if sample in chosen_samples])


def align_matrices(gene_expression, cpg_methylation, samples):
    """Align expression and methylation to shared-sample order."""
    return gene_expression.loc[samples], cpg_methylation.loc[samples]


def drop_incomplete_samples(gene_expression, cpg_methylation, sample_max_missing_share):
    """Remove samples if their missingness is more than the threshold allows."""
    missing_fractions = {}
    for data_type, data_matrix in (("gene expression", gene_expression), ("CpG methylation", cpg_methylation)):
        missing = ~np.isfinite(data_matrix.to_numpy())
        missing_fraction = missing.mean(axis=1)
        missing_fractions[data_type] = pd.Series(missing_fraction, index=data_matrix.index)
    missing_per_matrix = pd.DataFrame(missing_fractions)
    max_missing_fraction = missing_per_matrix.max(axis=1)
    within_threshold = max_missing_fraction <= sample_max_missing_share
    samples_to_keep = max_missing_fraction.index[within_threshold]
    return gene_expression.loc[samples_to_keep], cpg_methylation.loc[samples_to_keep]
