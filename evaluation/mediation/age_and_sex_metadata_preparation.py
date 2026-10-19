import numpy as np


def prepare_age_and_sex_metadata(prepared_samples, metadata):
    age_and_sex_metadata = metadata.loc[prepared_samples, ["age", "sex"]].copy()
    sample_has_complete_age_and_sex_metadata = (np.isfinite(age_and_sex_metadata["age"].to_numpy()) & age_and_sex_metadata["sex"].notna().to_numpy() & age_and_sex_metadata["sex"].ne("").to_numpy())
    complete_samples = prepared_samples[sample_has_complete_age_and_sex_metadata]
    age_and_sex_metadata = age_and_sex_metadata.loc[complete_samples]
    age_and_sex_metadata.index = complete_samples
    sex_metadata = age_and_sex_metadata["sex"]
    age_and_sex_metadata = age_and_sex_metadata.drop(columns="sex")
    if sex_metadata.nunique() > 1: # include sex only when both sexes are present 
        age_and_sex_metadata["sex_female"] = sex_metadata.eq("female") # true: female; false male
    age_and_sex_columns = age_and_sex_metadata.columns.tolist() 
    return complete_samples, age_and_sex_metadata, age_and_sex_columns
