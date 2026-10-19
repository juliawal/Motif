import pandas as pd


def read_cpg_methylation_blocks(cpg_methylation_file, cpg_ids):
    """Read the methylation matrix in blocks to save memory."""
    requested_cpg_ids = set(cpg_ids)
    selected_cpg_methylation_blocks = []
    cpg_methylation_blocks = pd.read_csv(cpg_methylation_file, sep="\t", index_col=0, chunksize=20_000)
    for cpg_methylation_block in cpg_methylation_blocks:
        cpgs_to_keep = cpg_methylation_block.index.isin(requested_cpg_ids)
        selected_cpg_methylation = cpg_methylation_block.loc[cpgs_to_keep]
        if not selected_cpg_methylation.empty:
            selected_cpg_methylation_blocks.append(selected_cpg_methylation)
    return selected_cpg_methylation_blocks