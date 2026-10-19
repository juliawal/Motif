"""Download ReMap peaks for MFs without a JASPAR PWM and record their availability."""

from argparse import ArgumentParser
import tempfile
from pathlib import Path
import pandas as pd
import requests
import yaml
from evaluation.binding_site.select_binding_motifs import find_motifs_for_mfs_in_jaspar_database
from data_processing.expression_data_preparation import load_symbol_aliases 


def download_peak_file_for_given_mf(mf_symbol, chipseq_dir):
    """Return the peak file, or None when ReMap has no entry for this MF."""
    peak_file = chipseq_dir / f"{mf_symbol}.bed.gz"
    if peak_file.exists(): # reuse peak files that were already downloaded for this gene symbol
        return peak_file

    peak_file_url = f"https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/TF/{mf_symbol}/remap2022_{mf_symbol}_nr_macs2_hg38_v1_0.bed.gz"
    response = requests.get(peak_file_url)
    if response.status_code == 404: # no ReMap entry for this MF
        return None
    response.raise_for_status()
    if not response.content:
        raise ValueError(f"{mf_symbol}: ReMap returned an empty peak file")

    chipseq_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir = chipseq_dir) as temporary_dir:
        temporary_path = Path(temporary_dir) / peak_file.name
        with open(temporary_path, "wb") as output_file:
            output_file.write(response.content)
        temporary_path.replace(peak_file)
    return peak_file
 


def build_chipseq_overview(mf_symbols, chipseq_dir, gene_symbol_aliases):
    overview_columns = ["MF", "track", "peak_is_available", "remap_symbol", "peak_file", "source_url", "release", "genome_version", "symbols_tried"]
    overview_rows = []
    for mf_symbol in sorted(set(mf_symbols)):
        candidate_symbols = list(dict.fromkeys([mf_symbol, *gene_symbol_aliases.get(mf_symbol, [])])) # try the mf symbol first, then the mf symbols aliases
        symbols_tried = []
        overview_row = {"MF": mf_symbol, "track": "ChIP_seq", "peak_is_available": False, "remap_symbol": pd.NA, "peak_file": pd.NA, "source_url": pd.NA, "release": "ReMap2022", "genome_version": "hg38",}
        for mf_symbol in candidate_symbols:
            symbols_tried.append(mf_symbol)
            peak_file = download_peak_file_for_given_mf(mf_symbol, chipseq_dir)
            if peak_file is not None:
                overview_row["peak_is_available"] = True
                overview_row["remap_symbol"] = mf_symbol
                overview_row["peak_file"] = peak_file.name
                overview_row["source_url"] = f"https://remap.univ-amu.fr/storage/remap2022/hg38/MACS2/TF/{mf_symbol}/remap2022_{mf_symbol}_nr_macs2_hg38_v1_0.bed.gz"
                break
        overview_row["symbols_tried"] = ";".join(symbols_tried)
        overview_rows.append(overview_row)
    return pd.DataFrame(overview_rows, columns = overview_columns)


def main():
    parser = ArgumentParser()
    parser.add_argument("--from-pipeline-results", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("config.yaml"))
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())

    # Select all MFs.
    mf_symbols = set()
    for selection_rank in range(1, config["n_selected_stable_settings"] + 1):
        for cohort in config["cohorts"]:
            mf_cpg_pairs_file = args.from_pipeline_results / "selected_settings" / f"rank_{selection_rank}" / cohort / "evaluation_input/mf_cpg_pairs.tsv"
            cohort_mf_symbols = pd.read_csv(mf_cpg_pairs_file, sep="\t", usecols=["MF"])["MF"]
            mf_symbols.update(cohort_mf_symbols.dropna())
    gene_symbol_aliases = load_symbol_aliases(config["gene_annotation_file"])
    _, mfs_without_pwm = find_motifs_for_mfs_in_jaspar_database(sorted(mf_symbols), config["jaspar_meme"], jaspar_sqlite_file = config["jaspar_sqlite"], gene_symbol_aliases = gene_symbol_aliases)

    # Download available peaks and note MFs with missing peaks.
    chipseq_overview = build_chipseq_overview(mfs_without_pwm, Path(config["chipseq_dir"]), gene_symbol_aliases)
    Path(config["chipseq_overview"]).parent.mkdir(parents=True, exist_ok=True)
    chipseq_overview.to_csv(config["chipseq_overview"], sep="\t", index=False)


if __name__ == "__main__":
    main()
