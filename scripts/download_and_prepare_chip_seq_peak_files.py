from argparse import ArgumentParser
import gzip
import json
from pathlib import Path
import shutil
from urllib.parse import quote
from urllib.request import urlopen
import pandas as pd
import yaml
from data_processing.expression_data_preparation import load_symbol_aliases
from evaluation.binding_site.select_binding_motifs import find_motifs_for_mfs_in_jaspar_database

remap_file = "remap2022_nr_macs2_hg38_v1_0.bed.gz"
remap_url = "https://zenodo.org/records/10527088/files/" + remap_file + "?download=1"


def download_remap_file(output_file):
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with urlopen(remap_url, timeout=120) as download, output_file.open("wb") as remap_file:
        shutil.copyfileobj(download, remap_file)


def select_mfs(config, mf_cpg_pair_files):
    mf_symbols = set()
    for mf_cpg_pair_file in mf_cpg_pair_files:
        mf_cpg_pairs = pd.read_csv(mf_cpg_pair_file, sep="\t", usecols=["MF"])
        for mf_symbol in mf_cpg_pairs["MF"]:
            mf_symbols.add(mf_symbol)
    gene_symbol_aliases = load_symbol_aliases(config["gene_annotation_file"])
    _, mfs_without_pwm = find_motifs_for_mfs_in_jaspar_database(mf_symbols, config["jaspar_meme"], jaspar_sqlite_file=config["jaspar_sqlite"], gene_symbol_aliases=gene_symbol_aliases)
    return mf_symbols, mfs_without_pwm, gene_symbol_aliases


def save_peaks_per_mf(remap_file, required_symbols, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    peak_files_by_symbol = {}
    with gzip.open(remap_file, "rt") as remap:
        for peak_line in remap:
            gene_symbol = peak_line.split("\t")[3].split(":")[0]
            if gene_symbol in required_symbols:
                peak_filename = quote(gene_symbol, safe="") + ".bed" # encode / (and other characters) for the filename
                peak_file = output_dir / peak_filename
                file_mode = "a" # append peak to the BED file
                if gene_symbol not in peak_files_by_symbol:
                    file_mode = "w" # write first peak the BED file
                with peak_file.open(file_mode) as local_file:
                    local_file.write(peak_line)
                peak_files_by_symbol[gene_symbol] = peak_filename
    return peak_files_by_symbol


def prepare_chipseq(config, mf_cpg_pair_files, output_dir):

    # Collect MFs with alternative gene names.
    mf_symbols, mfs_without_pwm, gene_symbol_aliases = select_mfs(config, mf_cpg_pair_files)
    symbols_by_mf = {}
    required_symbols = set()
    for mf_symbol in mfs_without_pwm:
        symbols_by_mf[mf_symbol] = [mf_symbol]
        if mf_symbol in gene_symbol_aliases:
            for alias in gene_symbol_aliases[mf_symbol]:
                if alias not in symbols_by_mf[mf_symbol]:
                    symbols_by_mf[mf_symbol].append(alias)
        for symbol in symbols_by_mf[mf_symbol]:
            required_symbols.add(symbol)

    peak_files_by_symbol = save_peaks_per_mf(config["remap_file"], required_symbols, Path(output_dir) / "peaks")

    # Save one ChIP peak file per MF if available.
    chipseq_peak_files_by_mf = {}
    for mf_symbol in mf_symbols:
        chipseq_peak_files_by_mf[mf_symbol] = None
    for mf_symbol in mfs_without_pwm:
        for symbol in symbols_by_mf[mf_symbol]:
            if symbol in peak_files_by_symbol:
                chipseq_peak_files_by_mf[mf_symbol] = "peaks/" + peak_files_by_symbol[symbol]
                break

    with (Path(output_dir) / "peak_files.json").open("w") as local_file:
        json.dump(chipseq_peak_files_by_mf, local_file)
    return chipseq_peak_files_by_mf


def load_chipseq_peak_files(peak_mapping_file, mf_symbols):
    peak_mapping_file = Path(peak_mapping_file)
    with peak_mapping_file.open() as peak_mapping_json:
        peak_files_by_mf = json.load(peak_mapping_json)
    chipseq_peak_files_by_mf = {}
    for mf_symbol in mf_symbols:
        peak_mapping_json = peak_files_by_mf[mf_symbol]
        if peak_mapping_json is None:
            chipseq_peak_files_by_mf[mf_symbol] = None
        else:
            chipseq_peak_files_by_mf[mf_symbol] = peak_mapping_file.parent / peak_mapping_json
    return chipseq_peak_files_by_mf


def main():
    parser = ArgumentParser()
    parser.add_argument("command", choices=["download", "prepare"])
    parser.add_argument("--config", type=Path, default=Path("config.yaml"))
    parser.add_argument("--pairs", type=Path, nargs="+")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    with args.config.open() as config_file:
        config = yaml.safe_load(config_file)
    if args.command == "download":
        download_remap_file(config["remap_file"])
    if args.command == "prepare":
        prepare_chipseq(config, args.pairs, args.output_dir)


if __name__ == "__main__":
    main()
