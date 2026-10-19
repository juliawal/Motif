from argparse import ArgumentParser
import urllib.request
from pathlib import Path
import yaml


def download_gene_library(library_name, gene_libraries_dir):
    """Save one enrichr library as a GMT file."""
    gene_library_url = f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={library_name}"
    with urllib.request.urlopen(gene_library_url, timeout=120) as response:
        gmt_content = response.read()
    gene_library_file = gene_libraries_dir / f"{library_name}.gmt"
    with open(gene_library_file, "wb") as output_file:
        output_file.write(gmt_content)
     

def main():
    parser = ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("resources/gene_sets"))
    parser.add_argument("--config", type=Path, default=Path("config.yaml"))
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text()) 
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for library_name in config["gene_set_enrichment_libraries"]:
        download_gene_library(library_name, args.output_dir)


if __name__ == "__main__":
    main()
