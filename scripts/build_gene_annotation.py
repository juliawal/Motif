"""Build gene annotations, transcript TSS annotations and coding-gene IDs."""

from argparse import ArgumentParser
import gzip
from pathlib import Path
import pandas as pd
from data_processing.expression_data_preparation import strip_ensembl_version


def read_gencode_genes(gencode_file):
    gene_rows = []
    with gzip.open(gencode_file, "rt") as gencode:
        for gencode_line in gencode:
            if gencode_line.startswith("#"):
                continue
            gencode_fields = gencode_line.rstrip("\n").split("\t")
            if gencode_fields[2] != "gene":
                continue
            gene_attributes = gencode_fields[8].split(";")
            gene_id = gene_attributes[0].split('"')[1]
            gencode_gene_name = gene_attributes[3].split('"')[1]
            chromosome = gencode_fields[0] 
            gene_type = gene_attributes[1].split('"')[1]
            gene_strand = gencode_fields[6]
            gene_start, gene_end = int(gencode_fields[3]), int(gencode_fields[4])
            gene_tss = gene_start if gene_strand == "+" else gene_end
            gene_rows.append({"ensembl_id": gene_id, "gencode_gene_name": gencode_gene_name, "chromosome": chromosome, "gene_tss": gene_tss, "gene_type": gene_type})
    gencode_genes = pd.DataFrame(gene_rows)
    gencode_genes["ensembl_id"] = strip_ensembl_version(gencode_genes["ensembl_id"])
    return gencode_genes


def match_gencode_genes_with_hgnc_gene_information(gencode_genes, hgnc_file, output_path):
    """Add gene symbols, aliases, previous symbols and Entrez IDs from HGNC to the GENCODE genes.
    If the Ensembl ID has no linked Entrez ID, look up the gene symbol in HGNC, where it may be linked to a different Ensembl ID.
    Use the Entrez ID from that HGNC entry, older gene symbols and aliases for the original GENCODE Ensembl ID.""" 

    hgnc_genes = pd.read_csv(hgnc_file, sep="\t", dtype=str, usecols=["symbol", "alias_symbol", "prev_symbol", "entrez_id", "ensembl_gene_id"])
    hgnc_genes = hgnc_genes.rename(columns={"ensembl_gene_id": "ensembl_id", "symbol": "selected_gene_symbol", "alias_symbol": "gene_symbol_aliases", "prev_symbol": "previous_gene_symbols"})
    gencode_gene_names_by_ensembl_id = dict(zip(gencode_genes["ensembl_id"], gencode_genes["gencode_gene_name"]))
    # For each Ensembl ID, use the matching GENCODE name, or the first HGNC entry if none matches.
    selected_hgnc_rows = []
    for gene_id, hgnc_entries in hgnc_genes.groupby("ensembl_id"):
        gencode_gene_name = gencode_gene_names_by_ensembl_id.get(gene_id)
        hgnc_entries_with_symbol_matching_gencode_name = hgnc_entries[hgnc_entries["selected_gene_symbol"] == gencode_gene_name]
        if not hgnc_entries_with_symbol_matching_gencode_name.empty:
            selected_hgnc_rows.append(hgnc_entries_with_symbol_matching_gencode_name.iloc[0])
        else:
            selected_hgnc_rows.append(hgnc_entries.iloc[0])
    selected_hgnc_entries = pd.DataFrame(selected_hgnc_rows)

    hgnc_metadata_columns = ["ensembl_id", "selected_gene_symbol", "entrez_id", "previous_gene_symbols", "gene_symbol_aliases"]
    hgnc_gene_information = selected_hgnc_entries[hgnc_metadata_columns]
    genes_with_hgnc_information = gencode_genes.merge(hgnc_gene_information, how="left", on="ensembl_id") # add HGNC information by Ensembl ID; left >> keep all GENCODE genes
    hgnc_symbols = genes_with_hgnc_information["selected_gene_symbol"]
    gencode_gene_names = genes_with_hgnc_information["gencode_gene_name"]
    genes_with_hgnc_information["selected_gene_symbol"] = hgnc_symbols.fillna(gencode_gene_names) # use the GENCODE name when the HGNC symbol is missing

    # If the Entrez ID is missing, look up the gene symbol in HGNC and see if the gene symbol has an Entrez ID in HGNC.
    hgnc_entries_by_symbol = hgnc_genes.set_index("selected_gene_symbol")
    genes_with_missing_entrez_id = genes_with_hgnc_information["entrez_id"].isna() & genes_with_hgnc_information["selected_gene_symbol"].isin(hgnc_entries_by_symbol.index)
    hgnc_entries_matched_by_symbol = hgnc_entries_by_symbol.loc[genes_with_hgnc_information.loc[genes_with_missing_entrez_id, "selected_gene_symbol"]]
    for metadata_column in ["entrez_id", "previous_gene_symbols", "gene_symbol_aliases"]:
        genes_with_hgnc_information.loc[genes_with_missing_entrez_id, metadata_column] = hgnc_entries_matched_by_symbol[metadata_column].to_numpy()

    column_names = ["ensembl_id", "selected_gene_symbol", "entrez_id", "chromosome", "gene_tss", "gene_type", "previous_gene_symbols", "gene_symbol_aliases"]
    gene_annotation = genes_with_hgnc_information[column_names]
    gene_annotation.to_csv(output_path, sep="\t", index=False)


def build_gene_tss_annotations(gencode_file, output_path): 
    transcript_rows = []
    with gzip.open(gencode_file, "rt") as gencode:
        for gencode_line in gencode:
            if gencode_line.startswith("#"):
                continue
            gencode_fields = gencode_line.rstrip("\n").split("\t")
            if gencode_fields[2] != "transcript":
                continue
            gene_id = gencode_fields[8].split('"')[1]
            chromosome = gencode_fields[0]
            gene_strand = gencode_fields[6] 
            transcript_start, transcript_end = int(gencode_fields[3]), int(gencode_fields[4])
            if gene_strand == "+":
                gene_tss = transcript_start
            if gene_strand == "-":
                gene_tss = transcript_end
            transcript_rows.append({"gene_id": gene_id, "chromosome": chromosome, "gene_tss": gene_tss, "gene_strand": gene_strand})
    gencode_transcripts = pd.DataFrame(transcript_rows)
    gencode_transcripts["gene_id"] = strip_ensembl_version(gencode_transcripts["gene_id"])
    gene_tss_annotations = gencode_transcripts.drop_duplicates() # we store only gene ID, chromosome, tss and strand; multiple transcripts of the same gene can start at the same tss; so we could have identical rows; >> .drop_duplicates()
    gene_tss_annotations.to_csv(output_path, sep="\t", index=False)


def build_coding_gene_list(gencode_genes, output_path):
    protein_coding_genes = gencode_genes.loc[gencode_genes["gene_type"] == "protein_coding", "ensembl_id"]
    protein_coding_genes.to_csv(output_path, index=False, header=False)


def main():
    parser = ArgumentParser()
    parser.add_argument("--gencode-file", type=Path, required=True, help="Gencode v22 GTF in gzip format.")
    parser.add_argument("--hgnc-file", type=Path, required=True)
    parser.add_argument("--gene-annotation-output", type=Path, required=True)
    parser.add_argument("--tss", type=Path, required=True)
    parser.add_argument("--coding-genes", type=Path, required=True)
    args = parser.parse_args()
 
    gencode_genes = read_gencode_genes(args.gencode_file)
    match_gencode_genes_with_hgnc_gene_information(gencode_genes, args.hgnc_file, args.gene_annotation_output)
    build_gene_tss_annotations(args.gencode_file, args.tss)
    build_coding_gene_list(gencode_genes, args.coding_genes)


if __name__ == "__main__":
    main()
