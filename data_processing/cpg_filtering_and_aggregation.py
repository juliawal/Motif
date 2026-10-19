"""Map CpGs to gene regions, then filter or aggregate them."""

import gzip
import re
from collections import defaultdict
import pandas as pd
from .expression_data_preparation import strip_ensembl_version


def build_cpg_gene_mapping(cpg_methylation, should_filter_cpgs, cpg_filter_genes_file, cpg_probe_coordinates_file, gene_regions_gtf_file, keep_promoters, keep_bodies, promoter_length, promoter_overlap):
    """Match CpGs to the selected promoter and gene-body regions."""
    gene_ids = None  # without a filter list: use all annotated genes
    if should_filter_cpgs:
        with open(cpg_filter_genes_file) as file:
            gene_ids = [line.strip() for line in file if line.strip()]
        gene_ids = strip_ensembl_version(gene_ids).unique().tolist()

    gene_region_coordinates = get_gene_coordinates(gene_ids, gene_regions_gtf_file, keep_promoters=keep_promoters, keep_bodies=keep_bodies, promoter_length=promoter_length, promoter_overlap=promoter_overlap)
    cpg_coordinates = get_cpg_coordinates(cpg_methylation, cpg_probe_coordinates_file)
    cpgs_by_chromosome = {chromosome: chromosome_cpgs for chromosome, chromosome_cpgs in cpg_coordinates.groupby("chrom")}

    cpg_gene_region_matches = []
    for gene_region in gene_region_coordinates.itertuples(index=False):
        chromosome_cpgs = cpgs_by_chromosome.get(gene_region.chrom)
        if chromosome_cpgs is None:
            continue
        in_gene_region = ((chromosome_cpgs["chromStart"] >= gene_region.chromStart) & (chromosome_cpgs["chromEnd"] <= gene_region.chromEnd))
        cpg_ids_in_gene_region = chromosome_cpgs.loc[in_gene_region, "name"]
        for cpg_id in cpg_ids_in_gene_region: 
            cpg_gene_region_matches.append({"aggregated_cpg": f"cpgs_agg_{gene_region.ensembl_gene_id}", "ensembl_gene_id": gene_region.ensembl_gene_id, "cpg_id": cpg_id, "region_type": gene_region.region_type})
    if not cpg_gene_region_matches:
        raise RuntimeError("No CpG gene-region overlaps were found.")

    cpg_gene_mapping = pd.DataFrame(cpg_gene_region_matches).drop_duplicates(ignore_index=True) # count each CpG once, also if it overlaps a promoter and a gene body
    cpg_gene_mapping["n_cpgs"] = cpg_gene_mapping.groupby("aggregated_cpg")["cpg_id"].transform("nunique")
    return cpg_gene_mapping


def get_cpg_coordinates(cpg_methylation, cpg_probe_coordinates_file):
    cpg_probe_coordinates = pd.read_csv(cpg_probe_coordinates_file, sep="\t", usecols=["#id", "chrom", "chromStart", "chromEnd"], low_memory=False)
    cpg_probe_coordinates = cpg_probe_coordinates.rename(columns={"#id": "cpg_id"})
    input_cpgs = pd.DataFrame({"cpg_id": list(cpg_methylation.columns), "_input_order": range(cpg_methylation.shape[1])})
    cpg_coordinates = input_cpgs.merge(cpg_probe_coordinates[["cpg_id", "chrom", "chromStart", "chromEnd"]], how="inner", on="cpg_id", validate="one_to_one").sort_values("_input_order")
    cpg_coordinates["chromStart"] = pd.to_numeric(cpg_coordinates["chromStart"], errors="raise").astype(int)
    cpg_coordinates["chromEnd"] = pd.to_numeric(cpg_coordinates["chromEnd"], errors="raise").astype(int)
    return cpg_coordinates[["chrom", "chromStart", "chromEnd", "cpg_id"]].rename(columns={"cpg_id": "name"})


def get_gene_coordinates(gene_ids, gene_regions_gtf_file, keep_promoters, keep_bodies, promoter_length, promoter_overlap):
    """Return promoter and gene-body intervals; gene_ids=None includes all genes."""
    gene_ids_to_keep = set(gene_ids) if gene_ids is not None else None
    gene_body_coordinates_by_gene = {}
    transcription_starts_by_gene = defaultdict(set)
    with gzip.open(gene_regions_gtf_file, "rt") as handle:
        for line in handle:
            if line.startswith("#"): # skip commentar fields
                continue
            fields = line.rstrip("\n").split("\t")
            chromosome, _, feature_type, start, end, _, strand, _, attributes = fields
            if feature_type not in ("gene", "transcript"):
                continue
            gene_id_search = re.search(r'(?:^|;\s*)gene_id "([^"]+)"', attributes)
            gene_id = strip_ensembl_version([gene_id_search.group(1)])[0]
            start = int(start) - 1 # convert interval to 0 based coordinates
            end = int(end)
            if feature_type == "gene":
                gene_body_coordinates_by_gene[gene_id] = (chromosome, start, end, strand)
            elif feature_type == "transcript":
                transcription_start = start if strand == "+" else end - 1
                transcription_starts_by_gene[gene_id].add(transcription_start)

    if gene_ids is None:
        gene_ids = list(gene_body_coordinates_by_gene)
    gene_regions = []
    for gene_id in gene_ids:
        chromosome, body_start, body_end, strand = gene_body_coordinates_by_gene[gene_id]
        if keep_promoters:
            for transcription_start in sorted(transcription_starts_by_gene[gene_id]):
                if strand == "+":
                    promoter_start = max(0, transcription_start - promoter_length)
                    promoter_end = transcription_start + promoter_overlap
                else:
                    promoter_start = max(0, transcription_start - promoter_overlap + 1)
                    promoter_end = transcription_start + promoter_length + 1
                gene_regions.append({"ensembl_gene_id": gene_id, "chrom": chromosome, "chromStart": promoter_start, "chromEnd": promoter_end, "region_type": "promoter"})
        if keep_bodies:
            gene_regions.append({"ensembl_gene_id": gene_id, "chrom": chromosome, "chromStart": body_start, "chromEnd": body_end, "region_type": "body"})
    return pd.DataFrame(gene_regions)



def aggregate_cpgs_per_gene(cpg_methylation, cpg_gene_mapping):
    aggregated_cpgs = {}
    for aggregated_cpg, mapping_rows_for_aggregated_cpg in cpg_gene_mapping.groupby("aggregated_cpg", sort=False):
        cpg_ids_in_gene_region = list(dict.fromkeys(mapping_rows_for_aggregated_cpg["cpg_id"]))
        aggregated_cpgs[aggregated_cpg] = cpg_methylation.loc[:, cpg_ids_in_gene_region].mean(axis=1, skipna=False)
    return pd.DataFrame(aggregated_cpgs, index=cpg_methylation.index)


def filter_cpgs(cpg_methylation, cpg_gene_mapping):
    mapped_cpg_ids = set(cpg_gene_mapping["cpg_id"])
    columns = [column for column in cpg_methylation.columns if column in mapped_cpg_ids]
    filtered_cpgs = cpg_methylation.loc[:, columns].copy()
    return filtered_cpgs
