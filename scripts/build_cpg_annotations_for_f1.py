"""Build CpG island, gene and nearest TSS annotations for all CpGs (for binding-site evaluation)."""

from argparse import ArgumentParser
import gzip
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd


def read_gene_regions_and_transcript_starts(gencode_file):
    gene_regions = defaultdict(list)
    transcript_starts = defaultdict(list)
    with gzip.open(gencode_file, "rt") as gencode:
        for line in gencode:
            if line.startswith("#"):
                continue
            gencode_fields = line.split("\t")
            if gencode_fields[2] not in ("gene", "transcript"):
                continue
            chromosome = gencode_fields[0]
            region_start = int(gencode_fields[3]) - 1 # GENCODE uses 1 based coordinates with start and end included; >> convert to 0 based coordinates, include start, exclude end (<< Zhou CpG annotation)
            region_end = int(gencode_fields[4])
            strand = gencode_fields[6]
            if gencode_fields[2] == "gene":
                gene_attributes = gencode_fields[8].split(";")
                gene_name = gene_attributes[3].split('"')[1]
                gene_regions[chromosome].append((region_start, region_end, gene_name))
            if gencode_fields[2] == "transcript":
                tss = region_start if strand == "+" else region_end - 1
                strand_direction = 1 if strand == "+" else -1
                transcript_starts[chromosome].append((tss, strand_direction))
    return gene_regions, transcript_starts


def calculate_distances_to_nearest_tss(cpg_positions, tss_positions, tss_strands):
    distances = []
    for cpg_position in cpg_positions:
        distances_to_tss = abs(cpg_position - tss_positions)
        shortest_distance = distances_to_tss.min()
        nearest_tss_positions = tss_positions[distances_to_tss == shortest_distance]
        nearest_tss_position = nearest_tss_positions[0]
        distance = cpg_position - nearest_tss_position
        strands_at_nearest_tss = tss_strands[tss_positions == nearest_tss_position]
        nearest_tss_strand = strands_at_nearest_tss[0]
        if nearest_tss_strand == -1: # if the TSS is on the minus strand, reverse the sign of the distance 
            distance = -distance
        distances.append(distance)
    return np.array(distances)


def main():
    parser = ArgumentParser()
    parser.add_argument("--cpg-annotations-file", type=Path, required=True, help="Zhou v36 table with CpGs and CGIposition (cpg_island_context).")
    parser.add_argument("--gencode-file", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cpg_annotations = pd.read_csv(args.cpg_annotations_file, sep="\t", usecols=["probeID", "CpG_chrm", "CpG_beg", "CGIposition"]).rename(columns={"probeID": "cpg_id"})
    gene_regions, transcript_starts = read_gene_regions_and_transcript_starts(args.gencode_file) 
    cpg_annotations["genes_that_contain_cpg"] = pd.NA
    cpg_tss_distances = pd.Series(pd.NA, index=cpg_annotations.index)
    chromosome_tables = []
    for chromosome, cpgs_on_chromosome in cpg_annotations.groupby("CpG_chrm"):
        cpg_indices = cpgs_on_chromosome.index 
        tss_positions = []
        tss_strands = []
        for tss_position, tss_strand in sorted(transcript_starts[chromosome]):
            tss_positions.append(tss_position)
            tss_strands.append(tss_strand)
        tss_positions = np.array(tss_positions)
        tss_strands = np.array(tss_strands)
        cpg_positions = cpgs_on_chromosome["CpG_beg"]
        distances_to_nearest_tss = calculate_distances_to_nearest_tss(cpg_positions, tss_positions, tss_strands)
        cpg_tss_distances[cpg_indices] = distances_to_nearest_tss
        gene_starts, gene_ends, gene_names = [np.array(values) for values in zip(*gene_regions[chromosome])]
        genes_with_cpg_inside = []
        for cpg in cpg_positions:
            cpg_is_inside_gene = (gene_starts <= cpg) & (cpg < gene_ends)
            genes_with_this_cpg_inside = gene_names[cpg_is_inside_gene]
            unique_gene_names = set(genes_with_this_cpg_inside)
            gene_names_text = ";".join(unique_gene_names)
            genes_with_cpg_inside.append(gene_names_text)
        cpgs_on_chromosome["genes_that_contain_cpg"] = genes_with_cpg_inside
        chromosome_tables.append(cpgs_on_chromosome)
    cpgs_without_chromosome = cpg_annotations[cpg_annotations["CpG_chrm"].isna()]
    chromosome_tables.append(cpgs_without_chromosome)
    cpg_annotations = pd.concat(chromosome_tables).sort_index()
    cpg_annotations["distance_to_nearest_tss"] = cpg_tss_distances
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cpg_annotations[["cpg_id", "CGIposition", "genes_that_contain_cpg", "distance_to_nearest_tss"]].to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
