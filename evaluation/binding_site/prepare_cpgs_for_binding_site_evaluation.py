import os
import tempfile
from pathlib import Path
import pandas as pd
from pyfaidx import Fasta


def load_cpg_coordinates(cpg_annotations_file):
    cpg_annotations = pd.read_csv(cpg_annotations_file, sep="\t").rename(columns={"probeID": "cpg_id"})
    cpg_annotations = cpg_annotations.dropna(subset=["CpG_chrm", "CpG_beg", "CpG_end"]) # remove CpGs with missing CpG_chrm, CpG_beg or CpG_end values
    cpg_annotations["CpG_beg"] = cpg_annotations["CpG_beg"].astype(int)
    cpg_annotations["CpG_end"] = cpg_annotations["CpG_end"].astype(int) 
    return cpg_annotations.set_index("cpg_id")[["CpG_chrm", "CpG_beg", "CpG_end"]]


def get_cpg_annotations(cpg_ids, cpg_gene_annotations, promoter_length=1500, promoter_overlap=200):
    cpg_ids = list(cpg_ids)
    cpg_gene_annotations = cpg_gene_annotations.reindex(cpg_ids)
    cpg_annotations = pd.DataFrame(index=cpg_ids) 
    cpg_annotations["island_context"] = cpg_gene_annotations["CGIposition"].fillna("OpenSea").replace(["N_Shore", "S_Shore", "N_Shelf", "S_Shelf"], "Shore_Shelf") 
    nearest_tss_distances = cpg_gene_annotations["distance_to_nearest_tss"]
    cpg_is_in_promoter = (nearest_tss_distances >= -promoter_length) & (nearest_tss_distances < promoter_overlap)
    cpg_is_in_gene = cpg_gene_annotations["genes_that_contain_cpg"].notna()
    cpg_annotations["gene_region"] = "Intergenic"
    cpg_annotations.loc[cpg_is_in_gene, "gene_region"] = "GeneBody"
    cpg_annotations.loc[cpg_is_in_promoter, "gene_region"] = "Promoter"
    return cpg_annotations


def extract_cpg_sequences(cpg_coordinates, hg38_genome, fasta_index, window_around_cpg=250):
    hg38_genome = Path(hg38_genome)
    genome_sequences = Fasta(str(hg38_genome), indexname = str(fasta_index), as_raw = True, sequence_always_upper = True, rebuild = False, build_index = False)

    # Exctract DNA sequence around each CpG and store it with its window coordinates.
    cpg_sequences = {}
    for cpg_id, cpg_coordinate in cpg_coordinates.iterrows():
        chromosome = cpg_coordinate["CpG_chrm"]
        cpg_start = int(cpg_coordinate["CpG_beg"])
        window_start = max(0, cpg_start - window_around_cpg)
        window_end = cpg_start + window_around_cpg + 2
        cpg_sequence = genome_sequences[chromosome][window_start:window_end]
        if "N" in cpg_sequence and cpg_sequence.count("N") / len(cpg_sequence) > 0.1: # skip CpG sequences with more than 10 % unknown bases
            continue
        cpg_sequences[cpg_id] = {"chrom": chromosome, "start": window_start, "end": window_end, "sequence": cpg_sequence}
    return cpg_sequences


def write_cpg_sequences_as_fasta(cpg_sequences, output_path):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as fasta_file:
        for cpg_id, cpg_sequence_window in cpg_sequences.items():
            fasta_file.write(f">{cpg_id}\n{cpg_sequence_window['sequence']}\n")
    return output_path


def compute_gc_fraction_by_cpg(cpg_sequences):
    gc_fraction_by_cpg = {}
    for cpg_id, cpg_sequence_window in cpg_sequences.items():
        cpg_sequence = cpg_sequence_window["sequence"]
        n_valid_bases = len(cpg_sequence) - cpg_sequence.count("N") # N are unknown bases
        if n_valid_bases == 0:
            gc_fraction_by_cpg[cpg_id] = float("nan")
        else:
            gc_fraction_by_cpg[cpg_id] = (cpg_sequence.count("G") + cpg_sequence.count("C")) / n_valid_bases
    return pd.Series(gc_fraction_by_cpg, name="gc_fraction")


def compute_base_frequencies_for_fimo(cpg_sequences, output_path):
    """Calculate how often A, C, G and T are in the CpG windows and save their frequencies for FIMO."""
    base_counts = {base: 0 for base in "ACGT"}
    for cpg_sequence_window in cpg_sequences.values():
        cpg_sequence = cpg_sequence_window["sequence"]
        for base in base_counts:
            base_counts[base] += cpg_sequence.count(base)
    n_bases = sum(base_counts.values())
    base_frequencies = {base: base_counts[base] / n_bases for base in base_counts}
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    background_text = ""
    for base, base_frequency in base_frequencies.items():
        background_text += f"{base}\t{base_frequency}\n"
    with open(output_path, "w", newline="") as output_file:
        output_file.write(background_text)
    return base_frequencies
