from argparse import ArgumentParser
from pathlib import Path
import tempfile
import pandas as pd
from data_processing.expression_data_preparation import load_symbol_aliases 
from evaluation.binding_site.select_background_cpgs import build_cpg_matching_categories, draw_random_matched_background_cpgs
from evaluation.binding_site.find_chipseq_overlaps import cpgs_with_peak_overlap, peak_file_for_mf
from evaluation.binding_site.prepare_cpgs_for_binding_site_evaluation import compute_gc_fraction_by_cpg, extract_cpg_sequences, get_cpg_annotations, load_cpg_coordinates, write_cpg_sequences_as_fasta, compute_base_frequencies_for_fimo
from evaluation.binding_site.compare_binding_hit_rates import create_binding_site_evaluation_results, calculate_hit_rate_sign_test, summarize_mf_hit_rates
from evaluation.binding_site.select_binding_motifs import find_motifs_for_mfs_in_jaspar_database, write_selected_motifs 
from evaluation.binding_site.find_motif_hits import map_fimo_hits_to_cpgs, run_fimo


def main():
    parser = ArgumentParser()
    parser.add_argument("--window-around-cpg", type=int, required=True)
    parser.add_argument("--bg-to-target-ratio", type=int, required=True)
    parser.add_argument("--n-gc-bins", type=int, required=True)
    parser.add_argument("--fimo-pvalue", type=float, required=True)
    parser.add_argument("--fimo-max-stored-scores", type=int, required=True)
    parser.add_argument("--promoter-length", type=int, required=True)
    parser.add_argument("--promoter-overlap", type=int, required=True) 
    parser.add_argument("--bg-sampling-seed", type=int, required=True)
    parser.add_argument("--jaspar-meme", type=Path, required=True)
    parser.add_argument("--jaspar-sqlite", type=Path, required=True)
    parser.add_argument("--chipseq-dir", type=Path, required=True)
    parser.add_argument("--mf-cpg-pairs", type=Path, required=True)
    parser.add_argument("--cpg-annotations-file", type=Path, required=True)  
    parser.add_argument("--cpg-names", type=Path, required=True)
    parser.add_argument("--cpg-gene-annotations-file", type=Path, required=True)
    parser.add_argument("--gene-annotation-file", type=Path, required=True)
    parser.add_argument("--hg38-genome", type=Path, required=True)
    parser.add_argument("--fasta-index", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--binding-sites-results-directory", type=Path, required=True)
    args = parser.parse_args()
    args.binding_sites_results_directory.mkdir(parents=True, exist_ok=True)

    # Prepare genes.
    mf_cpg_pairs = pd.read_csv(args.mf_cpg_pairs, sep="\t", usecols=["MF", "CpG"])
    mf_symbols = mf_cpg_pairs["MF"].unique().tolist()
    target_cpgs_by_mf = {mf_symbol: set(mf_cpg_pairs_for_mf["CpG"]) for mf_symbol, mf_cpg_pairs_for_mf in mf_cpg_pairs.groupby("MF")}
    gene_symbol_aliases = load_symbol_aliases(args.gene_annotation_file)
    pwm_ids_by_mf, mfs_without_pwm = find_motifs_for_mfs_in_jaspar_database(mf_symbols, args.jaspar_meme, args.jaspar_sqlite, gene_symbol_aliases = gene_symbol_aliases)

    # Prepare CpGs.
    cpg_coordinates = load_cpg_coordinates(args.cpg_annotations_file) # CpGs names have already been filtered in stage A
    cpg_names = pd.read_csv(args.cpg_names, header=None)[0].tolist()
    cpg_coordinates = cpg_coordinates.loc[cpg_coordinates.index.isin(cpg_names)]
    cpg_gene_annotations = pd.read_csv(args.cpg_gene_annotations_file, sep="\t").rename(columns={"probeID": "cpg_id", "probe_strand": "cpg_strand", "genesUniq": "genes_that_contain_cpg", "distToTSS_oriented": "distance_to_nearest_tss"}).set_index("cpg_id")
    cpg_annotations = get_cpg_annotations(cpg_coordinates.index, cpg_gene_annotations, args.promoter_length, args.promoter_overlap)
    cpg_sequences = extract_cpg_sequences(cpg_coordinates, args.hg38_genome, args.fasta_index, window_around_cpg=args.window_around_cpg)
    gc_fraction_by_cpg = compute_gc_fraction_by_cpg(cpg_sequences)
    fimo_background_file = args.binding_sites_results_directory / "fimo_background.txt"
    compute_base_frequencies_for_fimo(cpg_sequences, fimo_background_file)

    # Prepare Background.
    cpg_matching_categories = build_cpg_matching_categories(cpg_annotations, gc_fraction_by_cpg, args.n_gc_bins)
    cpgs_for_matching = set(cpg_matching_categories.index)
    background_cpgs_by_mf, n_original_target_cpgs_by_mf, n_target_cpgs_excluded_by_mf = {}, {}, {}
    for mf_symbol, target_cpgs in target_cpgs_by_mf.items():
        n_original_target_cpgs_by_mf[mf_symbol] = len(target_cpgs)
        n_target_cpgs_excluded_by_mf[mf_symbol] = len(target_cpgs - cpgs_for_matching)
        target_cpgs_by_mf[mf_symbol] = target_cpgs & cpgs_for_matching
        background_cpgs_by_mf[mf_symbol] = set(draw_random_matched_background_cpgs(target_cpgs_by_mf[mf_symbol], cpg_matching_categories, bg_to_target_ratio=args.bg_to_target_ratio, bg_sampling_seed=args.bg_sampling_seed, mf_symbol=mf_symbol))

    # FIMO.
    summarized_motif_hit_rates = pd.DataFrame()
    fimo_hits = pd.DataFrame(columns=("motif_id", "sequence_name"))
    fimo_hits_with_cpg_ids = pd.DataFrame(columns=("motif_id", "sequence_name", "cpg_id"))
    if pwm_ids_by_mf:
        cpgs_for_fimo = set() # target and background CpGs for all MFs
        for mf_symbol in pwm_ids_by_mf:
            for cpg_id in target_cpgs_by_mf[mf_symbol]:
                cpgs_for_fimo.add(cpg_id)
            for cpg_id in background_cpgs_by_mf[mf_symbol]:
                cpgs_for_fimo.add(cpg_id)

        cpg_sequences_for_fimo = {}
        for cpg_id in sorted(cpgs_for_fimo):
            if cpg_id in cpg_sequences:
                cpg_sequences_for_fimo[cpg_id] = cpg_sequences[cpg_id]
        cpg_sequences_fasta_file = args.binding_sites_results_directory / "fimo_sequences.fa"
        write_cpg_sequences_as_fasta(cpg_sequences_for_fimo, cpg_sequences_fasta_file)
        pwm_ids = [pwm_id for mf_pwm_ids in pwm_ids_by_mf.values() for pwm_id in mf_pwm_ids]
        selected_motifs_file = args.binding_sites_results_directory / "fimo_motifs.meme"
        write_selected_motifs(args.jaspar_meme, pwm_ids, selected_motifs_file)
        with tempfile.TemporaryDirectory() as fimo_output_dir: 
            fimo_hits = run_fimo(selected_motifs_file, cpg_sequences_fasta_file, fimo_output_dir, fimo_background_file, args.fimo_pvalue, args.fimo_max_stored_scores)
        fimo_hits_with_cpg_ids = map_fimo_hits_to_cpgs(fimo_hits)

        # Compare target and background hit rates per MF.
        motif_hit_cpgs_by_mf = {}
        for mf_symbol, mf_pwm_ids in pwm_ids_by_mf.items():
            mf_hits = fimo_hits_with_cpg_ids.loc[fimo_hits_with_cpg_ids["motif_id"].isin(mf_pwm_ids)]
            motif_hit_cpgs_by_mf[mf_symbol] = set(mf_hits["cpg_id"].dropna())
        motif_target_cpgs_by_mf = {mf_symbol: target_cpgs_by_mf[mf_symbol] for mf_symbol in pwm_ids_by_mf}
        motif_background_cpgs_by_mf = {mf_symbol: background_cpgs_by_mf[mf_symbol] for mf_symbol in pwm_ids_by_mf}
        summarized_motif_hit_rates = summarize_mf_hit_rates(motif_target_cpgs_by_mf, motif_background_cpgs_by_mf, motif_hit_cpgs_by_mf, "FIMO") 
        pwm_ids_text_by_mf = {}
        for mf_symbol, mf_pwm_ids in pwm_ids_by_mf.items():
            pwm_ids_text_by_mf[mf_symbol] = ",".join(mf_pwm_ids)
        summarized_motif_hit_rates["jaspar_ids"] = summarized_motif_hit_rates["MF"].map(pwm_ids_text_by_mf)
        summarized_motif_hit_rates["chipseq_peak_file"] = pd.NA
    fimo_hits.to_csv(args.binding_sites_results_directory / "fimo_hits.tsv", sep="\t", index=False)
    fimo_hits_with_cpg_ids.to_csv(args.binding_sites_results_directory / "fimo_hits_with_cpg_ids.tsv", sep="\t", index=False)

    # ChIP-seq.
    summarized_chipseq_hit_rates = pd.DataFrame()
    if mfs_without_pwm:
        cpgs_for_chipseq = set() # target and background CpGs for all MFs
        for mf_symbol in mfs_without_pwm:
            for cpg_id in target_cpgs_by_mf[mf_symbol]:
                cpgs_for_chipseq.add(cpg_id)
            for cpg_id in background_cpgs_by_mf[mf_symbol]:
                cpgs_for_chipseq.add(cpg_id)

        cpg_coordinates_for_chipseq = cpg_coordinates.loc[cpg_coordinates.index.isin(cpgs_for_chipseq)] # for the peak overlap calculation

        # Compare target and background hit rates per MF.
        chipseq_peak_files_by_mf = {}
        chipseq_hit_cpgs_by_mf = {}
        for mf_symbol in mfs_without_pwm:
            chipseq_peak_file = peak_file_for_mf(mf_symbol, args.chipseq_dir, gene_symbol_aliases = gene_symbol_aliases)
            chipseq_peak_files_by_mf[mf_symbol] = chipseq_peak_file
            if chipseq_peak_file is None:
                chipseq_hit_cpgs_by_mf[mf_symbol] = set()
            else:
                chipseq_hit_cpgs_by_mf[mf_symbol] = cpgs_with_peak_overlap(cpg_coordinates_for_chipseq, chipseq_peak_file, window_around_cpg=args.window_around_cpg)
        chipseq_target_cpgs_by_mf = {mf_symbol: target_cpgs_by_mf[mf_symbol] for mf_symbol in mfs_without_pwm}
        chipseq_background_cpgs_by_mf = {mf_symbol: background_cpgs_by_mf[mf_symbol] for mf_symbol in mfs_without_pwm}
        summarized_chipseq_hit_rates = summarize_mf_hit_rates(chipseq_target_cpgs_by_mf, chipseq_background_cpgs_by_mf, chipseq_hit_cpgs_by_mf, "ChIP_seq")
        chipseq_peak_filenames_by_mf = {}
        mfs_without_chipseq = []
        for mf_symbol, chipseq_peak_file in chipseq_peak_files_by_mf.items():
            if chipseq_peak_file is None:
                chipseq_peak_filenames_by_mf[mf_symbol] = pd.NA
                mfs_without_chipseq.append(mf_symbol)
            else:
                chipseq_peak_filenames_by_mf[mf_symbol] = chipseq_peak_file.name
        summarized_chipseq_hit_rates.loc[summarized_chipseq_hit_rates["MF"].isin(mfs_without_chipseq), "status"] = "no_chipseq"
        summarized_chipseq_hit_rates["chipseq_peak_file"] = summarized_chipseq_hit_rates["MF"].map(chipseq_peak_filenames_by_mf)


    # Save binding site evaluation results.
    motif_sign_test_results = calculate_hit_rate_sign_test(summarized_motif_hit_rates)
    chipseq_sign_test_results = calculate_hit_rate_sign_test(summarized_chipseq_hit_rates)
    binding_site_evaluation_summary = pd.DataFrame([motif_sign_test_results, chipseq_sign_test_results])
    binding_site_evaluation_summary.insert(0, "track", ["FIMO", "ChIP_seq"])
    binding_site_evaluation_results = pd.concat([summarized_motif_hit_rates, summarized_chipseq_hit_rates], ignore_index=True)
    binding_site_evaluation_results = create_binding_site_evaluation_results(binding_site_evaluation_results, n_original_target_cpgs_by_mf, n_target_cpgs_excluded_by_mf)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    binding_site_evaluation_summary.to_csv(args.summary, sep="\t", index=False)
    args.output.parent.mkdir(parents=True, exist_ok=True) 
    binding_site_evaluation_results.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main() 

  
