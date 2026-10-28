from pathlib import Path

configfile: "config.yaml"

if config["n_seeds"] % 2:
    raise ValueError("Uneven number of run! Split half stability needs an even number.")

cohorts = config["cohorts"]
run_seeds = list(range(1, config["n_seeds"] + 1))
selection_ranks = list(range(1, config["n_selected_stable_settings"] + 1))
split_cohorts = config["split"]["cohorts"]
use_cpg_regions = config["should_filter_cpgs"] or config["should_aggregate_cpgs"]
chipseq_resources_dir = Path(config["chipseq_dir"])

evaluation_folder = config["evaluation_folder"]
evaluation_cohort_folder = evaluation_folder + "/{cohort}"
evaluation_prepared_folder = "results/{cohort}/prepared"
if evaluation_folder == "results/manual_mfs":
    evaluation_cohort_folder = evaluation_folder
    evaluation_prepared_folder = evaluation_folder + "/prepared"


def get_primary_tumour_type(cohort):
    if cohort in split_cohorts:
        return split_cohorts[cohort]
    primary_tumour_types = config["primary_tumour_type"]
    return primary_tumour_types.get(cohort, primary_tumour_types["default"])


def get_excluded_tumour_type(cohort):
    """Return the other split sample type, or none."""
    if cohort not in split_cohorts:
        return "none"
    for other_cohort, other_tumour_type in split_cohorts.items():
        if other_cohort != cohort:
            return other_tumour_type


def require_cpg_level_targets(job):
    """Stop after stage B if the GRN targets are regional means."""
    if config["should_aggregate_cpgs"]:
        raise ValueError("The stages after GRN inference require individual CpGs. Disable should_aggregate_cpgs or run only stage A snf stage B.")
    return f"results/{job.cohort}/prepared/target_names.tsv"


rule all:
    input:
        # expand() creates a list of file paths for all combinations.
        expand("results/selected_settings/rank_{selection_rank}/{cohort}/selected_mfs/mfs.tsv", selection_rank=selection_ranks, cohort=cohorts),
        expand("results/selected_settings/rank_{selection_rank}/{cohort}/evaluation/binding_sites.tsv", selection_rank=selection_ranks, cohort=cohorts),
        expand("results/selected_settings/rank_{selection_rank}/{cohort}/evaluation/ppi_proximity.tsv", selection_rank=selection_ranks, cohort=cohorts),
        expand("results/selected_settings/rank_{selection_rank}/{cohort}/evaluation/gene_set_enrichment.tsv", selection_rank=selection_ranks, cohort=cohorts),
        expand("results/selected_settings/rank_{selection_rank}/{cohort}/evaluation/mediation.tsv", selection_rank=selection_ranks, cohort=cohorts),
        expand("results/selected_settings/rank_{selection_rank}/evaluation_summary/evaluation_summary.tsv", selection_rank=selection_ranks),
        expand("results/selected_settings/rank_{selection_rank}/evaluation_summary/enrichment_top_terms.tsv", selection_rank=selection_ranks),
        expand("results/{cohort}/stability/grid_stability.tsv", cohort=cohorts),
        "results/stability_comparison/selection.tsv",
        "results/stability_comparison/stability_results_with_enough_mfs.tsv",


# Check if all GRN Boost results are present and run the missing tasks using grn_boost_runs rule.
rule all_grn_results:
    input:
        expand("results/{cohort}/grn/seed_{seed}.tsv", cohort=cohorts, seed=run_seeds),


# Run MF selection and prepare the MF-CpG pairs for evaluation.
rule f0:
    input:
        expand("results/selected_settings/rank_{selection_rank}/{cohort}/evaluation_input/mf_cpg_pairs.tsv", selection_rank=selection_ranks, cohort=cohorts),


rule evaluation:
    input:
        expand("results/selected_settings/rank_{selection_rank}/evaluation_summary/evaluation_summary.tsv", selection_rank=selection_ranks),
        expand("results/selected_settings/rank_{selection_rank}/evaluation_summary/enrichment_top_terms.tsv", selection_rank=selection_ranks),


rule union_evaluation:
    input:
        expand("results/union/{cohort}/evaluation/{filename}", cohort=cohorts, filename=["binding_sites.tsv", "ppi_proximity.tsv", "gene_set_enrichment.tsv", "mediation.tsv"]),


rule manually_evaluation:
    input:
        "results/manual_mfs/evaluation/binding_sites.tsv",
        "results/manual_mfs/evaluation/ppi_proximity.tsv",
        "results/manual_mfs/evaluation/gene_set_enrichment.tsv",


rule prepare_data:
    input:
        cpg_methylation_file = "data/{cohort}/methylation.tsv",
        gene_expression_file = "data/{cohort}/expression.tsv",
        raw_counts = "data/{cohort}/raw_counts.tsv",
        shared_samples = "data/{cohort}/shared_samples.tsv",
        coding_gene_list_file = config["coding_gene_list_file"],
        cpg_annotations_file = config["cpg_annotations_file"],
        cpg_filter_genes_file = config["cpg_filter_genes_file"] if config["should_filter_cpgs"] else [],
        cpg_probe_coordinates_file = config["cpg_probe_coordinates_file"] if use_cpg_regions else [],
        gene_regions_gtf_file = config["gene_regions_gtf_file"] if use_cpg_regions else [],
    output:
        feature_matrix = "results/{cohort}/prepared/feature_matrix.tsv",
        gene_names = "results/{cohort}/prepared/regulator_names.tsv",
        cpg_names = "results/{cohort}/prepared/target_names.tsv",
        cpg_region_map = "results/{cohort}/prepared/cpg_region_map.tsv",
        samples = "results/{cohort}/prepared/samples.tsv",
        counts = "results/{cohort}/prepared/counts.tsv",
    params:
        only_coding_genes = config["only_coding_genes"],
        primary_tumour_type = lambda job: get_primary_tumour_type(job.cohort),
        exclude_patients_of_type = lambda job: get_excluded_tumour_type(job.cohort),
        n_hvgs = config["n_hvgs"],
        should_filter_cpgs = config["should_filter_cpgs"],
        should_aggregate_cpgs = config["should_aggregate_cpgs"],
        cpg_keep_promoters = config["cpg_keep_promoters"],
        cpg_keep_bodies = config["cpg_keep_bodies"],
        promoter_length = config["promoter_length"],
        promoter_overlap = config["promoter_overlap"],
        sample_max_missing_share = config["sample_max_missing_share"],
        cpg_filter_genes_file = config["cpg_filter_genes_file"],
        cpg_probe_coordinates_file = config["cpg_probe_coordinates_file"],
        gene_regions_gtf_file = config["gene_regions_gtf_file"],
    threads: config["prepare_threads"]
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = config["prepare_runtime_minutes"],
        mem_mb = config["prepare_mem_mb"],
        prepare_jobs = 1,
    shell:
        "python -m pipeline.prepare_feature_matrix "
        "--cpg-methylation-file {input.cpg_methylation_file} "
        "--gene-expression-file {input.gene_expression_file} "
        "--raw-counts {input.raw_counts} "
        "--cpg-annotations-file {input.cpg_annotations_file} "
        "--shared-samples {input.shared_samples} "
        "--coding-gene-list-file {input.coding_gene_list_file} "
        "--out-feature-matrix {output.feature_matrix} "
        "--out-gene-names {output.gene_names} "
        "--out-cpg-names {output.cpg_names} "
        "--out-cpg-gene-region-map {output.cpg_region_map} "
        "--out-samples {output.samples} "
        "--out-counts {output.counts} "
        "--only-coding-genes {params.only_coding_genes} "
        "--primary-tumour-type {params.primary_tumour_type} "
        "--exclude-patients-of-type {params.exclude_patients_of_type} "
        "--n-hvgs {params.n_hvgs} "
        "--sample-max-missing-share {params.sample_max_missing_share} "
        "--should-filter-cpgs {params.should_filter_cpgs} "
        "--should-aggregate-cpgs {params.should_aggregate_cpgs} "
        "--cpg-filter-genes-file {params.cpg_filter_genes_file} "
        "--cpg-probe-coordinates-file {params.cpg_probe_coordinates_file} "
        "--gene-regions-gtf-file {params.gene_regions_gtf_file} "
        "--cpg-keep-promoters {params.cpg_keep_promoters} "
        "--cpg-keep-bodies {params.cpg_keep_bodies} "
        "--promoter-length {params.promoter_length} "
        "--promoter-overlap {params.promoter_overlap}"


rule grn_boost_runs:
    input: 
        feature_matrix = "results/{cohort}/prepared/feature_matrix.tsv",
        gene_names = "results/{cohort}/prepared/regulator_names.tsv",
        cpg_names = "results/{cohort}/prepared/target_names.tsv",
    output:
        edges = "results/{cohort}/grn/seed_{seed}.tsv",
    params:
        seed = "{seed}",
        grn_limit = config["grn_limit"],
    threads: config["grn_threads"]
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = config["grn_runtime_minutes"],
        mem_mb = config["grn_mem_mb"] or 0, # 0 leaves grn_mem_mb to Slurm
        grn_jobs = 1,
    shell:
        r"""
        scratch=$(mktemp -d "${{TMPDIR:-/tmp}}/motif_grn.XXXXXX")
        partial={output.edges:q}.partial
        trap 'rm -rf -- "$scratch"; rm -f -- "$partial"' EXIT
        cp {input:q} "$scratch/"
        mkdir "$scratch/dask"
        export DASK_TEMPORARY_DIRECTORY="$scratch/dask"
        export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
        export PYTHONHASHSEED=0

        python -m pipeline.run_grn_inference \
            --feature-matrix "$scratch/feature_matrix.tsv" \
            --gene-names "$scratch/regulator_names.tsv" \
            --cpg-names "$scratch/target_names.tsv" \
            --seed {params.seed} --threads {threads} --grn-limit {params.grn_limit} \
            --out-edges "$partial"
        mv -- "$partial" {output.edges:q}
        """


# Apply the sex filter to the union of all GRN runs.
rule sex_filter:
    input:
        edges_per_run_files = expand("results/{{cohort}}/grn/seed_{seed}.tsv", seed=run_seeds),
        target_contract = require_cpg_level_targets,
        samples = "results/{cohort}/prepared/samples.tsv",
        counts = "results/{cohort}/prepared/counts.tsv",
        cpg_methylation_file = "data/{cohort}/methylation.tsv",
        metadata = "data/{cohort}/metadata.tsv",
    output:
        pairs_after_sex_filter = "results/{cohort}/filtered/after_sex.tsv",
        gene_results = "results/{cohort}/filtered/gene_sex_test.tsv",
        cpg_results = "results/{cohort}/filtered/cpg_sex_test.tsv",
    params:
        fdr_alpha = config["fdr_alpha"],
        m_value_boundary_shift = config["m_value_boundary_shift"],
    threads: 4
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = 240,
        mem_mb = 48000,
    shell:
        "python -m pipeline.apply_filters sex-filter "
        "--edges-per-run-files {input.edges_per_run_files} "
        "--counts {input.counts} "
        "--cpg-methylation {input.cpg_methylation_file} "
        "--metadata {input.metadata} "
        "--samples {input.samples} "
        "--fdr-alpha {params.fdr_alpha} "
        "--m-value-boundary-shift {params.m_value_boundary_shift} "
        "--n-cpus {threads} "
        "--output {output.pairs_after_sex_filter}"


rule reverse_causality_filter:
    input:
        pairs_after_sex_filter = "results/{cohort}/filtered/after_sex.tsv",
        samples = "results/{cohort}/prepared/samples.tsv",
        gene_expression_file = "data/{cohort}/expression.tsv",
        cpg_methylation_file = "data/{cohort}/methylation.tsv",
        cpg_coordinates_file = config["cpg_coordinates_file"],
        gene_tss_annotations_file = config["gene_tss_annotations_file"],
    output:
        pairs_after_filters = "results/{cohort}/filtered/retained_edges.tsv",
    params:
        reverse_tss_window = config["reverse_tss_window"],
    threads: 2
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = 180,
        mem_mb = 32000,
    shell:
        "python -m pipeline.apply_filters reverse-causality-filter "
        "--pairs-after-sex-filter {input.pairs_after_sex_filter} "
        "--gene-expression {input.gene_expression_file} "
        "--cpg-methylation {input.cpg_methylation_file} "
        "--cpg-coordinates-file {input.cpg_coordinates_file} "
        "--gene-tss-annotations-file {input.gene_tss_annotations_file} "
        "--samples {input.samples} "
        "--reverse-tss-window {params.reverse_tss_window} "
        "--output {output.pairs_after_filters}"


rule aggregate_runs:
    input:
        edges_per_run_files = expand("results/{{cohort}}/grn/seed_{seed}.tsv", seed=run_seeds),
        pairs_after_filters = "results/{cohort}/filtered/retained_edges.tsv",
        selected_setting_file = "results/selected_settings/rank_{selection_rank}/selection_applied.tsv",
    output:
        aggregated_edges = "results/selected_settings/rank_{selection_rank}/{cohort}/aggregated/edges.tsv",
    params:
        run_seeds = " ".join(str(run_seed) for run_seed in run_seeds),
    threads: 2
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = 60,
        mem_mb = 16000,
    shell:
        "python -m pipeline.cut_and_aggregate_runs "
        "--edges-per-run-files {input.edges_per_run_files} "
        "--run-seeds {params.run_seeds} "
        "--pairs-after-filters {input.pairs_after_filters} "
        "--selected-setting-file {input.selected_setting_file} "
        "--output {output.aggregated_edges}"


# Select MFs using the score from the chosen stability setting.
rule select_top_mfs:
    input:
        aggregated_edges = "results/selected_settings/rank_{selection_rank}/{cohort}/aggregated/edges.tsv",
        selected_setting_file = "results/selected_settings/rank_{selection_rank}/selection_applied.tsv",
    output:
        selected_mfs = "results/selected_settings/rank_{selection_rank}/{cohort}/selected_mfs/mfs.tsv",
    params:
        n_mfs = config["n_mfs"],
        run_seeds = " ".join(str(run_seed) for run_seed in run_seeds),
    threads: 1
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = 30,
        mem_mb = 4000,
    shell:
        "python -m pipeline.select_top_mfs "
        "--aggregated-edges {input.aggregated_edges} "
        "--selected-setting-file {input.selected_setting_file} "
        "--n-mfs {params.n_mfs} "
        "--run-seeds {params.run_seeds} "
        "--output {output.selected_mfs}"


rule stability:
    input:
        edges_per_run_files = expand("results/{{cohort}}/grn/seed_{seed}.tsv", seed=run_seeds),
        gene_names = "results/{cohort}/prepared/regulator_names.tsv",
        pairs_after_filters = "results/{cohort}/filtered/retained_edges.tsv",
    output:
        cohort_stability_results = "results/{cohort}/stability/grid_stability.tsv",
    params:
        cohort = "{cohort}",
        run_seeds = " ".join(str(run_seed) for run_seed in run_seeds),
        n_mfs = config["n_mfs"],
        k_edges = " ".join(str(k_edges) for k_edges in config["stability_k_edges"]),
        n_splits = config["stability_max_splits"],
        split_seed = config["stability_split_seed"],
    threads: 4
    resources:
        slurm_partition = config["slurm_partition"],
        mem_mb  = 16000,
        runtime = 1440,
    shell:
        "python -m pipeline.calculate_stability "
        "--edges-per-run-files {input.edges_per_run_files} "
        "--run-seeds {params.run_seeds} "
        "--gene-names {input.gene_names} "
        "--pairs-after-filters {input.pairs_after_filters} "
        "--cohort {params.cohort} "
        "--n-mfs {params.n_mfs} "
        "--k-edges {params.k_edges} "
        "--n-splits {params.n_splits} "
        "--split-seed {params.split_seed} "
        "--output {output.cohort_stability_results}"


# Rank settings by their stability across cohorts.
rule stability_comparison:
    input:
        grid_stability_files = expand("results/{cohort}/stability/grid_stability.tsv", cohort=cohorts),
    output:
        best_settings_file = "results/stability_comparison/selection.tsv",
        stability_results_with_enough_mfs = "results/stability_comparison/stability_results_with_enough_mfs.tsv",
        ranked_settings_with_stability = "results/stability_comparison/ranking.tsv",
        stability_matrix = "results/stability_comparison/kuncheva_by_cohort.tsv",
        cohorts_results_for_best_settings = "results/stability_comparison/selected_stability.tsv",
    params:
        output_dir = "results/stability_comparison",
        n_selected_top_settings = config["n_selected_stable_settings"],
    threads: 1
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = 60,
        mem_mb = 8000,
    shell:
        "python -m pipeline.select_most_stable_setting "
        "--grid-stability-files {input.grid_stability_files} "
        "--n-selected-top-settings {params.n_selected_top_settings} "
        "--output-dir {params.output_dir}"


# Save the selected setting for this rank.
rule stability_selection:
    input:
        best_settings_file = "results/stability_comparison/selection.tsv",
    output:
        selected_setting_file = "results/selected_settings/rank_{selection_rank}/selection_applied.tsv",
    params:
        selection_rank = "{selection_rank}",
    threads: 1
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = 5,
        mem_mb = 1000,
    shell:
        "python -m pipeline.apply_selection "
        "--best-settings-file {input.best_settings_file} "
        "--selection-rank {params.selection_rank} "
        "--output {output.selected_setting_file}"


# Map selected Ensembl IDs to symbols and keep their CpG pairs.
rule preprocess_evaluation:
    input:
        aggregated_edges_file = "results/selected_settings/rank_{selection_rank}/{cohort}/aggregated/edges.tsv",
        selected_mfs_file = "results/selected_settings/rank_{selection_rank}/{cohort}/selected_mfs/mfs.tsv",
        gene_annotation_file = config["gene_annotation_file"],
    output:
        evaluation_pairs = "results/selected_settings/rank_{selection_rank}/{cohort}/evaluation_input/mf_cpg_pairs.tsv",
        mapping = "results/selected_settings/rank_{selection_rank}/{cohort}/evaluation_input/mf_mapping.tsv",
    threads: 1
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = 30,
        mem_mb = 4000,
    shell:
        "python -m pipeline.prepare_evaluation "
        "--aggregated-edges-file {input.aggregated_edges_file} "
        "--selected-mfs-file {input.selected_mfs_file} "
        "--gene-annotation-file {input.gene_annotation_file} "
        "--output {output.evaluation_pairs} "
        "--output-mf-mapping {output.mapping}"


rule prepare_chipseq_resources:
    input:
        pairs = expand("results/selected_settings/rank_{selection_rank}/{cohort}/evaluation_input/mf_cpg_pairs.tsv", selection_rank=selection_ranks, cohort=cohorts),
        remap_file = config["remap_file"],
        gene_annotation_file = config["gene_annotation_file"],
        jaspar_meme = config["jaspar_meme"],
        jaspar_sqlite = config["jaspar_sqlite"],
    output:
        resources = directory(str(chipseq_resources_dir)),
    threads: 4
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = 120,
        mem_mb = 24000,
    shell:
        "python -m scripts.download_and_prepare_chip_seq_peak_files prepare "
        "--config config.yaml "
        "--pairs {input.pairs:q} --output-dir {output.resources:q}"


rule binding_site_validation:
    input:
        mf_cpg_pairs = evaluation_cohort_folder + "/evaluation_input/mf_cpg_pairs.tsv",
        cpg_annotations_file = config["cpg_annotations_file"],
        cpg_names = evaluation_prepared_folder + "/target_names.tsv",
        cpg_gene_annotations_file = config["cpg_gene_annotations_file"],
        gene_annotation_file = config["gene_annotation_file"],
        hg38_genome = config["genome_fasta"],
        fasta_index = config["genome_fasta_index"],
        jaspar_meme = config["jaspar_meme"],
        jaspar_sqlite = config["jaspar_sqlite"],
        chipseq_resources = str(chipseq_resources_dir),
    output:
        results = evaluation_cohort_folder + "/evaluation/binding_sites_per_mf.tsv",
        summary = evaluation_cohort_folder + "/evaluation/binding_sites.tsv",
        binding_sites_results_directory = directory(evaluation_cohort_folder + "/evaluation/binding_evidence"),
    params:
        chipseq_peak_files = str(chipseq_resources_dir / "peak_files.json"),
        window_around_cpg = config["window_around_cpg"],
        bg_to_target_ratio = config["bg_to_target_ratio"],
        n_gc_bins = config["n_gc_bins"],
        bg_sampling_seed = config["binding_background_seed"],
        fimo_pvalue = config["fimo_pvalue"],
        fimo_max_stored_scores = config["fimo_max_stored_scores"],
        promoter_length = config["promoter_length"],
        promoter_overlap = config["promoter_overlap"],
    threads: 8
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = 1080,
        mem_mb = 32000,
    shell:
        "python -m pipeline.evaluate_binding_sites "
        "--mf-cpg-pairs {input.mf_cpg_pairs} "
        "--cpg-annotations-file {input.cpg_annotations_file} "
        "--cpg-names {input.cpg_names} "
        "--cpg-gene-annotations-file {input.cpg_gene_annotations_file} "
        "--gene-annotation-file {input.gene_annotation_file} "
        "--hg38-genome {input.hg38_genome} "
        "--fasta-index {input.fasta_index} "
        "--chipseq-peak-files {params.chipseq_peak_files} "
        "--jaspar-meme {input.jaspar_meme} "
        "--jaspar-sqlite {input.jaspar_sqlite} "
        "--window-around-cpg {params.window_around_cpg} "
        "--bg-to-target-ratio {params.bg_to_target_ratio} "
        "--n-gc-bins {params.n_gc_bins} "
        "--bg-sampling-seed {params.bg_sampling_seed} "
        "--fimo-pvalue {params.fimo_pvalue} "
        "--fimo-max-stored-scores {params.fimo_max_stored_scores} "
        "--promoter-length {params.promoter_length} "
        "--promoter-overlap {params.promoter_overlap} "
        "--output {output.results} "
        "--summary {output.summary} "
        "--binding-sites-results-directory {output.binding_sites_results_directory}"


rule ppi_proximity:
    input:
        mf_cpg_pairs = evaluation_cohort_folder + "/evaluation_input/mf_cpg_pairs.tsv",
        hippie_network_file = config["hippie_network_file"],
        gene_annotation_file = config["gene_annotation_file"],
        gene_ids = evaluation_prepared_folder + "/regulator_names.tsv",
    output:
        results = evaluation_cohort_folder + "/evaluation/ppi_proximity.tsv",
    params:
        n_random_sets = config["n_random_sets"],
        degree_bin_min_size = config["degree_bin_min_size"],
        pagerank_alpha = config["proximity_pagerank_alpha"],
        random_seed = config["proximity_random_seed"],
        hippie_release = config["hippie_release"],
    threads: 4
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = 120,
        mem_mb = 16000,
    shell:
        "python -m pipeline.evaluate_ppi_network_proximity "
        "--mf-cpg-pairs {input.mf_cpg_pairs} "
        "--hippie-network-file {input.hippie_network_file} "
        "--hippie-release {params.hippie_release} "
        "--gene-annotation-file {input.gene_annotation_file} "
        "--gene-ids {input.gene_ids} "
        "--n-random-sets {params.n_random_sets} "
        "--degree-bin-min-size {params.degree_bin_min_size} "
        "--pagerank-alpha {params.pagerank_alpha} "
        "--random-seed {params.random_seed} "
        "--output {output.results}"


rule gene_set_enrichment:
    input:
        mf_cpg_pairs = evaluation_cohort_folder + "/evaluation_input/mf_cpg_pairs.tsv",
        background_gene_ids = evaluation_prepared_folder + "/regulator_names.tsv",
        gene_annotation_file = config["gene_annotation_file"],
        gene_set_files = [
            f"{config['gene_set_dir']}/{gene_set_library}.gmt" for gene_set_library in config["gene_set_enrichment_libraries"]],
    output:
        results = evaluation_cohort_folder + "/evaluation/gene_set_enrichment.tsv",
    params:
        min_overlap = config["min_overlap"],
        gene_set_dir = config["gene_set_dir"],
        gene_set_libraries = " ".join(config["gene_set_enrichment_libraries"]),
    threads: 2
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = 60,
        mem_mb = 8000,
    shell:
        "python -m pipeline.evaluate_gene_enrichment "
        "--mf-cpg-pairs {input.mf_cpg_pairs} "
        "--background-gene-ids {input.background_gene_ids} "
        "--gene-annotation-file {input.gene_annotation_file} "
        "--gene-set-dir {params.gene_set_dir} "
        "--gene-set-libraries {params.gene_set_libraries} "
        "--min-overlap {params.min_overlap} "
        "--output {output.results}"


rule mediation:
    input:
        mf_cpg_pairs = evaluation_folder + "/{cohort}/evaluation_input/mf_cpg_pairs.tsv",
        prepared_samples = "results/{cohort}/prepared/samples.tsv",
        metadata = "data/{cohort}/age_and_sex.tsv",
        gene_expression_file = "data/{cohort}/expression.tsv",
        cpg_methylation_file = "data/{cohort}/methylation.tsv",
        cpg_coordinates_file = config["cpg_coordinates_file"],
        gene_annotation_file = config["gene_annotation_file"],
    output:
        triplets = evaluation_folder + "/{cohort}/evaluation/mediation_triplets.tsv",
        summary = evaluation_folder + "/{cohort}/evaluation/mediation.tsv",
    params:
        cohort = "{cohort}",
        max_cpg_tss_distance_bp = config["max_cpg_tss_distance_bp"],
        min_samples = config["mediation_min_samples"],
        fdr_alpha = config["mediation_fdr_alpha"],
        m_value_boundary_shift = config["m_value_boundary_shift"],
    threads: 1
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = 1440,
        mem_mb = 32000,
    shell:
        "python -m pipeline.evaluate_mediation "
        "--mf-cpg-pairs {input.mf_cpg_pairs} "
        "--prepared-samples {input.prepared_samples} "
        "--metadata {input.metadata} "
        "--gene-expression-file {input.gene_expression_file} "
        "--cpg-methylation-file {input.cpg_methylation_file} "
        "--cpg-coordinates-file {input.cpg_coordinates_file} "
        "--gene-annotation-file {input.gene_annotation_file} "
        "--cohort {params.cohort} "
        "--output-triplets {output.triplets} "
        "--summary {output.summary} "
        "--max-cpg-tss-distance-bp {params.max_cpg_tss_distance_bp} "
        "--min-samples {params.min_samples} "
        "--fdr-alpha {params.fdr_alpha} "
        "--m-value-boundary-shift {params.m_value_boundary_shift}"


rule evaluation_summary:
    input:
        selected_setting_file = "results/selected_settings/rank_{selection_rank}/selection_applied.tsv",
        selected_mfs = expand("results/selected_settings/rank_{{selection_rank}}/{cohort}/selected_mfs/mfs.tsv", cohort=cohorts),
        mf_cpg_pairs = expand("results/selected_settings/rank_{{selection_rank}}/{cohort}/evaluation_input/mf_cpg_pairs.tsv", cohort=cohorts),
        binding = expand("results/selected_settings/rank_{{selection_rank}}/{cohort}/evaluation/binding_sites.tsv", cohort=cohorts),
        proximity = expand("results/selected_settings/rank_{{selection_rank}}/{cohort}/evaluation/ppi_proximity.tsv", cohort=cohorts),
        enrichment = expand("results/selected_settings/rank_{{selection_rank}}/{cohort}/evaluation/gene_set_enrichment.tsv", cohort=cohorts),
        mediation = expand("results/selected_settings/rank_{{selection_rank}}/{cohort}/evaluation/mediation.tsv", cohort=cohorts),
    output:
        summary = "results/selected_settings/rank_{selection_rank}/evaluation_summary/evaluation_summary.tsv",
        detail  = "results/selected_settings/rank_{selection_rank}/evaluation_summary/evaluation_by_cohort.tsv",
        enrichment_top_terms = "results/selected_settings/rank_{selection_rank}/evaluation_summary/enrichment_top_terms.tsv",
    params:
        selected_setting_folder = "results/selected_settings/rank_{selection_rank}",
        output_dir = "results/selected_settings/rank_{selection_rank}/evaluation_summary",
        cohorts = " ".join(cohorts),
        min_overlap = config["min_overlap"],
        n_enrichment_top_terms = config["n_enrichment_top_terms"],
    threads: 1
    resources:
        slurm_partition = config["slurm_partition"],
        runtime = 30,
        mem_mb = 4000,
    shell:
        "python -m pipeline.summarize_evaluation "
        "--selected-setting-folder {params.selected_setting_folder} "
        "--cohorts {params.cohorts} "
        "--min-overlap {params.min_overlap} "
        "--n-enrichment-top-terms {params.n_enrichment_top_terms} "
        "--selected-setting-file {input.selected_setting_file} "
        "--output-dir {params.output_dir}" 
