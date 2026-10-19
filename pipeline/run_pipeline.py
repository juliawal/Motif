from argparse import ArgumentParser
import subprocess
import sys
from pathlib import Path
import yaml

steps_by_stage = {"A": ("prepare_data",), "B": ("grn_boost_runs",), "C1": ("sex_filter",), "C2": ("reverse_causality_filter",), "D": ("stability", "stability_comparison"), "E": ("stability_selection", "aggregate_runs", "select_top_mfs"), "F0": ("preprocess_evaluation",), "F1": ("binding_site_validation",), "F2": ("ppi_proximity",), "F3": ("gene_set_enrichment",), "F4": ("mediation",), "G": ("evaluation_summary",)}

def get_stage_output_paths(config, stage):
    cohorts = config["cohorts"]
    selection_ranks = range(1, config["n_selected_stable_settings"] + 1)
    save_per_cohort = {"A": "prepared/feature_matrix.tsv", "C1": "filtered/after_sex.tsv", "C2": "filtered/retained_edges.tsv", "D": "stability/grid_stability.tsv"}
    save_per_rank = {"E": "selected_mfs/mfs.tsv", "F0": "evaluation_input/mf_cpg_pairs.tsv", "F1": "evaluation/binding_sites.tsv", "F2": "evaluation/ppi_proximity.tsv", "F3": "evaluation/gene_set_enrichment.tsv", "F4": "evaluation/mediation.tsv"}
    
    if stage == "B": # GRN inference has one output for each cohort and run
        return [f"results/{cohort}/grn/seed_{seed}.tsv" for cohort in cohorts for seed in range(1, config["n_seeds"] + 1)]
    
    if stage in save_per_cohort:
        targets = [f"results/{cohort}/{save_per_cohort[stage]}" for cohort in cohorts]
        if stage == "D": # add a shared file with the most stable settings across all cohorts
            targets.append("results/stability_comparison/selection.tsv")
        return targets
    if stage in save_per_rank:
        return [f"results/selected_settings/rank_{selection_rank}/{cohort}/{save_per_rank[stage]}" for selection_rank in selection_ranks for cohort in cohorts]
    return [f"results/selected_settings/rank_{selection_rank}/evaluation_summary/{filename}" for selection_rank in selection_ranks for filename in ("evaluation_summary.tsv", "enrichment_top_terms.tsv")] # stage G


def snakemake(output_paths, selected_steps, args, config):
    """Run selected steps and repeat interrupted jobs."""
    command = ["snakemake", *output_paths, "--latency-wait", "60", "--rerun-incomplete", "--resources", f"prepare_jobs={config['prepare_max_parallel']}", f"grn_jobs={config['grn_max_parallel']}"]
    if args.executor == "local":
        command.extend(("--cores", args.cores))
    else:
        command.extend(("--executor", "slurm", "--jobs", str(config["slurm_jobs"])))
        if config["slurm_account"]:
            command.extend(("--default-resources", f"slurm_account={config['slurm_account']}"))
    command.extend(("--allowed-rules", *selected_steps))
    subprocess.run(command, check=True)



 
def main():
    parser = ArgumentParser()
    parser.add_argument("stage")
    parser.add_argument("--executor")
    parser.add_argument("--cores")
    args = parser.parse_args()

    config = yaml.safe_load(Path("config.yaml").read_text())
    if args.stage != "all":
        snakemake(get_stage_output_paths(config, args.stage), steps_by_stage[args.stage], args, config)
        return
    selected_steps = [step for steps in steps_by_stage.values() for step in steps]
    downloaded = Path("results/chipseq.done")
    if not downloaded.is_file():
        snakemake(get_stage_output_paths(config, "F0"), selected_steps, args, config)
        subprocess.run([sys.executable, "-m", "scripts.download_chipseq_peaks", "--from-pipeline-results", "results"], check=True)
        downloaded.touch() 
    snakemake(get_stage_output_paths(config, "G"), selected_steps, args, config)


if __name__ == "__main__":
    main()
