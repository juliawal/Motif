import subprocess
from pathlib import Path
import pandas as pd


# Find DNA sequences in the windows around CpGs that match DNA binding motifs of the MF proteins (= hits).
def run_fimo(selected_motifs_file, cpg_sequences_fasta_file, fimo_output_dir, fimo_background_file, fimo_pvalue = 1e-4, fimo_max_stored_scores = 1000000):
    fimo_output_dir = Path(fimo_output_dir)
    fimo_output_dir.mkdir(parents=True, exist_ok=True)
    fimo_command = ["fimo", "--thresh", str(fimo_pvalue), "--max-stored-scores", str(fimo_max_stored_scores), "--bfile", str(fimo_background_file), "--no-pgc", "--oc", str(fimo_output_dir), str(selected_motifs_file), str(cpg_sequences_fasta_file)]
    fimo_process = subprocess.run(fimo_command, capture_output = True, text = True, check = True)
    if "Reached max stored scores" in fimo_process.stderr:
        raise RuntimeError("FIMO reached the maximum number of stored hits; >> increase fimo_max_stored_scores.")
    try:
        fimo_hits = pd.read_csv(fimo_output_dir / "fimo.tsv", sep = "\t", comment = "#", usecols = ("motif_id", "sequence_name"))
        return fimo_hits
    except pd.errors.EmptyDataError:
        return pd.DataFrame().reindex(columns = ("motif_id", "sequence_name"))


def map_fimo_hits_to_cpgs(fimo_hits):
    fimo_hits_with_cpg_ids = fimo_hits.copy()
    fimo_hits_with_cpg_ids["cpg_id"] = fimo_hits_with_cpg_ids["sequence_name"]
    return fimo_hits_with_cpg_ids
