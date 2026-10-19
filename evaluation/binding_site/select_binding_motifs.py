from pathlib import Path
from pyjaspar import jaspardb
 


def get_pwm_ids_in_meme_file(jaspar_meme_file):
    pwm_ids = set()
    with open(jaspar_meme_file) as jaspar_meme_file:
        for line in jaspar_meme_file:
            if line.startswith("MOTIF "):
                pwm_ids.add(line.split()[1])
    return pwm_ids
 

def find_motifs_for_mfs_in_jaspar_database(mf_symbols, jaspar_meme_file, jaspar_sqlite_file, gene_symbol_aliases): 
    pwm_ids_in_meme_file = get_pwm_ids_in_meme_file(jaspar_meme_file)
    jaspar_database = jaspardb(sqlite_db_path=jaspar_sqlite_file)
    pwm_ids_by_mf = {}
    mfs_without_pwm = []
    for mf_symbol in mf_symbols:
        mf_pwm_ids = []
        for symbol in [mf_symbol, *gene_symbol_aliases.get(mf_symbol, [])]:
            if not mf_pwm_ids:
                motifs_for_symbol = jaspar_database.fetch_motifs_by_name(symbol.replace("'", "''")) # sometimes gene names contain an apostrophe, like 3'HEXO. In the database query, 3'HEXO would make the database think the name ends after 3. For that we write 3''HEXO instead. In the database query '3''HEXO' is read as 3'HEXO
                for motif in motifs_for_symbol:
                    if motif.matrix_id in pwm_ids_in_meme_file and "9606" in motif.species: # 9606 is species ID for Homo sapiens
                        mf_pwm_ids.append(motif.matrix_id)
        if mf_pwm_ids:
            pwm_ids_by_mf[mf_symbol] = mf_pwm_ids
        else:
            mfs_without_pwm.append(mf_symbol)
    return pwm_ids_by_mf, mfs_without_pwm

 

def write_selected_motifs(jaspar_meme_file, pwm_ids, selected_motifs_meme_file):
    selected_pwm_ids = set(pwm_ids)
    current_motif_pwm_id_is_in_selected_ids = True 
    with open(jaspar_meme_file) as jaspar_meme_file:
        with open(selected_motifs_meme_file, "w", newline = "") as selected_motifs_meme_file:
            for meme_line in jaspar_meme_file:
                if meme_line.startswith("MOTIF "):
                    pwm_id = meme_line.split()[1]
                    current_motif_pwm_id_is_in_selected_ids = pwm_id in selected_pwm_ids
                if current_motif_pwm_id_is_in_selected_ids:
                    selected_motifs_meme_file.write(meme_line)
 