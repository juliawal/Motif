import pandas as pd


def select_pairs_with_cpg_coordinates(mf_cpg_pairs_with_weight, cpg_coordinates):
    mf_cpg_pairs = mf_cpg_pairs_with_weight[["MF", "MF_ensembl", "CpG"]]
    cpg_positions = cpg_coordinates[["chrom", "pos"]]
    cpg_positions = cpg_positions.rename(columns={"chrom": "cpg_chrom", "pos": "cpg_pos"})
    pairs_with_cpg_coordinates = mf_cpg_pairs.join(cpg_positions, on="CpG")
    pairs_with_cpg_coordinates = pairs_with_cpg_coordinates.dropna(subset=["cpg_chrom", "cpg_pos"]) # remove pairs when the CpG has no value for chromosome or position
    return pairs_with_cpg_coordinates


def select_triplets(pairs_with_cpg_coordinates, gene_annotation, expression_gene_ids, max_cpg_tss_distance_bp):
    """Find possible G' genes for each MF–CpG pair. Use protein-coding genes with expression data with it's TSS on the same chromosome and close enough (max_cpg_tss_distance_bp) to the CpG.
    Return one triplet for each gene G'."""
    # The pinned gene annotation already has unique, unversioned IDs and valid TSS positions.
    is_protein_coding = gene_annotation["gene_type"] == "protein_coding"
    has_expression_data = gene_annotation["ensembl_id"].isin(expression_gene_ids)
    g_prime_candidates = gene_annotation.loc[is_protein_coding & has_expression_data, ["ensembl_id", "symbol", "chromosome", "tss"]].copy()
    mediation_triplet_rows = []
    for mf_cpg_pair in pairs_with_cpg_coordinates.itertuples():
        g_primes_distances_to_cpg = g_prime_candidates["tss"] - mf_cpg_pair.cpg_pos
        g_primes_on_same_chromosome = g_prime_candidates["chromosome"] == mf_cpg_pair.cpg_chrom
        g_primes_within_max_cpg_tss_distance = g_primes_distances_to_cpg.abs() <= max_cpg_tss_distance_bp
        g_primes_near_cpg = g_prime_candidates.loc[g_primes_on_same_chromosome & g_primes_within_max_cpg_tss_distance].copy()

        # Remove the MF from g_primes_near_cpg.
        is_mf = g_primes_near_cpg["ensembl_id"] == mf_cpg_pair.MF_ensembl
        mf_row_indices = g_primes_near_cpg.index[is_mf]
        g_primes_near_cpg = g_primes_near_cpg.drop(index=mf_row_indices)
        
        g_primes_near_cpg["distance"] = g_primes_near_cpg["tss"] - mf_cpg_pair.cpg_pos
        g_primes_near_cpg["absolute_distance"] = g_primes_near_cpg["distance"].abs()
        g_primes_near_cpg = g_primes_near_cpg.sort_values(["absolute_distance", "ensembl_id"])
        for gene in g_primes_near_cpg.itertuples(index=False):
            mediation_triplet_rows.append([mf_cpg_pair.MF, mf_cpg_pair.MF_ensembl, mf_cpg_pair.CpG, gene.symbol, gene.ensembl_id, gene.distance])
    mediation_triplets = pd.DataFrame.from_records(mediation_triplet_rows, columns=["MF", "MF_ensembl", "CpG", "G_prime", "G_prime_ensembl", "G_prime_distance_bp"])
    return mediation_triplets
