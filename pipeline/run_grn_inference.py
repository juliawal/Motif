from argparse import ArgumentParser
import os
import numpy as np
import pandas as pd
from dask.distributed import Client, LocalCluster
from adapted_grnboost import algo


def close_dask_client(client, cluster):
    for close in (client.close, cluster.close):
        try:
            close() # close the dask client
        except Exception:
            pass


def run_grnboost2(feature_matrix, gene_names, cpg_names, output_path, seed, n_threads, grn_limit):
    cluster = LocalCluster(n_workers=n_threads, threads_per_worker=1, processes=True, dashboard_address=None)
    client = Client(cluster)
    try:
        run_edges = algo.grnboost2(expression_data=feature_matrix, tf_names=gene_names, target_genes=cpg_names, seed=seed, limit=grn_limit, client_or_address=client)
        run_edges = run_edges.rename(columns={"TF": "MF_ensembl", "target": "CpG"})
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        run_edges.to_csv(output_path, sep="\t", index=False)
    finally:
        close_dask_client(client, cluster)


def main():
    parser = ArgumentParser()
    parser.add_argument("--feature-matrix", required=True, help="Samples x features (genes and CpGs) .tsv file from prepare_feature_matrix.")
    parser.add_argument("--gene-names", required=True, help="Gene names.")
    parser.add_argument("--cpg-names", required=True, help="CpG or CpG regions.")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--threads", type=int, required=True, help="Dask worker processes.")
    parser.add_argument("--grn-limit", type=int, required=True, help="Keep only top K edges.")
    parser.add_argument("--out-edges",required=True)
    args = parser.parse_args()

    feature_matrix = pd.read_csv(args.feature_matrix, sep="\t", index_col=0).astype(np.float32)
    gene_names = pd.read_csv(args.gene_names, header=None)[0].tolist()
    cpg_names  = pd.read_csv(args.cpg_names, header=None)[0].tolist()

    run_grnboost2(feature_matrix, gene_names, cpg_names, args.out_edges, args.seed, n_threads=args.threads, grn_limit=args.grn_limit)


if __name__ == "__main__":
    main()
