"""
TCHRONET.py

This script constructs a correlation-based network from a feature matrix using 
Spearman correlations. Instead of stopping at the first connected network,
it now saves edges in multiple correlation threshold ranges (e.g., 0.4–0.5, 0.5–0.6, etc.).

Main workflow:
1. Load input matrix and preprocess for Spearman correlation
2. Compute pairwise correlations in parallel
3. Save separate HDF5 files for edges within defined correlation bins
"""

import argparse
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
import deepgraph as dg
import shutil


class ei_args:
    """Container class for edge computation arguments."""
    def __init__(self, *, in_graph, in_pos_array, in_matrix, in_n_samples, in_tmp_dir, in_th):
        self.graph = in_graph
        self.pos_array = in_pos_array
        self.data_matrix = in_matrix
        self.n_samples = in_n_samples
        self.tmp_dir = in_tmp_dir
        self.th = in_th


def reset_dir(path):
    """Deletes a directory if it exists, then recreates it."""
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


def spearman_corr(index_s, index_t, data_matrix, n_samples):
    """Compute pairwise Spearman correlations efficiently."""
    features_s = data_matrix[index_s].astype(np.float32)
    features_t = data_matrix[index_t].astype(np.float32)
    corr = np.einsum('ij,ij->i', features_s, features_t) / n_samples
    np.round(corr, 3, out=corr)
    return corr


def run_create_ei(args):
    """Wrapper for multiprocessing."""
    return create_ei(*args)


def create_ei(i, ei_args):
    """Compute correlations and save above-threshold edges to temporary pickles."""
    graph = dg.DeepGraph(ei_args.graph.v.copy(deep=True))
    pos_array = ei_args.pos_array
    step_size = 1e5
    chunk_size = 10000
    th = ei_args.th

    from_pos = pos_array[i]
    to_pos = pos_array[i + 1]

    def local_spearman_corr(index_s, index_t):
        """
        Local wrapper for spearman_corr that captures ei_args from closure.
        """
        corr = spearman_corr(index_s, index_t, ei_args.data_matrix, ei_args.n_samples)
        return corr

    edge_chunks = range(from_pos, to_pos, chunk_size)
    for chunk_start in edge_chunks:
        graph.create_edges(
            connectors=local_spearman_corr,
            step_size=step_size,
            from_pos=chunk_start,
            to_pos=min(chunk_start + chunk_size, to_pos)
        )

        graph.e = graph.e[graph.e['corr'] > th]

        graph.e['row_name_s'] = graph.e.index.get_level_values('s').map(
            graph.v.set_index('index')['row_name']
        )
        graph.e['row_name_t'] = graph.e.index.get_level_values('t').map(
            graph.v.set_index('index')['row_name']
        )

        graph.e = graph.e.reset_index().drop(['s', 't'], axis=1)
        graph.e = graph.e.set_index(['row_name_s', 'row_name_t'])

        tmp_save_pickle = os.path.join(ei_args.tmp_dir, f'{i}_{chunk_start}.pickle')
        graph.e.to_pickle(tmp_save_pickle)
        graph.e.drop(graph.e.index, inplace=True)


def main():
    parser = argparse.ArgumentParser(prog="TCHRONET")

    parser.add_argument("-m", "--matrix", type=str, required=True,
                        help="Path to the input matrix (TSV format, features as rows)")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Prefix for output HDF5 files (one per threshold range)")
    parser.add_argument("-s", "--stepsize", type=int, required=False, default=1e5,
                        help="Stepsize for RAM parameters (DeepGraph parameter)")
    parser.add_argument("-@", "--threads", type=int, required=False, default=8,
                        help="Number of threads to use for parallel processing")
    parser.add_argument("-t", "--tempdir", type=str, required=False,
                        help="Temporary directory to use. Will create a '/tmp' directory inside it")
    parser.add_argument("--min_th", type=float, default=0.4,
                        help="Minimum correlation threshold (start of first interval)")
    parser.add_argument("--max_th", type=float, default=1.0,
                        help="Maximum correlation threshold (end of last interval)")
    parser.add_argument("--step", type=float, default=0.1,
                        help="Size of each correlation interval (e.g., 0.1 for 0.4–0.5, 0.5–0.6, etc.)")

    args = parser.parse_args()

    # Setup temporary directory
    tempfolder_name = "tchronet_tmp"
    if args.tempdir is not None:
        abs_args_tempdir = os.path.abspath(args.tempdir)
        final_tempdir = os.path.join(abs_args_tempdir, tempfolder_name)
    else:
        final_tempdir = os.path.join(os.path.abspath('.'), tempfolder_name)
    reset_dir(final_tempdir)
    print(f"Using: {final_tempdir} as temporary folder")

    # Load data
    print("Loading input matrix...")
    data_matrix_df = pd.read_csv(args.matrix, sep="\t", index_col=0)
    row_names = list(data_matrix_df.index)
    data_matrix = data_matrix_df.to_numpy()
    n_features = data_matrix.shape[0]
    n_samples = data_matrix.shape[1]

    # Preprocess for Spearman
    print("Preprocessing matrix for Spearman correlation...")
    data_matrix = data_matrix.argsort(axis=1).argsort(axis=1)
    data_matrix = (data_matrix - data_matrix.mean(axis=1, keepdims=True)) / data_matrix.std(axis=1, keepdims=True)

    # Prepare DeepGraph vertex DataFrame
    v = pd.DataFrame({'index': range(data_matrix.shape[0]), 'row_name': row_names})

    # Divide edge space
    if args.threads == 1:
        pos_array = np.array([0, n_features * (n_features - 1) // 2], dtype=int)
        indices = np.array([0])
    else:
        pos_array = np.array(np.linspace(0, n_features * (n_features - 1) // 2, args.threads), dtype=int)
        indices = np.arange(0, args.threads - 1)

    # Create shared parameters
    ei_args_obj = ei_args(
        in_graph=dg.DeepGraph(v),
        in_pos_array=pos_array,
        in_matrix=data_matrix,
        in_n_samples=n_samples,
        in_tmp_dir=final_tempdir,
        in_th=args.min_th
    )

    ##############################################
    ##### Modified: save per threshold interval
    ##############################################
    thresholds = np.arange(args.min_th, args.max_th, args.step)

    for low_th in thresholds:
        high_th = round(low_th + args.step, 3)
        ei_args_obj.th = low_th

        print(f"\n>>> Processing threshold range: {low_th:.2f} - {high_th:.2f}")

        # Parallel correlation computation
        ei_args_list = [(i, ei_args_obj) for i in indices]
        with Pool(processes=args.threads) as pool:
            for _ in pool.imap_unordered(run_create_ei, ei_args_list):
                pass

        # Merge all temporary pickles into a single HDF5 file
        print("Consolidating temporary files into HDF5...")
        files = os.listdir(final_tempdir)
        files.sort()

        curr_output_file = f"{args.output}/edges_{low_th:.1f}_{high_th:.1f}.h5"
        store = pd.HDFStore(curr_output_file, mode='w')

        for f in files:
            et = pd.read_pickle(os.path.join(final_tempdir, f))
            et = et[(et['corr'] >= low_th) & (et['corr'] < high_th)]
            if not et.empty:
                store.append('e', et, format='t', data_columns=True, index=False,
                             min_itemsize={'row_name_s': 35, 'row_name_t': 35})
        store.close()
        print(f"Saved {curr_output_file}")

        # Clean temp dir for next run
        reset_dir(final_tempdir)

    print("\nAll threshold intervals processed successfully.")


def entry_point():
    main()


if __name__ == "__main__":
    entry_point()