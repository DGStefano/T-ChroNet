"""
TCHRONET.py

This script constructs a correlation-based network from a feature matrix using 
Spearman correlations. It iteratively lowers the correlation threshold until 
the resulting network is fully connected.

Main workflow:
1. Load input matrix and preprocess for Spearman correlation
2. Compute pairwise correlations in parallel
3. Filter edges by correlation threshold
4. Check if resulting graph is connected
5. If not connected, lower threshold and repeat
6. Save final connected graph

"""

import argparse
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
import deepgraph as dg
import sys
import shutil
import igraph as ig

class ei_args:
    """
    Container class for edge computation arguments.
    
    Stores all parameters needed for parallel edge creation to avoid
    passing multiple arguments to multiprocessing functions.
    
    Attributes:
        graph: DeepGraph object containing vertex information
        pos_array: Array defining position ranges for parallel processing
        data_matrix: Preprocessed data matrix for correlation computation
        n_samples: Number of samples (columns) in the original matrix
        tmp_dir: Directory for storing temporary pickle files
        th: Correlation threshold for edge filtering
    """
    def __init__(self, *, in_graph, in_pos_array, in_matrix, in_n_samples, in_tmp_dir, in_th):
        self.graph = in_graph
        self.pos_array = in_pos_array
        self.data_matrix = in_matrix
        self.n_samples = in_n_samples
        self.tmp_dir = in_tmp_dir
        self.th = in_th


def reset_dir(path):
    """
    Deletes a directory if it exists, then recreates it.
    
    WARNING: This permanently deletes all contents of the directory!
    
    Args:
        path: Path to the directory to reset
    """
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


def spearman_corr(index_s, index_t, data_matrix, n_samples):
    """
    Compute pairwise Spearman correlations between features.
    
    This function calculates correlations using a preprocessed rank-transformed
    matrix, making the computation efficient through vectorized operations.
    
    Args:
        index_s: Source feature indices
        index_t: Target feature indices
        data_matrix: Preprocessed (rank-transformed and standardized) data matrix
        n_samples: Number of samples for normalization
    
    Returns:
        Array of correlation values rounded to 3 decimal places
    """
    # Extract features for source and target indices
    features_s = data_matrix[index_s].astype(np.float32)
    features_t = data_matrix[index_t].astype(np.float32)
    
    # Compute correlation using Einstein summation (efficient dot product)
    # This calculates the Spearman correlation from the preprocessed ranks
    corr = np.einsum('ij,ij->i', features_s, features_t) / n_samples
    
    # Round to 3 decimal places to save memory
    np.round(corr, 3, out=corr)
    return corr


def run_create_ei(args):
    """
    Wrapper function for create_ei to enable use with multiprocessing.Pool.
    
    Args:
        args: Tuple of (index, ei_args_obj)
    
    Returns:
        Result of create_ei function
    """
    return create_ei(*args)


def create_ei(i, ei_args):
    """
    Create edges for a specific chunk of the correlation matrix in parallel.
    
    This function processes a portion of all possible pairwise correlations,
    filters by threshold, and saves results to temporary files. Processing
    is done in smaller chunks to manage memory usage.
    
    Args:
        i: Index of the processing chunk (for parallel execution)
        ei_args: ei_args object containing all necessary parameters
    """
    # Create a deep copy of the graph structure for this worker
    graph = dg.DeepGraph(ei_args.graph.v.copy(deep=True))
    pos_array = ei_args.pos_array
    step_size = 1e5  # Number of edges to process at once in DeepGraph
    
    chunk_size = 10000  # Process edges in chunks of this size
    th = ei_args.th  # Correlation threshold
    
    # Determine the range of edges this worker should process
    from_pos = pos_array[i]
    to_pos = pos_array[i + 1]
    
    def local_spearman_corr(index_s, index_t):
        """
        Local wrapper for spearman_corr that captures ei_args from closure.
        """
        corr = spearman_corr(index_s, index_t, ei_args.data_matrix, ei_args.n_samples)
        return corr
    
    # Process edges in smaller chunks to manage memory
    edge_chunks = range(from_pos, to_pos, chunk_size)
    for chunk_start in edge_chunks:
        # Create edges using DeepGraph's connector function
        graph.create_edges(
            connectors=local_spearman_corr,
            step_size=step_size,
            from_pos=chunk_start,
            to_pos=min(chunk_start + chunk_size, to_pos)
        )
        
        # Filter edges: keep only those above threshold (hard thresholding)
        graph.e = graph.e[graph.e['corr'] > th]
        
        # Map numeric indices back to original row names for interpretability
        graph.e['row_name_s'] = graph.e.index.get_level_values('s').map(graph.v.set_index('index')['row_name'])
        graph.e['row_name_t'] = graph.e.index.get_level_values('t').map(graph.v.set_index('index')['row_name'])
        
        # Reset index and use row names as new index
        graph.e = graph.e.reset_index()  # Remove the old MultiIndex
        graph.e.drop(['s','t'], axis=1, inplace=True)  # Drop numeric indices
        graph.e = graph.e.set_index(['row_name_s', 'row_name_t'])
        
        # Save this chunk to a temporary pickle file
        tmp_save_pickle = os.path.join(ei_args.tmp_dir, f'{i}_{chunk_start}.pickle')
        graph.e.to_pickle(tmp_save_pickle)
        
        # Clear the DataFrame to free memory for the next chunk
        graph.e.drop(graph.e.index, inplace=True)


def main():
    """
    Main execution function for TCHRONET.
    
    Workflow:
    1. Parse command-line arguments
    2. Load and preprocess input matrix
    3. Iteratively compute correlations with decreasing thresholds
    4. Check graph connectivity after each threshold
    5. Save final connected graph
    """
    
    #################################
    ##### Parsing args
    
    parser = argparse.ArgumentParser(
                        prog="TCHRONET")
    
    parser.add_argument("-m", "--matrix", type=str, required=True, 
                        help="Path to the input matrix (TSV format, features as rows)")
    parser.add_argument("-o", "--output", type=str, required=True, 
                        help="Path to the output file (HDF5 format)")
    parser.add_argument("-s", "--stepsize", type=int, required=False, default=1e5, 
                        help="Stepsize for RAM parameters (DeepGraph parameter)")
    parser.add_argument("-@", "--threads", type=int, required=False, default=8, 
                        help="Number of threads to use for parallel processing")
    parser.add_argument("-t", "--tempdir", type=str, required=False, 
                        help="Temporary directory to use. Will create a '/tmp' directory inside it")
    parser.add_argument("-r", "--threshold", type=float, required=False, default=0.8, 
                        help="Starting correlation threshold (will decrease if graph not connected)")
    
    args = parser.parse_args()
    
    
    #################################
    ##### Arg managing
    
    tempfolder_name = "tchronet_tmp"
    
    # Determine temporary directory location
    if args.tempdir is not None:
        # User specified a custom temp directory
        abs_args_tempdir = os.path.abspath(args.tempdir)
        final_tempdir = os.path.join(abs_args_tempdir, tempfolder_name)
        final_tempdir = os.path.abspath(final_tempdir)
        reset_dir(final_tempdir)
    else:
        # Use current directory for temp files
        final_tempdir = os.path.abspath('.')
        final_tempdir = os.path.join(final_tempdir, tempfolder_name)
        final_tempdir = os.path.abspath(final_tempdir)
        reset_dir(final_tempdir)
    
    print(f"Using: {final_tempdir} as temporary folder")
    
    
    #################################
    ##### Main method
    
    # Load the input matrix (features x samples)
    print("Loading input matrix...")
    data_matrix_df = pd.read_csv(args.matrix, sep="\t", index_col=0)
    row_names = list(data_matrix_df.index)
    data_matrix = data_matrix_df.to_numpy()
    n_features = data_matrix.shape[0]
    n_samples = data_matrix.shape[1]
    
    # Preprocess matrix for efficient Spearman correlation computation
    # Step 1: Convert values to ranks (argsort twice gives ranks)
    print("Preprocessing matrix for Spearman correlation...")
    data_matrix = data_matrix.argsort(axis=1).argsort(axis=1)
    
    # Step 2: Standardize ranks (mean=0, std=1) for correlation computation
    data_matrix = (data_matrix - data_matrix.mean(axis=1, keepdims=True)) / data_matrix.std(axis=1, keepdims=True)
    
    # Create vertex DataFrame with feature indices and names
    v = pd.DataFrame({'index': range(data_matrix.shape[0]), 'row_name': row_names})
    
    # Set up position array for dividing work among parallel processes
    if args.threads == 1:
        # Single-threaded: process all edges at once
        pos_array = np.array([0, n_features * (n_features - 1) // 2], dtype=int)
        indices = np.array([0])
    else:
        # Multi-threaded: divide edge space evenly among threads
        pos_array = np.array(np.linspace(0, n_features * (n_features - 1) // 2, args.threads), dtype=int)
    indices = np.arange(0, args.threads - 1)
    
    # Create parameter object to pass to worker processes
    ei_args_obj = ei_args(
        in_graph=dg.DeepGraph(v),
        in_pos_array=pos_array,
        in_matrix=data_matrix,
        in_n_samples=n_samples,
        in_tmp_dir=final_tempdir,
        in_th=args.threshold
    )
    
    # Track files created at different thresholds
    iteration_th_filelist = []
    last_connected = False
    last_th = ei_args_obj.th
    
    # Iterative threshold reduction loop
    # Continue until we find a threshold that produces a connected graph
    while not last_connected:
        print(f"Starting multiprocess with threshold: {ei_args_obj.th}")
        
        # Prepare arguments for each parallel worker
        ei_args_list = [(i, ei_args_obj) for i in indices]
        
        # Execute parallel edge computation
        with Pool(processes=args.threads) as pool:
            for _ in pool.imap_unordered(run_create_ei, ei_args_list):
                pass
        
        # Create output filename for this threshold iteration
        curr_output_file = args.output
        curr_output_file += f".th.{ei_args_obj.th}"
        iteration_th_filelist.append(curr_output_file)
        
        # Consolidate all temporary pickle files into a single HDF5 file
        print("Consolidating temporary files into HDF5...")
        files = os.listdir(f'{final_tempdir}/')
        files.sort()
        
        store = pd.HDFStore(curr_output_file, mode='w')
        for f in files:
            # Read each temporary pickle file
            et = pd.read_pickle(os.path.join(final_tempdir, f))
            # Append to HDF5 store with specified column sizes
            store.append('e', et, format='t', data_columns=True, index=False,
                        min_itemsize={'row_name_s': 35, 'row_name_t': 35})
        store.close()
        
        # Clean up temporary directory for next iteration
        reset_dir(final_tempdir)
        
        
        # Check if the resulting graph is connected
        print("Loading thresholded file...")
        e = pd.read_hdf(curr_output_file)
        e_df = e.reset_index()
        
        # Free memory
        del e
        
        # Build graph using igraph for connectivity testing
        print("Building graph and checking connectivity...")
        G = ig.Graph.TupleList(e_df.itertuples(index=False), directed=False, 
                               weights=None, edge_attrs='corr')
        
        last_th = ei_args_obj.th
        
        if G.is_connected():
            # Success! Graph is connected at this threshold
            print(f"Graph connected with threshold: {ei_args_obj.th}")
            print("Stopping here!")
            last_connected = True
            
            # Remove this file from the cleanup list and rename to final output
            iteration_th_filelist.remove(curr_output_file)
            os.rename(curr_output_file, args.output)
            
            # Clean up all intermediate threshold files
            for thfile in iteration_th_filelist:
                os.remove(thfile)
        else:
            # Graph not connected, need to lower threshold
            print(f"Graph not connected at threshold: {ei_args_obj.th}")
            print("Going down by 0.1")
            last_connected = False
            ei_args_obj.th = last_th - 0.1
            
            # Safety check: don't go below 0.1 threshold
            if ei_args_obj.th < 0.1:
                print("Lowest threshold reached, still not connected. Exiting.")
                exit(-1)


if __name__ == '__main__':
    main()
