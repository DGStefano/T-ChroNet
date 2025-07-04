import argparse
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
import deepgraph as dg
import sys
import shutil

class ei_args:
    def __init__(self, *, in_graph, in_pos_array, in_matrix, in_n_samples, in_tmp_dir):
        self.graph = in_graph
        self.pos_array = in_pos_array
        self.data_matrix = in_matrix
        self.n_samples = in_n_samples
        self.tmp_dir = in_tmp_dir

# Deletes a directory, if it exists, and recreates it. Watch out!
def reset_dir(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


# connector function to compute pairwise Speraman correlations
def spearman_corr(index_s, index_t, data_matrix, n_samples):
    features_s = data_matrix[index_s].astype(np.float32)
    features_t = data_matrix[index_t].astype(np.float32)
    corr = np.einsum('ij,ij->i', features_s, features_t) / n_samples
    np.round(corr, 3, out=corr)
    return corr


def run_create_ei(args):
    return create_ei(*args)


# parallel computation
def create_ei(i, ei_args):
    
    print(f"Starting create_ei with index {i}")
    graph = dg.DeepGraph(ei_args.graph.v.copy(deep=True))
    pos_array = ei_args.pos_array
    step_size = 1e5
    # n_processes = 40

    chunk_size = 10000
    th = 0.8

    from_pos = pos_array[i]
    to_pos = pos_array[i + 1]

    def local_spearman_corr(index_s, index_t):
        corr = spearman_corr(index_s, index_t, ei_args.data_matrix, ei_args.n_samples)
        return corr

    # Process in smaller chunks
    edge_chunks = range(from_pos, to_pos, chunk_size)
    for chunk_start in edge_chunks:
        graph.create_edges(
            connectors=local_spearman_corr,
            step_size=step_size,
            from_pos=chunk_start,
            to_pos=min(chunk_start + chunk_size, to_pos)
        )
        
        # hard thresholding
        graph.e = graph.e[graph.e['corr'] > th]
        
        # Add row names to the edges DataFrame
        graph.e['row_name_s'] = graph.e.index.get_level_values('s').map(graph.v.set_index('index')['row_name'])
        graph.e['row_name_t'] = graph.e.index.get_level_values('t').map(graph.v.set_index('index')['row_name'])
        
        graph.e = graph.e.reset_index()  # Remove the old MultiIndex
        graph.e.drop(['s','t'] , axis = 1 , inplace = True)
        graph.e = graph.e.set_index(['row_name_s', 'row_name_t'])

        tmp_save_pickle = os.path.join(ei_args.tmp_dir, f'{i}_{chunk_start}.pickle')

        # Save the result
        graph.e.to_pickle(tmp_save_pickle)
        graph.e.drop(graph.e.index, inplace=True)  # Clear the DataFrame to free memory


def main():

    #################################
    ##### Parsing args

    parser = argparse.ArgumentParser(
                        prog="TRACES",
                        description="TempoRally Aware Co-accESsibility network")

    parser.add_argument("-m", "--matrix", type=str, required=True, help="Path to the input matrix")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output file")
    parser.add_argument("-s", "--stepsize", type=int, required=False, default=1e5, help="Stepsize for RAM parameters")
    parser.add_argument("-@", "--threads", type=int, required=False, default=8, help="Number of threads to use")
    parser.add_argument("-t", "--tempdir", type=str, required=False, help="Temporary directory to use. Will create a \"/tmp\" directory inside it.")

    args = parser.parse_args()


    #################################
    ##### Arg managing

    tempfolder_name = "traces_tmp"
    # If the user hasn't requested a particular tempdir, we create one in this directory
    if args.tempdir is not None:
        # Delete and recreate tmp folder
        abs_args_tempdir = os.path.abspath(args.tempdir)
        final_tempdir = os.path.join(abs_args_tempdir, tempfolder_name)
        final_tempdir = os.path.abspath(final_tempdir)
        reset_dir(final_tempdir)
    else:
        # Same here
        final_tempdir = os.path.abspath('.')
        final_tempdir = os.path.join(final_tempdir, tempfolder_name)
        final_tempdir = os.path.abspath(final_tempdir)
        reset_dir(final_tempdir)

    print(f"Using: {final_tempdir} as temporary folder")

    #################################
    ##### Main method

    # Loading counts matrix
    data_matrix_df = pd.read_csv( args.matrix , sep = "\t" , index_col = 0 )
    row_names = list(data_matrix_df.index)
    data_matrix = data_matrix_df.to_numpy()
    n_features = data_matrix.shape[0]
    n_samples = data_matrix.shape[1]
    
    # Organize matrix for further spearman correlation
    data_matrix = data_matrix.argsort(axis=1).argsort(axis=1)
    data_matrix = (data_matrix - data_matrix.mean(axis=1, keepdims=True)) / data_matrix.std(axis=1, keepdims=True)

    # Store the nodes names 
    v = pd.DataFrame({'index': range(data_matrix.shape[0]), 'row_name': row_names})

    # Multiprocessing paramenter
    if args.threads == 1:
        pos_array = np.array([0, n_features * (n_features - 1) // 2], dtype=int)
        indices = np.array([0])
    else:
        pos_array = np.array(np.linspace(0, n_features * (n_features - 1) // 2, args.threads), dtype=int)
    indices = np.arange(0, args.threads - 1)

    # Make the parameter object
    ei_args_obj = ei_args(
        in_graph=dg.DeepGraph(v),
        in_pos_array=pos_array,
        in_matrix=data_matrix,
        in_n_samples=n_samples,
        in_tmp_dir=final_tempdir
    )

    print("Starting multiprocess..")
    # Prepare iterable of argument tuples
    ei_args_list = [(i, ei_args_obj) for i in indices]
    with Pool(processes=args.threads) as pool:
        for _ in pool.imap_unordered(run_create_ei, ei_args_list):
            pass
    
    # Read files in tmp folder
    files = os.listdir(f'{final_tempdir}/')
    files.sort()
    # Read files in tmp folder
    store = pd.HDFStore( args.output, mode='w')
    for f in files:
        et = pd.read_pickle(os.path.join(final_tempdir, f))
        store.append('e', et, format='t', data_columns=True, index=False ,  min_itemsize={'row_name_s': 35 , 'row_name_t': 35})
    store.close()
    shutil.rmtree(final_tempdir)


if __name__ == '__main__':
    main()
