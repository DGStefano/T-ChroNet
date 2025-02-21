import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
import deepgraph as dg
import sys

# parallel computation
def create_ei(i ):
    
    from_pos = pos_array[i]
    to_pos = pos_array[i + 1]
    g = dg.DeepGraph(v)

    # Process in smaller chunks
    edge_chunks = range(from_pos, to_pos, chunk_size)
    for chunk_start in edge_chunks:
        g.create_edges(
            connectors=corr,
            step_size=step_size,
            from_pos=chunk_start,
            to_pos=min(chunk_start + chunk_size, to_pos)
        )
        
        # hard thresholding
        g.e = g.e[g.e['corr'] > th]
        
        # Add row names to the edges DataFrame
        g.e['row_name_s'] = g.e.index.get_level_values('s').map(v.set_index('index')['row_name'])
        g.e['row_name_t'] = g.e.index.get_level_values('t').map(v.set_index('index')['row_name'])
        
        g.e = g.e.reset_index()  # Remove the old MultiIndex
        g.e.drop(['s','t'] , axis = 1 , inplace = True)
        g.e = g.e.set_index(['row_name_s', 'row_name_t'])

        tmp_save_pickle =  f'{tmp_dir}/{i}_{chunk_start}.pickle'

        # Save the result
        g.e.to_pickle(tmp_save_pickle)
        g.e.drop(g.e.index, inplace=True)  # Clear the DataFrame to free memory

# connector function to compute pairwise Speraman correlations
def corr( index_s, index_t):
    features_s = data_matrix[index_s].astype(np.float32)  # Convert to float32
    features_t = data_matrix[index_t].astype(np.float32)  # Convert to float32
    corr = np.einsum('ij,ij->i', features_s, features_t) / n_samples
    np.round(corr, 3, out=corr)  # In-place rounding
    # return correlation
    return corr

# collecting paramenters
matrix_path = sys.argv[1]
out_filename = sys.argv[2]
step_size = int(sys.argv[3])
n_processes =  int(sys.argv[4])
tmp_dir = sys.argv[5]

# computation
if __name__ == '__main__':

    chunk_size = 10000
    th = 0.8

    # Loading counts matrix
    data_matrix_df = pd.read_csv( matrix_path , sep = "\t" , index_col = 0 )
    row_names = list(data_matrix_df.index)
    data_matrix = data_matrix_df.to_numpy()
    
    # RAM usage parameters
    # step_size = 1e5
    # n_processes = 40
    n_features = data_matrix.shape[0]
    n_samples = data_matrix.shape[1]
    

    # Organize matrix for further spearman correlation
    data_matrix = data_matrix.argsort(axis=1).argsort(axis=1)
    data_matrix = (data_matrix - data_matrix.mean(axis=1, keepdims=True)) / data_matrix.std(axis=1, keepdims=True)

    # Save binary numpy Matrix
#   np.save('/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/ClusteredPatients_Ball_Multicov', data_matrix)    

    # load samples as memory-map
#    X = np.load('/mnt/nas-safu/analysis/PhDsdigiove/method_coAcces/data/BALL/ClusteredPatients_Ball_Multicov.npy', mmap_mode='r')

    # Store the nodes names 
    v = pd.DataFrame({'index': range(data_matrix.shape[0]), 'row_name': row_names})

    # Create tmp folder
    os.makedirs(tmp_dir, exist_ok=True)
    
    # Multiprocessing paramenter
    pos_array = np.array(np.linspace(0, n_features * (n_features - 1) // 2, n_processes), dtype=int)
    indices = np.arange(0, n_processes - 1)
    p = Pool()

    # Calcualte correlation and start multi precesses
    for _ in p.imap_unordered(create_ei, indices):
        pass
    # store correlation values
    
    # Read files in tmp folder
    files = os.listdir(f'{tmp_dir}/')
    files.sort()
    # Read files in tmp folder
    store = pd.HDFStore( out_filename, mode='w')
    for f in files:
        et = pd.read_pickle(f'{tmp_dir}/{f}')
        store.append('e', et, format='t', data_columns=True, index=False ,  min_itemsize={'row_name_s': 35 , 'row_name_t': 35})
    store.close()
