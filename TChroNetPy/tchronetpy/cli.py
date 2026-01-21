import argparse
import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
import deepgraph as dg
import shutil
from scipy import stats
import pyarrow as pa
import pyarrow.parquet as pq

class ei_args:
    def __init__(self, *, in_graph, in_matrix, in_n_samples, in_th, in_max_th, in_output_dir, in_pos_array, in_step_size):
        self.graph = in_graph
        self.data_matrix = in_matrix
        self.n_samples = in_n_samples
        self.th = in_th
        self.max_th = in_max_th
        self.output_dir = in_output_dir
        self.pos_array = in_pos_array
        self.step_size = in_step_size

def reset_dir(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)

def corr_from_pvalue(p_value, n_samples):
    dof = n_samples - 2
    t_stat = stats.t.ppf(1 - p_value/2, df=dof)
    r = round(np.sqrt(t_stat**2 / (t_stat**2 + dof)), 3)
    return r

def spearman_corr(index_s, index_t, data_matrix, n_samples, min_th):
    features_s = data_matrix[index_s].astype(np.float32)
    features_t = data_matrix[index_t].astype(np.float32)
    corr = np.einsum('ij,ij->i', features_s, features_t) / n_samples
    corr = np.clip(corr, -1.0, 1.0)
    
    # Apply threshold
    significant_mask = corr < min_th
    corr[significant_mask] = 0
    
    # --- ROUNDING TO 3 DECIMAL PLACES ---
    np.round(corr, 3, out=corr)
    return corr

def run_create_ei(args):
    return create_ei(*args)

def create_ei(i, ei_args):
    graph = dg.DeepGraph(ei_args.graph.v.copy(deep=True))
    pos_array = ei_args.pos_array
    chunk_size = 10000 
    
    bins = np.round(np.arange(ei_args.th, ei_args.max_th + 0.1, 0.1), 1)
    v_lookup = ei_args.graph.v.set_index('index')['row_name']
    
    worker_buffer = {round(b, 1): [] for b in bins[:-1]}
    from_pos, to_pos = pos_array[i], pos_array[i + 1]

    def connector_func(index_s, index_t):
        corr = spearman_corr(index_s, index_t, ei_args.data_matrix, ei_args.n_samples, ei_args.th)
        return corr

    for chunk_start in range(from_pos, to_pos, chunk_size):
        graph.create_edges(
            connectors=connector_func,
            step_size=ei_args.step_size,
            from_pos=chunk_start,
            to_pos=min(chunk_start + chunk_size, to_pos)
        )

        if not graph.e.empty:
            if ei_args.max_th < 1.0:
                graph.e = graph.e[graph.e['corr'] <= ei_args.max_th]
            
            if not graph.e.empty:
                # Clean names and force string type
                graph.e['row_name_s'] = graph.e.index.get_level_values('s').map(v_lookup).astype(str).str.strip()
                graph.e['row_name_t'] = graph.e.index.get_level_values('t').map(v_lookup).astype(str).str.strip()

                for low in worker_buffer.keys():
                    high = round(low + 0.1, 1)
                    
                    # If this is the last bin (e.g., 0.9 to 1.0), make the high boundary inclusive
                    if high >= 1.0:
                        mask = (graph.e['corr'] >= low) & (graph.e['corr'] <= 1.0)
                    else:
                        mask = (graph.e['corr'] >= low) & (graph.e['corr'] < high)
                    
                    if mask.any():
                        worker_buffer[low].append(graph.e[mask].reset_index(drop=True))

        graph.e = graph.e.iloc[0:0]

    for low, data_list in worker_buffer.items():
        if data_list:
            final_df = pd.concat(data_list, ignore_index=True)
            table = pa.Table.from_pandas(final_df, preserve_index=False)
            bin_path = os.path.join(ei_args.output_dir, f"bin_{low:.1f}")
            os.makedirs(bin_path, exist_ok=True)
            out_file = os.path.join(bin_path, f"worker_{i}.parquet")
            pq.write_table(table, out_file, compression='zstd')

def main():
    parser = argparse.ArgumentParser(prog="tchronetpy")
    parser.add_argument("-m", "--matrix", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, required=True)
    parser.add_argument("-s", "--stepsize", type=int, default=1e5)
    parser.add_argument("-@", "--threads", type=int, default=1)
    parser.add_argument("--max_t_val", type=float, default=1.0)
    parser.add_argument("--min_t_type", type=str, default="pval", choices=["pval", "cor"])
    parser.add_argument("--min_t_val", type=float, default=0.1)
    args = parser.parse_args()

    reset_dir(args.output)
    df_input = pd.read_csv(args.matrix, sep="\t", index_col=0)
    row_names = list(df_input.index)
    matrix = df_input.to_numpy()
    
    # Pre-rank
    matrix_rank = stats.rankdata(matrix, axis=1)
    matrix_rank = (matrix_rank - matrix_rank.mean(axis=1, keepdims=True)) / matrix_rank.std(axis=1, keepdims=True)

    if args.min_t_type == "pval":
        actual_min_th = corr_from_pvalue(args.min_t_val, matrix.shape[1])
    else:
        actual_min_th = args.min_t_val

    v = pd.DataFrame({'index': range(matrix.shape[0]), 'row_name': row_names})
    pos_array = np.linspace(0, matrix.shape[0]*(matrix.shape[0]-1)//2, args.threads + 1, dtype=int)
    
    ei_args_obj = ei_args(
        in_graph=dg.DeepGraph(v), in_matrix=matrix_rank, in_n_samples=matrix.shape[1],
        in_th=actual_min_th, in_max_th=args.max_t_val, in_output_dir=args.output,
        in_pos_array=pos_array, in_step_size=args.stepsize
    )

    with Pool(processes=args.threads) as pool:
        pool.map(run_create_ei, [(i, ei_args_obj) for i in range(args.threads)])

if __name__ == "__main__":
    main()