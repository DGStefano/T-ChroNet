import subprocess
import time
import pandas as pd
import os

# --- CONFIGURATION ---
NODE_INTERVALS = [10000, 20000, 40000, 60000, 80000, 100000]
THREAD_LIST = [1,2,5, 15 ,20,25,30,35,40 ]  # Multi-threading will apply to tchronetpy and WGCNA # 5, 15, 20,25,30,35,40
RESULTS_FILE = "/home/sdigiove/T-ChroNet/paper_analysis/data/banchmark/final_performance_comparison.csv"
OUT_DIR_TCN = "/mnt/nas-safu02/sdigiove_workspace/check_th_TCHRONET/new_tchroent_ties/benchmark/parquet"
def get_command(base_path , tool_name, outdir, threads):
    norm_path = os.path.join(base_path, "norm_counts.tsv")
    raw_path = os.path.join(base_path, "raw_counts.tsv")
    bed_path = os.path.join(base_path, "subset.bed")

    if tool_name == "tchronetpy":
        # Python: Multi-threaded, Normalized data
        return ["tchronetpy",
         "-m" ,str(norm_path),
         "-o" ,str(outdir),
         "-@" , str(threads),
         "--min_t_type" , "pval" ,
         "--min_t_val", str(0.1)]
    
    if tool_name == "TchronetR":
        # TchronetR: Single-threaded, Normalized data
        return ["Rscript", "/home/sdigiove/T-ChroNet/paper_analysis/scripts/benchmark/tchronetr_bench.R",norm_path, OUT_DIR_TCN ]
    
    if tool_name == "WGCNA":
        # WGCNA: Multi-threaded, Raw data
        return ["Rscript", "/home/sdigiove/T-ChroNet/paper_analysis/scripts/benchmark/wgcna_bench.R", str(threads), raw_path]
    
    if tool_name == "TCseq_KM":
        # TCseq: Single-threaded, BED/BAM data logic
        return ["Rscript", "/home/sdigiove/T-ChroNet/paper_analysis/scripts/benchmark/tcseq_km.R", bed_path]

    if tool_name == "TCseq_HC":
        # TCseq: Single-threaded, BED/BAM data logic
        return ["Rscript", "/home/sdigiove/T-ChroNet/paper_analysis/scripts/benchmark/tcseq_hc.R",bed_path]

    return None

def run_benchmark():
    all_results = []
    tools = [ "TchronetR"  ] #"tchronetpy", "TCseq_KM" , "TCseq_HC" , "TCseq_KM" , "TchronetR", "WGCNA" 

    for n in NODE_INTERVALS:
        base_path = "/home/sdigiove/T-ChroNet/paper_analysis/data/banchmark/counts/data_" + str(n)+"/"
        for t in THREAD_LIST:
                for tool in tools:
                    
                    # --- CONSTRAINTS FOR SINGLE-THREADED TOOLS ---
                    # If tool is TCseq or TchronetR, only run when thread count is 1
                    if tool in ["TCseq_KM" , "TCseq_HC", "TchronetR"] and t > 1:
                        continue
                    
                    if tool == "tchronetpy" :
                        conda_wrapper = ["conda", "run", "-n", "tchronet_env" , "--no-capture-output"]
                    else :
                        conda_wrapper = ["conda", "run", "-n", "R" , "--no-capture-output"]

                    if not os.path.exists(base_path):
                        print(f"Skipping interval {n}: Data directory not found.")
                        continue

                    cmd_args = get_command(base_path, tool,OUT_DIR_TCN, t)
                    
                    if not cmd_args : 
                        continue

                    stats_file = "stats.tmp"
                    full_cmd = ["/usr/bin/time", "-v", "-o", stats_file] + conda_wrapper + cmd_args
                    
                    print(f"RUNNING: {tool} | Nodes: {n} | Threads: {t}")

                    start_wall = time.time()
                    try:
                        subprocess.run(full_cmd, check=True)
                        end_wall = time.time()
                        
                        with open(stats_file, "r") as f:
                            stats_output = f.read()
                        
                        mem_line = [line for line in stats_output.split('\n') if "Maximum resident set size" in line][0]
                        peak_kb = int(mem_line.split(":")[-1].strip())
                        peak_mb = peak_kb / (1024) #add * 1024 for Gb
                        
                        all_results.append({
                            "Tool": tool,
                            "Nodes": n,
                            "Threads": t if tool not in ["TCseq", "TchronetR"] else 1,
                            "Runtime_Sec": round(end_wall - start_wall, 4),
                            "Peak_RAM_MB": round(peak_mb, 4),
                            "Peak_RAM_GB": round(peak_mb/ (1024), 4)
                        })
                        
                        # Save progress incrementally
                        pd.DataFrame(all_results).to_csv(RESULTS_FILE, index=False)

                    except Exception as e:
                        print(f"!!! Error running {tool}: {e}")

    print(f"Benchmark complete. Results saved to {RESULTS_FILE}")

if __name__ == "__main__":
    run_benchmark()