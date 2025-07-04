# TRACES
TRACES: **T**emo**R**al **A**ware **C**o-acc**ES**sibility network

## Requirments
Create a conda environment : 
```
conda create -n TRACES os pandas numpy multiprocessing deepgraph sys igraph matplotlib seaborn leidenalg scipy
```
R packages for ontologies analysis :
- rGREAT
- tidyverse

## Execution
To run the TRACES.py script, there are some positional arguments to be set :
- path to the normalized counts matrix
- the name of the output edge list file
- number of elements in each chunk per thread (tested with 100000)
- number of threads to use
- path to a temp dir to be used during the execution

The other two scripts are useful for analyzing the network :
- **TRACES_regnet_Analysis.ipynb** loads the network and uses leiden algorithm to find communities
- **TRACES_rGREAT.R** allows the interrogation of GREAT to define enriched pathways for each community
