# T-ChroNet
T-ChroNet (**T**ime-aware **Chro**matin **Net**work)

## Abstract
Networks are widely applied to investigate relationships among the individual components of complex biological systems. Recent application of biological networks, such as gene co-expression networks and gene regulatory networks, has illuminated the principles underlying transcriptional modulation in development and diseases. However, computational methods that can embed the activity of cis-regulatory elements (CRE) into a network are still limited. Capturing temporal CRE activity within a network could help revealing regulatory programs involved in cell fate commitment and disease development. To address this, we present T-ChroNet (Time-aware Chromatin Network), a network-based method that models CRE as nodes and their temporal co-accessibility as edges. Through the detection of CRE sharing similar accessibility patterns over time, T-ChroNet allows the inference of upstream regulons and downstream biological pathways. We applied T-ChroNet to temporally-resolved CRE datasets, from both human and mouse, including both chromatin accessibility (ATACseq) and histone post-translational modifications (H3K27ac ChIPseq). T-ChroNet successfully recovered known regulators and enriched pathways for both modalities and species, while also uncovering novel putative factors and mechanisms regulating cell identity, organ development and disease progression.
<img width="2481" height="2234" alt="GraphicalAbstract" src="https://github.com/user-attachments/assets/e03c44e3-df8c-485f-a93e-d0e313a05a31" />

## Requirments
Create a conda environment : 
```
conda create -n TCHRONET os pandas numpy multiprocessing deepgraph sys igraph matplotlib seaborn leidenalg scipy
```
R packages for ontology analysis :
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
