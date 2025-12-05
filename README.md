# T-ChroNet
T-ChroNet (**T**ime-aware **Chro**matin **Net**work)

## Abstract
Networks are widely applied to investigate relationships among the individual components of complex biological systems. Recent application of biological networks, such as gene co-expression networks and gene regulatory networks, has illuminated the principles underlying transcriptional modulation in development and diseases. However, computational methods that can embed the activity of cis-regulatory elements (CRE) into a network are still limited. Capturing temporal CRE activity within a network could help revealing regulatory programs involved in cell fate commitment and disease development. To address this, we present T-ChroNet (Time-aware Chromatin Network), a network-based method that models CRE as nodes and their temporal co-accessibility as edges. Through the detection of CRE sharing similar accessibility patterns over time, T-ChroNet allows the inference of upstream regulons and downstream biological pathways. We applied T-ChroNet to temporally-resolved CRE datasets, from both human and mouse, including both chromatin accessibility (ATACseq) and histone post-translational modifications (H3K27ac ChIPseq). T-ChroNet successfully recovered known regulators and enriched pathways for both modalities and species, while also uncovering novel putative factors and mechanisms regulating cell identity, organ development and disease progression.
<img width="2481" height="2234" alt="GraphicalAbstract" src="https://github.com/user-attachments/assets/e03c44e3-df8c-485f-a93e-d0e313a05a31" />

## Install
To install TChronetPy:
```
cd ~
git clone https://github.com/DGStefano/T-ChroNet.git
cd ~/T-ChroNet
conda env create -f environment.yml
conda activate T-CHRONET
pip install git+https://github.com/DGStefano/T-ChroNet.git#subdirectory=TChroNetPy
```
To install TChronetR
```
remotes::install_github("DGStefano/T-ChroNet",subdir = "TChroNetR")
```

## Tutorial

This short tutorial explains how to run **T-ChroNet** using the included **toy dataset** composed of only chromosome 1 from the THP-1 dataset.

1. Run tchronet tool to evalute correlation parameters among all the peaks: 
```
conda activate T-CHRONET
tchronet -m ~/T-ChroNet/toy_data/thp1_test.tsv -o ~/T-ChroNet/toy_data/results/th/ -@ 1 -t /tmp/thp1_test/ --step 0.1 --pval 0.1
```
2. Use TChronetR for netowrk analysis
The full TChronetR vignette is avaliable [here](./TChroNetR/vignettes/TchronetR_Vignette.Rmd)



## Execution
To run the tchronet tool, there are some positional arguments to be set :
- -m or --matrix : Path to input matrix
- -o or --output : Path to the output file
- -s or --stepsize : Stepsize for RAM parameters
- -@ or --threads : Number of threads to use
- -t or --tempdir : Temporary directory to use. Will create a \"/tmp\" directory inside it
- --min_th : Minimum correlation threshold (start of first interval)
- --max_th : Maximum correlation threshold (end of last interval)
- --step : Size of each correlation interval (e.g., 0.1 for 0.4–0.5, 0.5–0.6, etc.)
- --pval : Pvalue threshold

[![DOI](https://zenodo.org/badge/933292065.svg)](https://doi.org/10.5281/zenodo.16737392)