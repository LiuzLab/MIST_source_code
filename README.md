# This repository contains necessary scirpts to reproduce all the results and figures in the MIST manuscript.

## Dependencies
We recommend using a conda environment to automatically install all required dependencies. Conda installation guide can be found at https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html. After installing conda, run the following command to install a ReST environment:

  `conda env create -f environment.yml`

Or, users could manually install all required dependencies as below:

  * pandas=0.25.3
  * numpy=1.18.5
  * matplotlib=3.3.4
  * statsmodels=0.12.0
  * scipy=1.6.1
  * tqdm=4.56.0
  * imageio
  * alphashape
  * descartes
  * joblib
  * gseapy
  * python=3.9

## Reproduce figures in the manuscript
All the code to reproduce the figures and quantification results are in folder `analysis_scripts`. Title of each jupyter notebook indicates which figures it generates.

## Data for manuscripts
Data sets are too large to be shared. All data sets are publicly available with links specified in the supplementary table 1 in the manuscript.

A sample  melanoma data is shared for a demonstration at MIST github repo: https://github.com/linhuawang/MIST.

### Raw count data (without imputation): 
Raw count data for 4 relatively small datasets are provided at our previous repo: https://github.com/LiuzLab/MIST_manuscript.

* ./data/MouseWT/raw.csv
* ./data/MouseAD/raw.csv
* ./data/Melanoma/raw.csv
* ./data/Prostate/raw.csv

The other nine data sets used in the manuscript could be downloaded from [10X](https://www.10xgenomics.com/resources/datasets?query=&page=1&configure%5Bfacets%5D%5B0%5D=chemistryVersionAndThroughput&configure%5Bfacets%5D%5B1%5D=pipeline.version&configure%5BhitsPerPage%5D=500&menu%5Bproducts.name%5D=Spatial%20Gene%20Expression).

### External validation data
Another single-cell RNA-sequencing cohort, which includes 16 mice brains with expressiondata for 37089 single cells published by Methodios Ximerakis in 2019 (GSE129788)was used to further validate the co-expression patterns of the two pairs of genes, Cldn11-Arhgef10and Gfap-Aqp4.

Ximerakis, M. et al.Single-cell transcriptomic profiling of the aging mouse brain. Nat. Neurosci.22, 1696–1708 (2019).

## Holdout experiments
### Run experiments
Resulted data will be available by running `scripts/run_holdout_experiments.sh`. This bash scripts will impute all the 13 datasets using all benchmarked methods.

### Evaluation
Run `evalHO.sh` to generate evaluation results for further analysis and visualization.

## Clustering benchmarking
Notebook `analysis_scripts/2-Figure2f.ipynb` was used to generate benchmarking results from different clustering methods to compare region detection accuracy. This script was also used to evaluate the performance. R script `analysis_scripts/runBayesSpace.R` was used to BayesSpace with different parameters (2 modes). 

## Figures for manuscript

1. Figure 2 and extended (run in the listed order):
* 1-Fig2a-e.ipynb
* 2-Figure2f.ipynb

2. Figure 3 and extended:
* 3-Figure3a-c-holdout-performance.ipynb
* 4-Figure3d-h.ipynb
* 5-Figure3i-and-supplementary.ipynb

3. Figure 4 and extended:
* 6-Figure4a-c&supp.ipynb
* 7-Figure4d-f.ipynb

4. Figure 5 and extended:
* 8-Figure5a&h-ABA.ipynb
* 9-Figure5b-g&i-n.ipynb

5. Other extended and supplementary figures:
* 11-Ext-Fig 7-augmentation-improvement.ipynb

"""
Note to authors: the path to the data on server is at: /houston_20t/alexw/ST/MIST_additional_holdout.
"""
