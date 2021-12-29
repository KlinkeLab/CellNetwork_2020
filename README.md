# **Data-driven learning how oncogenic gene expression locally alters heterocellular networks**

This repository supplies the code developed in the study of D.J. Klinke, A. Fernandez, W. Deng, and A.C. Pirkey **_"Data-driven learning how oncogenic gene expression locally alters heterocellular networks"_**. The corresponding pre-print can be found on bioRxiv ( doi: https://doi.org/10.1101/2020.05.04.077107). It can be used to reproduce the results of the study and investigate the methodology to be used for other datasets.

## **Requirements**

* Computer Operating System: Windows 10 or MAC OS ...
* R version 3.6.3.
* R libraries: SummarizedExperiment, TCGAbiolinks, colorspace, MASS, bnlearn, BiocManager, Rgraphviz, parallel, stats, lattice

## **Data**

All necessary data is provided via links to public data provided in the scripts.

For more information about data used in this study please refer to [https://portal.gdc.cancer.gov] for Files: the breast cancer (BRCA) and Skin Cutaneous Melanoma (SKCM) arms of the Cancer Genome Atlas downloaded via `TCGAbiolinks' (V package in R accessed 6/27/2017; and [https://www.ncbi.nlm.nih.gov] entry GSE98394.

## **Quick start**

To reproduce the results, download the relevant script and data and load the corresponding data into the R workspace. Running the scripts will generate all relevant figures and data tables for the given portion of the study.

# General notes

The code provided in this repository reproduces the main results of the study of D.J. Klinke, A. Fernandez, W. Deng and A.C. Pirkey **_"Data-driven learning how oncogenic gene expression locally alters heterocellular networks"_** but it is not meant as a self-contained module for use in further analysis.

## Citation

Klinke, D.J., Fernandez, A., Deng, W., & Pirkey, A.C.,  **_"Data-driven learning how oncogenic gene expression locally alters heterocellular networks."_** bioAxiv (2020). [https://doi.org/10.1101/2020.05.04.077107]
