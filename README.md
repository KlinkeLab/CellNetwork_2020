# **Data-driven learning how oncogenic gene expression locally alters heterocellular networks**

This repository supplies the code developed in the study of D.J. Klinke, A. Fernandez, W. Deng, H. Latifizadeh, and A.C. Pirkey **_"Data-driven learning how oncogenic gene expression locally alters heterocellular networks"_**. The corresponding pre-print can be found on bioRxiv ( doi: https://doi.org/10.1101/2020.05.04.077107). It can be used to reproduce the results of the study and investigate the methodology to be used for other datasets.

# **System Requirements**
## **Hardware Recommendations:**
The code contained in this repository requires only a standard computer to run, however, the amount of RAM required will be defined by the desires of the use. For low level performance, 2 GB of RAM should suffice. For optimal performance, a minimum of 16 GB of RAM and 8 Cores with 2.9 GHz/core is recommended. These specifications match the laboratory computer specifications that the code was tested on.

The runtimes below were generated with a computer with the optimal specifications (16 GB RAM and 8 Cores @ 2.9 GHz/core) and an internet download speed of 791 Mbps.
## **Software Recommendations:**
* Computer Operating System: Tested on Windows 10, version 21H2
* The code was originally written using R version 3.6.3. along with a corresponding version of RStudio. Most recent testing done using R version 4.1.0 and RStudio version 1.4.1717
* To use the provided code, you will need to install the following R libraries: BiocManager (v. 1.30.16), SummarizedExperiment (v. 1.22.0), TCGAbiolinks (v. 2.20.1), colorspace (v. 2.0.2), MASS (v. 7.3.54), bnlearn (v. 4.7), Rgraphviz (v. 2.36.0), parallel (v. 4.1.0), stats (v. 4.1.0), lattice (v. 0.20.44)

# **Data**

All necessary data is provided via links to public data provided in the scripts.

For more information about data used in this study please refer to [https://portal.gdc.cancer.gov] for Files: the breast cancer (BRCA) and Skin Cutaneous Melanoma (SKCM) arms of the Cancer Genome Atlas downloaded via `TCGAbiolinks' (V package in R accessed 6/27/2017; and [https://www.ncbi.nlm.nih.gov] entry GSE98394.

# **Instructions**
## **Installation Guide**
* The latest version of R can be downloaded at https://www.r-project.org/ and should download in under one minute
* The latest version of RStudio can be downloaded at https://www.rstudio.com/products/rstudio/download/ and should install in under one minute
* Download the files from the repository and extract the files from the zip folder to your desired working folder. The extraction process should be near instantaneous.
* Open the six provided R scripts and change the address in the setwd() command on each file to reflect the folder you saved the files to. The scripts should open instanteously.
* Ensure all R Libraries listed in the Software Requirements are installed. Code for the installation of all packages is provided below.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
    BiocManager::install("SummarizedExperiment")
    BiocManager::install("TCGAbiolinks")
    BiocManager::install("Rgraphviz")
    
    install.packages(c('colorspace','MASS', 'lattice'))
    install.packages("https://www.bnlearn.com/releases/bnlearn_latest.tar.gz", repos = NULL, type = "source")
 ```
## **Demo/Instructions for Use**
Running the scripts in the order provided will generate all relevant figures and data tables for the given portion of the study. Below you'll find a list of each script, it's run time, the files added to the working directory, and the number of figures generated.

* 1_BRCA_PreProcess.R
  * Run time:  14.5 minutes
  * Files created: 
    * GDCdata file folder
    * BRCA_RNAseq_HK.txt
    * BRCA_TPM_HK.Rda
    * brca-RNAseq-Counts.rda
  * Figures generated: 1

* 2_BRCA-EMT-SM.R
  * Run time:  < 10 seconds
  * Files created: 
     * BRCA_EMT_SM.csv
  * Figures generated: 2
  
* 3_BRCA-Proliferation-SM.R
  * Run time:  < 10 seconds
  * Files created: 
    * BRCA_proliferation_Zscore.csv
  * Figures generated: 1

* 4_Merge_Datasets.R
  * Run time:  < 10 seconds
  * Files created: 
    * BRCA_digitalcytometry_tot_noise.csv
  * Figures generated: 7

* 5_BRCA_BN_SeedEdges_Interval.R
  * Run time:  23 minutes
  * Files created: 
    * BRCA_HYBRID_CIBERSORT Bootstrap 10000_**[insert algorithm]**.csv 
        * algorithms: gs, rsmax2, mmhc, tabu, hc, pc_stable, iamb_fdr, si_hiton_pc, mmpc, iamb
    * BRCA-Hybrid-CorrelationStructure.pdf
    * BRCA-HybridVariableDist.pdf
    * BRCA_hybrid_CIBERSORT_Seed_Summary_Int.csv 
  * Figures generated: 1

* 6_BRCA_BN_EnsembleDAG.R
  * Run time: 2 minutes
  * Files created: 
    * BRCA_CIBERSORT_NKsubset_L1_mmhc.csv
    * BRCA_CIBERSORT_NKsubset_mmhc_Hybrid_Merge.csv
    * BRCA_ModelvsData-hybrid
  * Figures generated: 1
 
# **General notes**

The code provided in this repository reproduces the main results of the study of D.J. Klinke, A. Fernandez, W. Deng, H. Latifizadeh, and A.C. Pirkey **_"Data-driven learning how oncogenic gene expression locally alters heterocellular networks"_** but it is not meant as a self-contained module for use in further analysis.

# **Citation**

Klinke, D.J., Fernandez, A., Deng, W., Latifizadeh, H. & Pirkey, A.C.,  **_"Data-driven learning how oncogenic gene expression locally alters heterocellular networks."_** bioRxiv (2020). [https://doi.org/10.1101/2020.05.04.077107]
