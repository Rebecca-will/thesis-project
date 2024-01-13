# A Consensus Gene Analysis of The Cancer Genome Atlas
##### In partnership with ALMAC

This project focuses on analyzing transcriptomic data from The Cancer Genome Atlas (TCGA) with the overarching goal of identifying a pan-cancer molecular signal associated with the cancer hallmark of invasion and metastasis across various solid cancer types. The approach employed was unsupervised and non-hypothesis-driven, allowing for the unbiased discovery of patterns within the data.

Initial analysis began with the generation of a heatmap with hierarchical clustering. This visualisation technique facilitated the identification of molecular subgroups based on transcriptomic data. Subsequently, enrichment and survival analyses were conducted on these molecular subgroups to gain insights into specific characteristics associated with each group.

To indetify the important genes to distinguish the groups, a gene signature was developed using XGBOOST, a  gradient-boosted decision tree machine learning algorithm. This algorithm was utilised due to its ability to handle nonlinear relationships, feature importance analysis, robustness to overfitting, capacity to handle missing data, ensemble learning, scalability, and regularisation.

###### Heatmap 3
This script contains the code used to:
- Scale the genetic data.
- Identifying the genes with the strongest signal intensity to reduce noise. 
- Perform hierarchical clustering using wards linkage.
- A legend was added to identify the different cancer types.
- Enrichment analysis of the subgroups to understand the breakdown of cancer type within each group.
- Creation of pie charts to illustrate the subgroup's contents.

###### Survival_Final/ 2 cluster analysis
These scripts contain the code used to: 
- Add the corresponding subgroup identifier to the samples in the survival dataset.
- Investigation into why the datasets have differing numbers of samples (not every sample had follow-up data so it was not included).
- General data cleaning and preparation for analysis.
- Survival and disease-free progression analysis using Cox regression analysis.
- Plotting the Kaplan-Meier survival and disease-free progression curves.
- Breakdown of survival and characteristics of individual diseases.
- Analysis of TNM characteristics of the subgroups.

###### XGBOOST/gene_sig_dev/final model
These scripts contain the code used to:
- Different iterations of testing and training the xgboost model to develop a gene signature as the project progressed.
- Prediction accuracy and feature importance analysis of the models.
- Hyperparameter tuning and cross-fold validation of the models.
- Analysis of the characteristics of the classified samples.
