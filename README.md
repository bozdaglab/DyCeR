# DyCeR
DyCeR is a tool to infer dysregulated ceRNA crosstalks in cancer.  

Scripts to run DyCeR are in **run_DyCeR.R**. Required inputs for DyCeR include: 
* Four dataframes contains gene-level information of mRNA and miRNA expression, copy number alteration, and DNA methylation. The corresponding variable names used in **run_DyCeR.R** are _mRNA_, _miRNA_, _cna_, and _methyl_. 
  - Each column in those dataframes refer to a sample. The column names and order in all four dataframes must be the same. 
  - Each row in the _mRNA_, _cna_, and _methyl_ refers to mRNA expression level, copy number, and methylation beta value of a differenitally expressed (DE) gene, respectively. The row names and order in _mRNA_, _cna_, and _methyl_ must be the same. Each row in _miRNA_ refers to the expression of a DE miRNA.  
* A dataframe named _tf.target.interactions_. Each row in _tf.target.interactions_ refers to a putative interaction between a transcription factor (TF) and a gene. 
* A dataframe named _miRNA.target.interactions_. Each row in _miRNA.target.interactions_ refers to a putative interaction between a miRNA and a gene. 
  
Sample input data for DyCeR could be found in _sampleData.rda_. 

