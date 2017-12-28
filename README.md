
# DyCeR
DyCeR is a tool to infer dysregulated ceRNA crosstalk in cancer.  

Scripts to run DyCeR are in **run_DyCeR.R**. Required input for DyCeR include: 
* Four dataframes for gene-level information of mRNA and miRNA expression, copy number alteration, and DNA methylation. The corresponding variable names used in **run_DyCeR.R** are _mRNA_, _miRNA_, _cna_, and _methyl_. 
  - Each column in those dataframes refer to a sample. The column name and order in all four dataframes must be the same. Each row in the mRNA expression, copy number alternation
  - 