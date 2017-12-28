# Set directory --------------------------------------------------------------------------
# dir = "..."
setwd(dir)

setwd("~/Desktop/DyCeR")
# Load libraries  ------------------------------------------------------------------------
require(data.table)
require(parallel)
require(bnlearn)
require(WGCNA)
require(glmnet)

# Load utility functions to run DyCeR ----------------------------------------------------
source("helper_functions.R")

# Load required input for DyCeR  ---------------------------------------------------------
load("sampleData.rda")

# run bootstrapped LASSO procedure for each DE mRNAs -------------------------------------
regression.data = list(mRNA=mRNA, 
                       miRNA=miRNA, 
                       methyl=methyl, 
                       cna=cna,
                       miRNA.target.interactions=miRNA.target.interactions,
                       tf.target.interactions=tf.target.interactions)


# perform the bootstrapped LASSO procedure
lasso.list = get_lasso_result(regression.data = regression.data)

# keep predictors selected more than a threshold
chosen.threshold = 70
lasso.frequent.predictor.list = get_frequently_selected_predictors(lasso.list = lasso.list, 
                                                                   threshold = chosen.threshold)
# get list of predictors from multiple regression 
mRNA.miRNA.list = get_miRNA_mRNA_multiple_regression(regression.data = regression.data, 
                                                     frequent.predictor.list = lasso.frequent.predictor.list)

# filter ceRNA pairs ---------------------------------------------------------------------
# create pairs
pair.dt = create_pair_dt(mRNA.miRNA.list = mRNA.miRNA.list)

# compute hypergeometric p-value to test the signifiance of number of shared miRNAs for each pair
pair.dt = get_hypergeometric_pvalue(pair.dt = pair.dt, 
                                    mRNA.miRNA.list = mRNA.miRNA.list)
# compute correlation for each pair
pair.dt = get_correlation(pair.dt = pair.dt, 
                          mRNA.expression = regression.data$mRNA)
# compute partial correlation and senstivity for each pair based on their shared mRNA regulators
pair.dt = get_partial_correlation(pair.dt = pair.dt, 
                                  mRNA.miRNA.list = mRNA.miRNA.list,
                                  mRNA.regression = regression.data$mRNA,
                                  miRNA.regression = regression.data$miRNA)

# compute empirical pvalue of sensitivity correlation for each pair
pair.dt = get_empirical_pvalue_sensitivity_correlation(pair.dt = pair.dt, 
                                                       mRNA.miRNA.list = mRNA.miRNA.list,
                                                       mRNA = regression.data$mRNA,
                                                       miRNA = regression.data$miRNA)




