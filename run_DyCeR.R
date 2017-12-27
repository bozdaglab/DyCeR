### --- Load libraries ---- #####
suppressMessages(require(data.table))
suppressMessages(require(parallel))
suppressMessages(require(bnlearn))
suppressMessages(require(WGCNA))

source("helper_functions.R")

# # to edit -------------------------------------------------------------------------------
# setwd("~/Desktop/New_ceRNA/Test")
# load("rda/Gene.Info/TF.target.rda")
# cancer.type = "BRCA"
# cancer.dir = paste("processed",cancer.type,sep=".")
# load(paste("rda",cancer.dir,"variable_selection_data","variable.selection.data.rda",sep="/"))
# subsetDataByType = function(list.data, sample.type){
#   samples = grep(colnames(list.data$mRNA), pattern = sample.type)
#   return(list(mRNA = list.data$mRNA[,samples],
#               cna = list.data$cna[,samples],
#               methyl = list.data$methyl[,samples],
#               miRNA = list.data$miRNA[,samples],
#               interactions = list.data$interactions))
# }
# tumor.regression.data = subsetDataByType(list.data = variable.selection.data, sample.type = "Tumor")
# mRNA = tumor.regression.data$mRNA
# miRNA = tumor.regression.data$miRNA
# methyl = tumor.regression.data$methyl
# cna = tumor.regression.data$cna

# run bootstrapped LASSO procedure for each DE mRNAs ------------------------------------
regression.data = list(mRNA=mRNA, 
                       miRNA=miRNA, 
                       methyl=methyl, 
                       cna=cna,
                       miRNA.mRNA.interactions=miRNA.mRNA.interactions,
                       tf.mRNA.interactions=tf.mRNA.interactions)
# perform the bootstrapped LASSO procedure
lasso.list = get_lasso_result(regression.data = regression.data)
# get list of predictor summary 
lasso.dt = get_list_of_lasso_dt(lasso.list = lasso.list)
# keep predictors selected more than a threshold
lasso.frequent.predictors = get_frequently_selected_predictors(lasso.dt.list = lasso.dt, threshold = 70)
# get list of predictors and residual expression matrix 
mRNA.miRNA.list = get_miRNA_mRNA_multiple_regression(regression.data = regression.data, 
                                                  frequent.predictor.list = lasso.frequent.predictors)

# filter ceRNA pairs --------------------------------------------------------------------
# create pairs
pair.dt = create_pair_dt(mRNA.miRNA.list = mRNA.miRNA.list)
# compute correlation for each pair
pair.dt = get_correlation(pair.dt = pair.dt, 
                          mRNA.expression = regression.data$mRNA)
# compute partial correlation and senstivity for each pair based on their shared mRNA regulators
pair.dt = get_partial_correlation(pair.dt = pair.dt, 
                                  mRNA.miRNA.list = mRNA.miRNA.list,
                                  mRNA.regression = regression$mRNA,
                                  miRNA.regression = regression$miRNA)
# compute empirical pvalue of sensitivity correlation for each pair
pair.dt = get_empirical_pvalue_sensitivity_correlation(pair.dt = pair.dt, 
                                                       mRNA.miRNA.list = mRNA.miRNA.list,
                                                       mRNA.regression = regression$mRNA,
                                                       miRNA.regression = regression$miRNA)
# compute hypergeometric p-value to test the signifiance of number of shared miRNAs for each pair
pair.dt = get_hypergeometric_pvalue(pair.dt = pair.dt, 
                                    mRNA.miRNA.list = mRNA.miRNA.list)



