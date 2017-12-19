### --- Change directory and set miRanda directory ---- #####
# main.dir = "~/Desktop/New_ceRNA/Test" # for local machine 
main.dir = "~/Work/encRNA/Test" # for server
setwd(main.dir) 

### --- Load libraries ---- #####
suppressMessages(require(data.table))
suppressMessages(require(parallel))
suppressMessages(require(bnlearn))
suppressMessages(require(WGCNA))

source("helper_function/lasso_helper_functions.R")
source("helper_function/mre_functions.R")
source("helper_function/variable_selection_functions.R")

load("rda/Gene.Info/TF.target.rda")

### --- Global variables ---- #####
args = commandArgs(trailingOnly = T)
if (!exists("cancer.type")){
  cancer.type = ifelse(length(args) !=0 , args[1], "HNSC")
  cancer.dir = paste("processed",cancer.type,sep=".")
}
save.path = paste("rda",cancer.dir,"variable_selection_output",sep="/")

# load: variable.selection.data
load(paste("rda",cancer.dir,"variable_selection_data","variable.selection.data.rda",sep="/"))
tumor.regression.data = subsetDataByType(list.data = variable.selection.data, sample.type = "Tumor")
# tumor.regression.data = subsetDataByType(list.data = lumA.variable.selection.data, sample.type = "Tumor")

# load: lasso.list depending on if lasso was performed with or without TFs
lasso.with.TF = T
if(lasso.with.TF){
  load(paste("rda",cancer.dir,"variable_selection_output","lasso2.result.rda",sep="/"))
  save.object.name = "july22.lasso2.ceRNA.result.rda"
}else{
  load(paste("rda",cancer.dir,"variable_selection_output","lasso.result.rda",sep="/"))
  save.object.name = "ceRNA.result.rda"
}

### --- Variables need to extract ceRNA ----------- #####
if (exists("lasso.list") & exists("tumor.regression.data")){
  lasso.list = lasso.list  
  regression.data = tumor.regression.data
}

### --- 1. Get list of predictor summary ---- #####
cat(paste("Get list of predictor summary","\n"))
source("code_variable_selection/ceRNA_helper/get_list_of_lasso_dt.R")
# output: lasso.dt.list

### --- 2. Keep predictors selected more than a threshold ---- #####
cat(paste("Keep predictors selected more than a threshold","\n"))
threshold = 70
source("code_variable_selection/ceRNA_helper/get_frequently_selected_predictors.R")
# output: lasso.frequent.predictors
assign(paste(tolower(cancer.type), "lasso.frequent.predictors", sep="."), value = lasso.frequent.predictors)


### --- 3. Get list of predictors and residual expression matrix ---- #####
cat(paste("Get list of predictors and residual expression matrix","\n"))
source("code_variable_selection/ceRNA_helper/get_final_list_and_residual_expression.R")
# output: lasso.result

### --- 4. Create pair dt ---- #####
cat(paste("Create pair dt","\n"))
source("code_variable_selection/ceRNA_helper/create_pair.R")
# output: pair.dt

### --- 5. Get number of common miRs for each pair ---- #####
cat(paste("Get number of common miRs for each ceRNA pair","\n"))
source("code_variable_selection/ceRNA_helper/get_common_miRs.R")
# output: pair.dt

### --- 6. Compute correlation and its pvalue for each ceRNA pair ---- #####
cat(paste("Compute correlation and its pvalue for each ceRNA pair","\n"))
source("code_variable_selection/ceRNA_helper/get_correlation.R")
# output: pair.dt

### --- 7. Compute partial correlation and its pvalue for each ceRNA pair ---- #####
cat(paste("Compute partial correlation and its pvalue for each ceRNA pair","\n"))
source("code_variable_selection/ceRNA_helper/get_partial_correlation.R")
# output: pair.dt

### --- 8. Compute hypergeometric pvalue for each ceRNA pair ---- #####
cat(paste("Compute hypergeometric pvalue for each ceRNA pair","\n"))
source("code_variable_selection/ceRNA_helper/get_hypergeometric_pvalue.R")
# output: pair.dt

# ### --- 9. Filter pairs by correlation-based filtering ---- #####
# cat(paste("Filter pairs by correlation-based filtering","\n"))
# cor.threshold = 0.5; adj.pvalue.threshold = 0.05
# source("code_variable_selection/ceRNA_helper/filter_pair_dt.R")
# # output: filtered.pair.dt

save(lasso.dt.list, lasso.result, pair.dt,
     file = paste(save.path, save.object.name,sep="/"))

### --- Get resampling-based sensitivity correlation pvalue ---- #####
load(paste(save.path, save.object.name,sep="/"))
cat(paste("Get resampling-based sensitivity correlation pvalue","\n"))
source("code_variable_selection/ceRNA_helper/get_sensitivity_correlation_emperical_pvalues.R")
# output: filtered.pair.dt

save(lasso.dt.list, lasso.result, pair.dt, filtered.pair.dt,
     file = paste(save.path, save.object.name,sep="/"))
















