### --- Change directory and set miRanda directory ---- #####
# main.dir = "~/Desktop/New_ceRNA/Test" # for local machine 
main.dir = "~/Work/encRNA/Test" # for server
setwd(main.dir) 

### --- Load libraries ---- #####
suppressMessages(require(data.table))
suppressMessages(require(parallel))
suppressMessages(require(glmnet))
source("helper_function/lasso_helper_functions.R")
load("rda/Gene.Info/TF.target.rda")

### --- Global variables ---- #####
args = commandArgs(trailingOnly = T)
if (!exists("cancer.type")){
  cancer.type = ifelse(length(args) !=0 , args[1], "BRCA")
  cancer.dir = paste("processed",cancer.type,sep=".")
}

save.path = paste("rda",cancer.dir,"variable_selection_output",sep="/")
if(!file.exists(save.path)){
  dir.create(path=file.path(paste(main.dir,"rda",cancer.dir,sep="/"),"variable_selection_output"))
}

tumor.debug.file.name = paste("Tumor.Lasso",cancer.type,Sys.Date(),"txt",sep=".")
# tumor.debug.file.name = paste("Tumor.LumA.Lasso",cancer.type,Sys.Date(),"txt",sep=".")

data.to.load = "variable.selection.data.rda"
load(paste("rda",cancer.dir,"variable_selection_data",data.to.load,sep="/"))

#### --------- main codes ---------------------------------------------
tumor.regression.data = subsetDataByType(list.data=variable.selection.data, sample.type="Tumor")

mRNA = tumor.regression.data$mRNA
miRNA = tumor.regression.data$miRNA
methyl = tumor.regression.data$methyl
cna = tumor.regression.data$cna

mRNA.targets = rownames(tumor.regression.data$mRNA)
bt.size  = ncol(tumor.regression.data$mRNA)

lasso.list = get_lasso_result(mRNA.targets = mRNA.targets,
                              regression.data = tumor.regression.data,
                              TF.target.pairs = TF.target,
                              bt.size = bt.size,
                              outfile.name = tumor.debug.file.name)

save(lasso.list, file = paste(save.path,"lasso2.result.rda",sep="/"))
