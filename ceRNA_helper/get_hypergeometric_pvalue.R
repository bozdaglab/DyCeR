### ----- variables  ------#####
lasso.list = lasso.result$final.lasso.list
pair.dt = pair.dt

### ----- analysis  ------#####
pop.size = length(unique(unlist(lasso.list))) # total number of miRNAs
cluster = parallel::makeCluster(parallel::detectCores()-1)
parallel::clusterExport(cl = cluster,
                        varlist = c("lasso.list", "pair.dt","pop.size","getPvalueHypergeometric"),
                        envir = environment())
invisible(parallel::clusterEvalQ(cl = cluster,expr = library(data.table)))

pair.dt$hypergeometric.pvalue = pbapply::pbsapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(index){
  #index = 10000
  RNAi = pair.dt$RNAi[index]; RNAj = pair.dt$RNAj[index]
  miRs.RNAi = lasso.list[[RNAi]]
  miRs.RNAj = lasso.list[[RNAj]]
  success.pop.size = length(miRs.RNAi)
  sample.size = length(miRs.RNAj)
  success.sample.size = pair.dt$num.common.miRs[index]
  getPvalueHypergeometric(pop.size = pop.size, success.pop.size = success.pop.size,
                          sample.size = sample.size, success.sample.size = success.sample.size)
})

### ----- clean up  ------#####
rm(lasso.list)