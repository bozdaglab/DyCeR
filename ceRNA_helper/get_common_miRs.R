### ----- variables  ------#####
lasso.list = lasso.result$final.lasso.list
pop.size = length(unique(unlist(lasso.list))) # total number of miRNAs

### ----- analysis  ------#####
cluster = parallel::makeCluster(parallel::detectCores()-1)
parallel::clusterExport(cl = cluster,
                        varlist = c("lasso.list", "pair.dt","pop.size","getPvalueHypergeometric"),
                        envir = environment())
invisible(parallel::clusterEvalQ(cl = cluster,expr = library(data.table)))

pair.dt$num.common.miRs = pbapply::pbsapply(cl = cluster, X = 1:nrow(pair.dt), function(index){
  RNAi = pair.dt$RNAi[index]; RNAj = pair.dt$RNAj[index]
  miRs.RNAi = lasso.list[[RNAi]]
  miRs.RNAj = lasso.list[[RNAj]]
  success.pop.size = length(miRs.RNAi)
  sample.size = length(miRs.RNAj)
  success.sample.size = length(intersect(miRs.RNAi, miRs.RNAj))
  return(success.sample.size)
}); stopCluster(cluster)
pair.dt = pair.dt[-which(pair.dt$num.common.miRs < 1),,]

### ----- cleanup  ------#####
rm(cluster, lasso.list, pop.size)