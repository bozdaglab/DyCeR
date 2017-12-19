### ----- variables  ------#####
cat(paste("    compute correlation and pvalue","\n"))
#mRNA.expression = lasso.result$residual.expression
mRNA.expression = regression.data$mRNA
cor.info = WGCNA::corAndPvalue(t(mRNA.expression),use = "pairwise.complete.obs")
pair.dt = pair.dt

### ----- analysis  ------#####
cluster = parallel::makeCluster(parallel::detectCores()-1)
parallel::clusterExport(cl = cluster,varlist = c("pair.dt", "cor.info"),envir = environment())

cat(paste("   assign correlation","\n"))
pair.dt$cor = pbapply::pbsapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(index){
  RNAi = pair.dt$RNAi[index]; RNAj = pair.dt$RNAj[index]
  return(cor.info$cor[RNAi,RNAj])
})

cat(paste("   assign correlation pvalues","\n"))
pair.dt$cor.pvalue = pbapply::pbsapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(index){
  RNAi = pair.dt$RNAi[index]; RNAj = pair.dt$RNAj[index]
  return(cor.info$p[RNAi,RNAj])
})
stopCluster(cluster)

### ----- clean-up  ------#####
rm(mRNA.expression, cor.info)