### ----- function and variables ----------###
getPermutedData = function(RNAi, RNAj, all.miRs, common.miRs, mRNA.regression.df, miRNA.regression.df){
  RNAi.expression = as.matrix(mRNA.regression.df[RNAi,]); colnames(RNAi.expression) = RNAi
  RNAj.expression = as.matrix(mRNA.regression.df[RNAj,]); colnames(RNAj.expression) = RNAj
  # construct data expression
  random.miRs = sample(all.miRs, size = length(common.miRs), replace = F)
  permutated.miRs.expression = as.matrix(miRNA.regression.df[random.miRs, ])
  if (length(random.miRs) == 1){
    rownames(permutated.miRs.expression) = random.miRs
  }
  permutated.miRs.expression = t(miRNA.regression.df[random.miRs,])
  data.expression = as.data.frame(scale(cbind(RNAi.expression, RNAj.expression, permutated.miRs.expression)))
  return(data.expression)
}
 
miRNA.expression = regression.data$miRNA
all.miRs = rownames(miRNA.expression)
mRNA.expression = lasso.result$residual.expression
lasso.list = lasso.result$final.lasso.list
filtered.pair.dt = filtered.pair.dt

### ----- analysis ----------###
cluster = parallel::makeCluster(detectCores() - 1)
parallel::clusterExport(cl = cluster, varlist = c("all.miRs","filtered.pair.dt",
                                                  "lasso.list","mRNA.expression","miRNA.expression",
                                                  "getPermutedData"))
invisible(parallel::clusterEvalQ(cl = cluster,expr = library(bnlearn,stats)))

permutated.sensitivity.cor.list = pbapply::pblapply(cl = cluster, X = 1:nrow(filtered.pair.dt), FUN = function(row.index){
  RNAi = filtered.pair.dt$RNAi[row.index];  RNAj = filtered.pair.dt$RNAj[row.index]
  common.miRs = intersect(lasso.list[[RNAi]], lasso.list[[RNAj]])

  ceRNA.sensitivity.permuted = sapply(1:1000, function(bootstrap_iter){
    data.expression = getPermutedData(RNAi = RNAi, RNAj = RNAj,
                                      all.miRs = all.miRs, common.miRs = common.miRs,
                                      mRNA.regression.df =  mRNA.expression, miRNA.regression.df = miRNA.expression)
    # compute permutated correlation
    ceRNA.cor = as.numeric(WGCNA::corFast(data.expression[,1], data.expression[,2]))
    ceRNA.partial.cor = ci.test(x = colnames(data.expression)[1],
                                y = colnames(data.expression)[2],
                                z = colnames(data.expression)[3:ncol(data.expression)],
                                data = data.expression,test = "cor")$statistic
    return(unname(ceRNA.cor - ceRNA.partial.cor))
  })
})

cluster = parallel::makeCluster(detectCores() - 1)
parallel::clusterExport(cl = cluster, varlist = c("permutated.sensitivity.cor.list","filtered.pair.dt"))
filtered.pair.dt$sensitivity.cor.emperical.pvalue = pbapply::pbsapply(cl = cluster, X = 1:nrow(filtered.pair.dt), FUN = function(index){
  return(sum(permutated.sensitivity.cor.list[[index]] >= filtered.pair.dt$sensitivity.cor[index])/1000)
})

### ----- clean up ----------###
rm(miRNA.expression); rm(mRNA.expression); rm(all.miRs)







