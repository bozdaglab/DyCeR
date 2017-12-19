### ----- variables  ------#####
#mRNA.regression = lasso.result$residual.expression
mRNA.regression = regression.data$mRNA
miRNA.regression = regression.data$miRNA
lasso.list = lasso.result$final.lasso.list
pair.dt = pair.dt

### ----- analysis  ------#####
cluster = parallel::makeCluster(parallel::detectCores()-1)
parallel::clusterExport(cl = cluster, 
                        varlist = c("mRNA.regression","miRNA.regression","lasso.list", "pair.dt"),
                        envir = environment())
invisible(parallel::clusterEvalQ(cl = cluster,expr = library(bnlearn,stats)))

# compute partial correlation

pc.list = pbapply::pblapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(index){
  RNAi = pair.dt$RNAi[index]; RNAj = pair.dt$RNAj[index]; ceRNAs = c(RNAi, RNAj)
  common.miRs = intersect(lasso.list[[RNAi]], lasso.list[[RNAj]])
  
  #RNA.RNA.df = getDataForConditionalTest(regression.data = regression.data, ceRNAs = ceRNAs, common.miRs = common.miRs)
  mRNA.df = mRNA.regression[ceRNAs,] # since mRNA.name has only one element
  miRNA.df = miRNA.regression[common.miRs,]
  RNA.RNA.df = as.data.frame(scale(cbind(t(mRNA.df), t(miRNA.df))))
  
  ci.test.result = ci.test(x = colnames(RNA.RNA.df)[1], y = colnames(RNA.RNA.df)[2], z = colnames(RNA.RNA.df)[3:ncol(RNA.RNA.df)],
                           data = RNA.RNA.df,test = "cor")
  
  return(list(pcor = ci.test.result$statistic,
              pvalue = ci.test.result$p.value))
}); stopCluster(cluster)

# get partial correlation and p.value for the test
pair.dt$partial.corr = sapply(pc.list, function(result) result$pcor)
#pair.dt$pvalue.condTest = p.adjust(p = sapply(pc.list, function(result) result$pvalue),method = "BH")
pair.dt$sensitivity.cor = pair.dt$cor - pair.dt$partial.corr
pair.dt$partial.cor.pvalue = sapply(pc.list, function(result) result$pvalue);
return(pair.dt)

### ----- clean up  ------#####
rm(pc.list); rm(mRNA.regression); rm(miRNA.regression); rm(lasso.list)
