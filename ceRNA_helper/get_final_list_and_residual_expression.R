### ----- variables  ------#####
regression.data = regression.data
lasso.list = lasso.list
TF.target = TF.target

### ----- analysis  ------#####
cluster = makeCluster(detectCores()-1)
vars = c("regression.data","lasso.frequent.predictors", "TF.target")
funcs =  c("getRegressionDF","getDF","getOLSResultList")
varlist = c(vars, funcs)
parallel::clusterExport(cl = cluster,varlist = varlist, envir = environment())
invisible(parallel::clusterEvalQ(cl = cluster,expr = library(data.table)))
mRNA.names = names(lasso.frequent.predictors)
lasso.list = pbapply::pblapply(cl = cluster, X = 1:length(lasso.frequent.predictors), function(index){
  selected.mRNA = names(lasso.frequent.predictors[index])
  predictors = as.character(lasso.frequent.predictors[[index]]$name)
  normal.result = getOLSResultList(original.regression.df = 
                                     getDF(regression.data = regression.data,mRNA = selected.mRNA,TF.target.pairs = TF.target), 
                                   mRNA.target = selected.mRNA, 
                                   predictors = predictors)
}); stopCluster(cluster)  
names(lasso.list) = mRNA.names

residual.expression = NULL
final.lasso.list = lapply(lasso.list, function(list.item){
  residual.expression <<- cbind(residual.expression, list.item[[1]])
  list.item[[2]]
})
residual.expression = t(residual.expression)
# remove item with no predictors
final.lasso.list = final.lasso.list[-which(sapply(final.lasso.list, function(item) (is.null(item))))]
final.lasso.list = final.lasso.list[sort(names(final.lasso.list))]
lasso.result = list(final.lasso.list=final.lasso.list, residual.expression=residual.expression)

### ----- clean up ------#####
rm(cluster, mRNA.names, residual.expression, final.lasso.list)

