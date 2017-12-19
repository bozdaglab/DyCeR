### ----- variables  ------#####
lasso.dt.list = lasso.dt.list
threshold = threshold

### ----- analysis  ------#####
cluster = makeCluster(detectCores()-1)
parallel::clusterExport(cl = cluster,varlist = c("lasso.dt.list","threshold"),envir = environment())
invisible(parallel::clusterEvalQ(cl = cluster,expr = library(data.table)))

lasso.frequent.predictors = pbapply::pblapply(cl = cluster, X = lasso.dt.list, FUN = function(dt){
  dt = dt[which(dt$count>=threshold)]
  if(nrow(dt)==0) {return(NULL)}
  else {return(dt)}
})
stopCluster(cluster)
unkept.predictor.indicies = which(sapply(lasso.frequent.predictors, function(predictors){ is.null(predictors) }))
if (length(unkept.predictor.indicies) > 0){
  lasso.frequent.predictors = lasso.frequent.predictors[-unkept.predictor.indicies] 
}
### ----- clean up  ------#####
rm(cluster) ; rm(unkept.predictor.indicies)