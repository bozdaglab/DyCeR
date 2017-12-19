### ----- variables  ------#####
lasso.list = lasso.list

### ----- analysis  ------#####
cluster = parallel::makeCluster(detectCores()-1)
parallel::clusterExport(cl = cluster,varlist = c("lasso.list"),envir = environment())
invisible(parallel::clusterEvalQ(cl = cluster,expr = library(data.table)))

# convert make each element in the list to be a data table
lasso.dt.list = pbapply::pbsapply(cl = cluster, X = lasso.list, function(binders){
  predictors = unlist(binders)
  print(length(predictors))
  if (length(predictors) == 0) return(NULL)
  predictor.name = names(predictors)
  predictor.coef = unname(predictors)
  df = data.frame(name = predictor.name, coef = predictor.coef)
  dt = as.data.table(df)
  dt = dt[, .(count = .N, avg.coef = mean(coef)), by = .(name)]
  dt = dt[order(count,decreasing = T)]
  return(dt)
})
stopCluster(cluster)
lasso.dt.list = lasso.dt.list[-which(sapply(lasso.dt.list, function(predictors){return(is.null(predictors))}))]

### ----- clean up  ------#####
rm(cluster)

