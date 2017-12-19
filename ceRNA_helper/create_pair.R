### ----- variables  ------#####
lasso.list = lasso.result$final.lasso.list

### ----- analysis  ------#####
pair.dt = as.data.table(t(combn(names(lasso.list),2)))
colnames(pair.dt) = c("RNAi", "RNAj"); setkeyv(pair.dt, c("RNAi","RNAj"))

### ----- cleanup  ------#####
rm(lasso.list)