# get lasso result with TF as potential regulators
get_lasso_result = function(regression.data){
  mRNA.targets = rownames(regression.data$mRNA)
  bt.size = ncol(regression.data$mRNA)

  if(file.exists(file.path)){file.remove(file.path)}
  if(file.exists(temp.path)){file.remove(temp.path)}
  
  vars = c("mRNA.targets", "regression.data", "bt.size", "TF.target.pairs")
  functions = c("getDF","getLassoCoefficients")
  varlist = c(vars, functions)
  
  cluster = parallel::makeCluster(parallel::detectCores()-1, outfile = temp.path)
  parallel::clusterExport(cl = cluster,varlist = varlist, envir = environment())
  invisible(parallel::clusterEvalQ(cl = cluster,expr = library(glmnet)))
  
  lasso.list = pbapply::pblapply(cl = cluster, X = 1:length(mRNA.targets), FUN = function(mRNA.index){
    mRNA = mRNA.targets[mRNA.index]
    bt.selected.coefs = lapply(1:100, function(bt.iter){
      # regression.df = getRegressionDF(regression.data = regression.data, mRNA = mRNA)
      regression.df = getDF(regression.data = regression.data, mRNA = mRNA)
      bt.regression.df = regression.df[sample(rownames(regression.df), size = bt.size, replace = T),]
      nonzero.lasso.coefs = getLassoCoefficients(regression.df = bt.regression.df)
      return(nonzero.lasso.coefs)
    })
  });
  
  names(lasso.list) = mRNA.targets
  stopCluster(cluster); gc()

  return(lasso.list)
}

# get regression df for each mRNA, TF included as candidate regulators
getDF = function(regression.data, mRNA){
  miRNA.candidates = as.character(regression.data$interactions[regression.data$interactions$mRNA == mRNA,c("miRNA")])
  tf.mRNA.interactions = regression.data$tf.mRNA.interactions
  TF.candidates = tf.mRNA.interactions$TF[which(tf.mRNA.interactions$target == mRNA)] 
  TF.candidates = TF.candidates[TF.candidates %in% rownames(regression.data$mRNA)]
  
  rna.expression.vector = regression.data$mRNA[mRNA,]
  rna.cna.vector = regression.data$cna[mRNA,]
  rna.methyl.vector = regression.data$methyl[mRNA,]
  predicted.TF.df = regression.data$mRNA[TF.candidates,]
  predicted.miRNA.df = regression.data$miRNA[miRNA.candidates,]
  
  regression.df = as.data.frame(scale(cbind(t(rna.expression.vector),t(rna.cna.vector),t(rna.methyl.vector),
                                            t(predicted.miRNA.df),t(predicted.TF.df))))
  colnames(regression.df)[1:3] = c("mRNA","CNA","Methyl")
  regression.df$CNA[which(is.na(regression.df$CNA))] = 0
  regression.df$Methyl[which(is.na(regression.df$Methyl))] = 0
  return(regression.df)
}

# Input: nested list of lasso-selected predictors for each DE mRNA
# Output: number of times a predictor is selected for each DE mRNA
get_list_of_lasso_dt = function(lasso.list){
  cluster = parallel::makeCluster(detectCores()-1)
  varlist = c("lasso.list")
  parallel::clusterExport(cl = cluster, varlist = varlist, envir = environment())
  invisible(parallel::clusterEvalQ(cl = cluster,expr = library(data.table)))
  # convert make each element in the list to be a data table
  lasso.dt.list = pbapply::pbsapply(cl = cluster, X = lasso.list, FUN = function(binders){
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
  return(lasso.dt.list)
}

# Input: number of times a predictor is selected for each DE mRNA
# Output: a predictor is kept if it is selected more than user-given threshold (default 70)
get_frequently_selected_predictors = function(lasso.dt.list,threshold){
  cluster = parallel::makeCluster(detectCores()-1)
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
  return(lasso.frequent.predictors)
}

get_miRNA_mRNA_multiple_regression = function(regression.data, frequent.predictor.list){
  cluster = makeCluster(detectCores()-1)
  varlist = c("regression.data","frequent.predictor.list","getRegressionDF","getOLSResultList")
  parallel::clusterExport(cl = cluster,varlist = varlist, envir = environment())
  invisible(parallel::clusterEvalQ(cl = cluster,expr = library(data.table)))
  mRNA.names = names(frequent.predictor.list)
  lasso.list = pbapply::pblapply(cl = cluster, X = 1:length(frequent.predictor.list), function(index){
    selected.mRNA = names(frequent.predictor.list[index])
    predictors = as.character(frequent.predictor.list[[index]]$name)
    normal.result = getOLSResultList(original.regression.df = getRegressionDF(regression.data,selected.mRNA),
                                     mRNA.target = selected.mRNA, predictors = predictors)
  }); stopCluster(cluster)
  names(lasso.list) = mRNA.names
  
  lasso.list = final.lasso.list[-which(sapply(final.lasso.list, function(item) (is.null(item))))]
  lasso.list = final.lasso.list[sort(names(final.lasso.list))]
  return(lasso.list)
}

getOLSResultList = function(original.regression.df, mRNA.target, predictors){
  if (length(predictors) == 0){
    returned.list = list(mRNA.expression = expression.matrix, seleted.miRs = NULL)
    names(returned.list)[2] = mRNA.target
    return(NULL)
  } else {
    predictors = colnames(original.regression.df)[which(colnames(original.regression.df) %in% predictors)]
    predictors = paste("`",predictors,"`",sep = "")
    # perform linear regression of the mRNA on the predicted miRNAs
    regression.formula = as.formula(paste("mRNA ~ ", paste(predictors, collapse= "+",sep=" ")))
    fit = lm(formula = regression.formula, data = original.regression.df)
    # get coefficients
    coefs = fit$coefficients[-1]; 
    pvalues = summary(fit)$coefficients[,4][-1]
    names(coefs) = gsub(names(coefs), pattern = "`",replacement = "") 
    names(pvalues) = gsub(names(pvalues), pattern = "`",replacement = "") 
    selected.index = which(pvalues<0.05)
    coefs = coefs[selected.index] 
    # extract non.miRs factors (Methyl and CNA), pos.miRs, and neg.miRs
    non.miR.coefs = coefs[-grep("hsa", names(coefs))]
    miR.coefs = coefs[grep("hsa", names(coefs))]
    pos.miR.coefs = miR.coefs[which(miR.coefs > 0)]
    neg.miR.coefs = miR.coefs[which(miR.coefs < 0)]
    non.neg.miR.factors = c(non.miR.coefs, pos.miR.coefs)
    # if there is no negative microRNA coefficients
    if (length(neg.miR.coefs) == 0){
      return(NULL)
    }else{
      return(names(neg.miR.coefs))
    }  
  }
}


# Input: list of miRNAs regulators for each mRNA
# Output: table of mRNA pairs such that the mRNAs in a pair share at least one miRNA
create_pair_dt = function(mRNA.miRNA.list){
  pair.dt = as.data.table(t(combn(names(mRNA.miRNA.list),2)))
  colnames(pair.dt) = c("RNAi", "RNAj"); setkeyv(pair.dt, c("RNAi","RNAj"))
  return(pair.dt)
}

# Input: table of mRNA pairs such that the mRNAs in a pair share at least one miRNA
# Output: compute correlation and correlation's p-value for each mRNA pairs
get_correlation = function(pair.dt, mRNA.expression){
  corMatrix = WGCNA::corFast(t(mRNA.expression),use = "pairwise.complete.obs")
  
  cluster = parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl = cluster,varlist = c("pair.dt", "corMatrix"),envir = environment())
  
  pair.dt$cor = pbapply::pbsapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(index){
    RNAi = pair.dt$RNAi[index]; RNAj = pair.dt$RNAj[index]
    return(corMatrix[RNAi,RNAj])
  }); stopCluster(cluster)
  
  return(pair.dt)
}

# Input: table of mRNA pairs such that the mRNAs in a pair share at least one miRNA
# Output: compute partial correlation and partial correlation's p-value for each mRNA pairs
get_partial_correlation = function(mRNA.regression, miRNA.regression, mRNA.miRNA.list, pair.dt){
  
  cluster = parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl = cluster, 
                          varlist = c("mRNA.regression","miRNA.regression","mRNA.miRNA.list", "pair.dt"),
                          envir = environment())
  invisible(parallel::clusterEvalQ(cl = cluster,expr = library(bnlearn,stats)))
  
  # compute partial correlation
  # pc.list = pbapply::pblapply(cl = cluster,X = 1:nrow(pair.dt), FUN = function(index)
  time = proc.time()
  pc.list = pbapply::pblapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(index){
    RNAi = pair.dt$RNAi[index]; RNAj = pair.dt$RNAj[index]; ceRNAs = c(RNAi, RNAj)
    common.miRs = intersect(mRNA.miRNA.list[[RNAi]], mRNA.miRNA.list[[RNAj]])
    
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
}

# Input: table of mRNA pairs such that the mRNAs in a pair share at least one miRNA
# Output: compute hypergeometric p-value for each mRNA pairs
get_hypergeometric_pvalue = function(pair.dt, mRNA.miRNA.list){
  pop.size = length(unique(unlist(mRNA.miRNA.list))) # total number of miRNAs
  cluster = parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl = cluster,
                          varlist = c("mRNA.miRNA.list", "pair.dt","pop.size","getPvalueHypergeometric"),
                          envir = environment())
  invisible(parallel::clusterEvalQ(cl = cluster,expr = library(data.table)))
  
  pair.dt$hypergeometric.pvalue = pbapply::pbsapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(index){
    RNAi = pair.dt$RNAi[index]; RNAj = pair.dt$RNAj[index]
    miRs.RNAi = mRNA.miRNA.list[[RNAi]]
    miRs.RNAj = mRNA.miRNA.list[[RNAj]]
    success.pop.size = length(miRs.RNAi)
    sample.size = length(miRs.RNAj)
    success.sample.size = pair.dt$num.common.miRs[index]
    getPvalueHypergeometric(pop.size = pop.size, success.pop.size = success.pop.size,
                            sample.size = sample.size, success.sample.size = success.sample.size)
  })
  return(pair.dt)
}

# Input: table of mRNA pairs such that the mRNAs in a pair share at least one miRNA
# Output: compute emperical p-value for each mRNA pair's sensitivity correlation
get_empirical_pvalue_sensitivity_correlation = function(pair.dt, mRNA.miRNA.list, mRNA.expression, miRNA.expression){
  permutated.sensitivity.cor.list = pbapply::pblapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(row.index){
    RNAi = pair.dt$RNAi[row.index];  RNAj = pair.dt$RNAj[row.index]
    common.miRs = intersect(mRNA.miRNA.list[[RNAi]], mRNA.miRNA.list[[RNAj]])
    
    ceRNA.sensitivity.permuted = sapply(1:1000, function(bootstrap_iter){
      data.expression = getDataWithResampledMiRNAs(RNAi = RNAi, RNAj = RNAj,
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
  parallel::clusterExport(cl = cluster, varlist = c("permutated.sensitivity.cor.list","pair.dt"))
  pair.dt$sensitivity.cor.pvalue = pbapply::pbsapply(cl = cluster, X = 1:nrow(pair.dt), FUN = function(index){
    return(sum(permutated.sensitivity.cor.list[[index]] >= pair.dt$sensitivity.cor[index])/1000)
  })
  return(pair.dt)
}

# compute pvalue for over-enrichment test using hypergeometric test
getPvalueHypergeometric = function(pop.size, success.pop.size, sample.size, success.sample.size){
  return(phyper(success.sample.size - 1, success.pop.size, pop.size - success.pop.size, sample.size, lower.tail = F))
}

# get expression data by resampling miRNAs without replacement
getDataWithResampledMiRNAs = function(RNAi, RNAj, all.miRs, common.miRs, mRNA.regression.df, miRNA.regression.df){
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



