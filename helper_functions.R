# get lasso result with TF as potential regulators
get_lasso_result = function(regression.data){
  mRNA.targets = rownames(regression.data$mRNA)
  bt.size = ncol(regression.data$mRNA)
  
  vars = c("mRNA.targets", "regression.data", "bt.size")
  functions = c("getDF","getLassoCoefficients")
  varlist = c(vars, functions)
  
  cluster = parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl = cluster,varlist = varlist, envir = environment())
  invisible(parallel::clusterEvalQ(cl = cluster,expr = library(glmnet)))
  
  lasso.list = pbapply::pblapply(cl = cluster, X = 1:length(mRNA.targets), FUN = function(mRNA.index){
    mRNA = mRNA.targets[mRNA.index]
    bt.selected.coefs = lapply(1:100, function(bt.iter){
      regression.df = getDF(regression.data = regression.data, mRNA = mRNA)
      bt.regression.df = regression.df[sample(rownames(regression.df), size = bt.size, replace = T),]
      nonzero.lasso.coefs = getLassoCoefficients(regression.df = bt.regression.df)
      return(nonzero.lasso.coefs)
    })
  });
  
  names(lasso.list) = mRNA.targets
  stopCluster(cluster); gc()
  
  lasso.list = get_list_of_lasso_dt(lasso.list)
  return(lasso.list)
}

# get regression df for each mRNA, TF included as candidate regulators
getDF = function(regression.data, mRNA){
  miRNA.target.interactions = regression.data$miRNA.target.interactions
  miRNA.candidates = as.character(miRNA.target.interactions[miRNA.target.interactions$target == mRNA,c("miRNA")])
  tf.target.interactions = regression.data$tf.target.interactions
  TF.candidates = tf.target.interactions$tf[which(tf.target.interactions$target == mRNA)] 
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
  lasso.dt.list = lapply(lasso.list, function(predictors){
    sort(table(unlist(sapply(predictors, function(binders) names(binders)))),decreasing = T)
  })
  no.predictor.indices = which(sapply(lasso.dt.list, function(predictors){return(is.null(predictors))}))
  if (length(no.predictor.indices) > 0){
    lasso.dt.list = lasso.dt.list[-no.predictor.indices]
  }
  return(lasso.dt.list)
}

# Input: number of times a predictor is selected for each DE mRNA
# Output: a predictor is kept if it is selected more than user-given threshold (default 70)
get_frequently_selected_predictors = function(lasso.list,threshold){
  selected.predictors = lapply(lasso.list, function(predictors){
    kept = which(predictors > threshold)
    if (length(kept) == 0){
      return(NULL)
    }else{
      return(names(predictors)[kept])
    }
  })
  no.predictor.indices = which(sapply(selected.predictors, function(predictors){return(is.null(predictors))}))
  if (length(no.predictor.indices) > 0){
    selected.predictors = selected.predictors[-no.predictor.indices]
  }
  return(selected.predictors)
}

get_miRNA_mRNA_multiple_regression = function(regression.data, frequent.predictor.list){
  cluster = makeCluster(detectCores()-1)
  varlist = c("regression.data","frequent.predictor.list","getOLSResultList")
  parallel::clusterExport(cl = cluster,varlist = varlist, envir = environment())
  invisible(parallel::clusterEvalQ(cl = cluster,expr = library(data.table)))
  mRNA.names = names(frequent.predictor.list)
  mRNA.miRNA.list = pbapply::pblapply(X = 1:length(frequent.predictor.list), function(index){
    selected.mRNA = names(frequent.predictor.list[index])
    predictors = as.character(frequent.predictor.list[[index]])
    original.regression.df = getDF(regression.data = regression.data,mRNA = selected.mRNA)
    normal.result = getOLSResultList(original.regression.df = original.regression.df,
                                     mRNA.target = selected.mRNA, 
                                     predictors = predictors)
  })
  names(mRNA.miRNA.list) = mRNA.names
  stopCluster(cluster)
  
  names(lasso.list) = mRNA.names
  no.predictor.indices = which(sapply(mRNA.miRNA.list, function(item) (is.null(item))))
  if (length(no.predictor.indices) > 0){
    mRNA.miRNA.list = mRNA.miRNA.list[-no.predictor.indices]
  }
  if (length(mRNA.miRNA.list) > 0){
    mRNA.miRNA.list = mRNA.miRNA.list[sort(names(mRNA.miRNA.list))]  
    return(mRNA.miRNA.list)
  }else{
    cat("No miRNA regulator found for any mRNA")
    return(NULL)
  }
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
  pair.dt$num.common.miRNAs = apply(pair.dt, 1, function(pair){
    RNAi = pair[1]
    RNAj = pair[2]
    num.common.miRs = length(intersect(mRNA.miRNA.list[[RNAi]],mRNA.miRNA.list[[RNAj]]))
    return(num.common.miRs)
  })
  pair.dt = pair.dt[pair.dt$num.common.miRNAs>0]
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
    success.sample.size = pair.dt$num.common.miRNAs[index]
    getPvalueHypergeometric(pop.size = pop.size, success.pop.size = success.pop.size,
                            sample.size = sample.size, success.sample.size = success.sample.size)
  })
  return(pair.dt)
}

# Input: table of mRNA pairs such that the mRNAs in a pair share at least one miRNA
# Output: compute emperical p-value for each mRNA pair's sensitivity correlation
get_empirical_pvalue_sensitivity_correlation = function(pair.dt, mRNA.miRNA.list, mRNA, miRNA){
  permutated.sensitivity.cor.list = pbapply::pblapply(X = 1:nrow(pair.dt), FUN = function(row.index){
    RNAi = pair.dt$RNAi[row.index];  RNAj = pair.dt$RNAj[row.index]
    
    ceRNA.sensitivity.permuted = sapply(1:1000, function(bootstrap_iter){
      data.expression = getDataWithResampledMiRNAs(RNAi = RNAi, RNAj = RNAj,
                                                   num.common.miRNAs = pair.dt$num.common.miRNAs[row.index],
                                                   mRNA.regression.df =  mRNA, 
                                                   miRNA.regression.df = miRNA)
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
getDataWithResampledMiRNAs = function(RNAi, RNAj, num.common.miRNAs, mRNA, miRNA){
  RNAi.expression = t(as.matrix(mRNA[RNAi,])); colnames(RNAi.expression) = RNAi
  RNAj.expression = t(as.matrix(mRNA[RNAj,])); colnames(RNAj.expression) = RNAj
  # construct data expression
  random.miRs = sample(rownames(miRNA), size = num.common.miRNAs, replace = F)
  permutated.miRs.expression = as.matrix(miRNA.regression.df[random.miRs, ])
  if (length(random.miRs) == 1){
    rownames(permutated.miRs.expression) = random.miRs
  }
  permutated.miRs.expression = t(miRNA.regression.df[random.miRs,])
  data.expression = as.data.frame(scale(cbind(RNAi.expression, RNAj.expression, permutated.miRs.expression)))
  return(data.expression)
}


getLassoCoefficients = function(regression.df){
  x = model.matrix(mRNA~.,regression.df)[,-1]; y=regression.df$mRNA
  # perform Lasso
  cv.out= cv.glmnet(x,y,alpha=1,nfolds = 10)
  coefs = coef.cv.glmnet(cv.out,s = cv.out$lambda.1se)[-1,]
  selected.coefs = coefs[coefs !=0]
  if (length(selected.coefs) == 0){
    return(NULL)
  }
  names(selected.coefs) = gsub(names(selected.coefs),pattern="`",replacement="")
  return(selected.coefs)
}



