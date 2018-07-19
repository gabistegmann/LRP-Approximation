lrp2_ALL = function (method, nlme.model = NULL, mvpartParams = NULL, 
                     data, start, group, rPartFormula, weight = NULL, R = NULL,
                     control.rp = rpart.control(cp = .01, xval = 10, minsplit = 17),
                     control.ct = ctree_control(minprob = 0.01),
                     control.et = evtree.control()) 
{
  
  groupingName = attr(terms(splitFormula(group, "~")[[1]]), 
                      "term.labels")
  responseName = attr(terms(getResponseFormula(nlme.model)), 
                      "term.labels")
  groupingFactor = data[, names(data) == groupingName]
  
  
  # Extract individual coeffs
  
  NLMEmodel = nlme(model = nlme.model, fixed = mvpartParams, 
                   data = data, random = mvpartParams, correlation = R, 
                   na.action = na.omit, start = start, group = group)
  
  data.2 = data
  data.2$id = groupingFactor
  
  data.2 =    data.2[!duplicated(data.2$id),]
  
  rand.coeffs = as.data.frame(NLMEmodel$coefficients$random[[1]])
  
  #scale coeffs
  stdd.rand.coef = scale(rand.coeffs)
  
  
  model = list()
  
  model$NLME = NLMEmodel
  
  #######
  # EFA #
  #######
  
  EFA.fit <- fa(r = stdd.rand.coef,
                nfactors = 1, fm = "ml")
  
  EFA.scores = as.data.frame(factor.scores(stdd.rand.coef, 
                                           EFA.fit, 
                                           method = "Thurstone")$scores)
  
  
  colnames(EFA.scores) = c("EFA.scores")
  
  EFA.scores$id = rownames(EFA.scores)
  
  # Merge EFA and data
  data.3 = merge(data.2, EFA.scores, by = "id")
  
 
  #CART

  model$EFA.rpart = list()
  

  model.rpart.1 = rpart(as.formula(paste("EFA.scores",
                                       rPartFormula)), 
                      data = data.3, 
                      control = control.rp)
  
  minCP = model.rpart.1$cptable[which.min(model.rpart.1$cptable[,4]),1]
  
  model.rpart = prune.rpart(model.rpart.1, cp = minCP)
  
  model$EFA.rpart$rpart_out <- model.rpart
  
  
  model$EFA.rpart$leaf_node <- model.rpart$where
  summary = list()
  
  for (j in 1:length(table(model.rpart$where))) {
    
    
    
    id <- names(table(model.rpart$where))[j] == model.rpart$where
    
    
    coeffs=NLMEmodel$coefficients$fixed
    
    nparams = length(coeffs)
    
    for(beta in 1:nparams){
      coeffs[beta]   =   NLMEmodel$coefficients$fixed[[beta]] + mean(rand.coeffs[id, beta])
    }
    
    
    summary[[as.numeric(names(table(model.rpart$where)))[j]]] <- coeffs
    
  }
  model$EFA.rpart$summary <- summary
  model$EFA.rpart$fixed_effects <- coeffs

  
  #CTREE
  model.ctree = ctree(as.formula(paste("EFA.scores",
                                       rPartFormula)), 
                      data = data.3, 
                      control = control.ct)
  
  
  model$EFA.ctree$ctree_out <- model.ctree
  
  
  summary = list()
  
  for (j in 1:length(model.ctree)) {
    
    if(j %in% nodeids(model.ctree, terminal = TRUE)){
      coef.node =rand.coeffs[rownames(model.ctree[j]$data),]
      
      coeffs=NLMEmodel$coefficients$fixed
      
      nparams = length(coeffs)
      
      for(beta in 1:nparams){
        coeffs[beta]   =   NLMEmodel$coefficients$fixed[[beta]] + mean(coef.node[, beta])
      }
      
      
      summary[[j]] <- coeffs
    }
    
  }
  model$EFA.ctree$summary <- summary
  
  
  
  #EVTREE
  
  model$EFA.evtree = list()
  
  model.evtree = evtree(as.formula(paste("EFA.scores",
                                         rPartFormula)), 
                        data = data.3, 
                        control = control.et)
  
  
  model$EFA.evtree$evtree_out <- model.evtree
  
  
  summary = list()
  
  for (j in 1:length(model.evtree)) {
    
    if(j %in% nodeids(model.evtree, terminal = TRUE)){
      coef.node =rand.coeffs[rownames(model.evtree[j]$data),]
      
      coeffs=NLMEmodel$coefficients$fixed
      
      nparams = length(coeffs)
      
      for(beta in 1:nparams){
        coeffs[beta]   =   NLMEmodel$coefficients$fixed[[beta]] + mean(coef.node[, beta])
      }
      
      
      summary[[j]] <- coeffs
    }
    
  }
  model$EFA.evtree$summary <- summary
  
  
  #######
  # PCA #
  #######
  
  #PCA
  PC.fit = prcomp(stdd.rand.coef)
  PC.scores = as.data.frame(PC.fit$x[,1])
  colnames(PC.scores) = c("PC.scores")
  
  PC.scores$id = rownames(PC.scores)
  
  # Merge PCA and data
  data.3 = merge(data.2, PC.scores, by = "id")
  
  
  #CART
  model$PCA.rpart = list()
  
  model.rpart.1 = rpart(as.formula(paste("PC.scores",
                                         rPartFormula)), 
                        data = data.3, 
                        control = control.rp)
  
  minCP = model.rpart.1$cptable[which.min(model.rpart.1$cptable[,4]),1]
  
  model.rpart = prune.rpart(model.rpart.1, cp = minCP)
  
  model$PCA.rpart$rpart_out <- model.rpart
  
  
  model$PCA.rpart$leaf_node <- model.rpart$where
  summary = list()
  
  for (j in 1:length(table(model.rpart$where))) {
    
    
    
    id <- names(table(model.rpart$where))[j] == model.rpart$where
    
    
    coeffs=NLMEmodel$coefficients$fixed
    
    nparams = length(coeffs)
    
    for(beta in 1:nparams){
      coeffs[beta]   =   NLMEmodel$coefficients$fixed[[beta]] + mean(rand.coeffs[id, beta])
    }
    
    
    summary[[as.numeric(names(table(model.rpart$where)))[j]]] <- coeffs
    
  }
  model$PCA.rpart$summary <- summary
  
  
  #CTREE
  model$PCA.ctree = list()
  
  model.ctree = ctree(as.formula(paste("PC.scores",
                                       rPartFormula)), 
                      data = data.3, 
                      control = control.ct)
  
  
  model$PCA.ctree$ctree_out <- model.ctree
  
  
  model$PCA.ctree$leaf_node <- model.ctree$where
  summary = list()
  
  for (j in 1:length(model.ctree)) {
    
    if(j %in% nodeids(model.ctree, terminal = TRUE)){
      coef.node =rand.coeffs[rownames(model.ctree[j]$data),]
      
      coeffs=NLMEmodel$coefficients$fixed
      
      nparams = length(coeffs)
      
      for(beta in 1:nparams){
        coeffs[beta]   =   NLMEmodel$coefficients$fixed[[beta]] + mean(coef.node[, beta])
      }
      
      
      summary[[j]] <- coeffs
    }
    
  }
  model$PCA.ctree$summary <- summary
  
  
  # EVTREE
  
  model$PCA.evtree = list()
  
  model.evtree = evtree(as.formula(paste("PC.scores",
                                         rPartFormula)), 
                        data = data.3, 
                        control = control.et)
  
  
  model$PCA.evtree$evtree_out <- model.evtree
  
  
  model$PCA.evtree$leaf_node <- model.evtree$where
  summary = list()
  
  for (j in 1:length(model.evtree)) {
    
    if(j %in% nodeids(model.evtree, terminal = TRUE)){
      coef.node =rand.coeffs[rownames(model.evtree[j]$data),]
      
      coeffs=NLMEmodel$coefficients$fixed
      
      nparams = length(coeffs)
      
      for(beta in 1:nparams){
        coeffs[beta]   =   NLMEmodel$coefficients$fixed[[beta]] + mean(coef.node[, beta])
      }
      
      
      summary[[j]] <- coeffs
    }
    
  }
  model$PCA.evtree$summary <- summary
  
  model$EFA.PCA.cor = abs(cor(EFA.scores$EFA.scores, PC.scores$PC.scores))
  
    return(model)
}
