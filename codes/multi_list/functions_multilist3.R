library(drpop)

risklim <- function(q1, q2, q12, yi, yj, risk_alpha, eps = 0.005){
  
  psiinvmat = rep(NA, 6)
  names(psiinvmat) = paste(rep(c("PI", "DR", "TMLE"), each = 2), c("l", "u"), sep = '.')
  varmat = psiinvmat
  
  for(lim in c('l', 'u')){
    if(lim == 'u'){
      delta = rep(risk_alpha, length(yi))
      delta = pmax(delta, (1-q2)*q12/(q1-q12)*q2)
    }else{
      delta = rep(1/risk_alpha, length(yi))
      delta = pmax(delta, (1-q2)*q12/(q1-q12)*q2)
    }
    gammainvhat = (delta*(q1 - q12) + q12)*q2/q12
    
    gammainvhat2 = gammainvhat
    gammainvhat2[gammainvhat2 < 1] = 1
    psiinvhat = mean(gammainvhat2, na.rm = TRUE)
    
    phihat = gammainvhat*(yj/q2 + (delta*(yi - yi*yj) + yi*yj)/(delta*(q1 - q12) + q12) - yi*yj/q12)
    phihat[gammainvhat < 1] = 1
    Qnphihat = mean(phihat, na.rm = TRUE)
    
    psiinvhatq = max(Qnphihat, 1)
    
    psiinvmat[paste(c('PI', "DR"), lim, sep = '.')] = c(psiinvhat, psiinvhatq)
    varmat[paste(c('PI', "DR"), lim, sep = '.')] = c(var(phihat), var(phihat))
    
    #print(c(psiinvhat, psiinvhatq))
    
    datmat = as.data.frame(cbind(yi, yj, yi*yj, q1 - q12, q2 - q12, q12, yi*(1 - yj), yj*(1 - yi)))
    datmat[,4:6] = cbind(apply(datmat[,4:6], 2, function(u) {return(pmin(pmax(u, eps), 1 - eps))}))
    colnames(datmat) = c("yi", "yj", "yij", "q10", "q02", "q12", "yi0", "y0j")
    
    epsilon_error = 0
    cnt = 0
    iter = 10; eps = 0.005; twolist = FALSE    
    while (abs(epsilon_error) > 0.00001){
      cnt = cnt + 1
      datmat_old = datmat
      if (cnt > iter){break}
      
      ########################### model 1 for q12
      dat1 = cbind(datmat$yij, logit(datmat$q12), delta*(datmat$q10 + datmat$q12)/datmat$q12
                   + delta*(datmat$q02 + datmat$q12)/datmat$q12
                   - delta*(datmat$q10 + datmat$q12)*(datmat$q02 + datmat$q12)/datmat$q12^2 + 1 - delta)
      colnames(dat1) = c("yij", "logitq12", "ratio")
      dat1 = as.data.frame(dat1)
      mod1 = try(glm(yij ~ -1 + offset(logitq12) + ratio
                     , family = binomial(link = logit), data = dat1, na.action = na.omit))
      if (class(mod1) != "try-error"){
        datmat[,"q12"] = predict(mod1, newdata = dat1, type = "response")
        
      }
      summary(cbind(datmat[,4:6], datmat_old[,4:6]))
      datmat$q12 = pmax(pmin(datmat$q12, 1), eps)
      
      ########################### model 2 for q1
      dat2 = cbind(datmat$yi0, logit(datmat$q10), delta*(datmat$q02 + datmat$q12)/datmat$q12)
      colnames(dat2) = c("yi0", "logitq10", "ratio")
      dat2 = as.data.frame(dat2)
      mod2 = try(glm(yi0 ~ -1 + offset(logitq10) + ratio, family = binomial(link = logit), data = dat2, na.action = na.omit))
      if (class(mod2) != "try-error"){
        datmat$q10 = predict(mod2, newdata = dat2, type = "response")
        datmat[,"q10"] = pmin(datmat[,"q10"], 1 - datmat$q12)
      }
      summary(cbind(datmat[,4:6], datmat_old[,4:6]))
      datmat$q10 = pmax(pmin(datmat$q10, 1), eps)
      
      ########################### model 3 for q2
      if (K > 2 | twolist == FALSE){
        dat3 = cbind(datmat$y0j, logit(datmat$q02), delta*(datmat$q10 + datmat$q12)/datmat$q12 + 1 - delta)
        colnames(dat3) = c("y0j", "logitq02", "ratio")
        dat3 = as.data.frame(dat3)
        mod3 = try(glm(y0j ~ -1 + offset(logitq02) + ratio, family = binomial(link = logit), data = dat3, na.action = na.omit))
        if (class(mod3) != "try-error"){
          datmat$q02 = predict(mod3, newdata = dat3, type = "response")
          datmat[,"q02"] = pmin(datmat[,"q02"], 1 - datmat$q10 - datmat$q12)
        }
      }else{
        mod3 = mod2
        datmat[,"q02"] = pmax(0, 1 - datmat$q10 - datmat$q12)
      }
      
      datmat$q02 = pmax(pmin(datmat$q02, 1), eps)
      
      epsilon_error = max(abs(c(mod2$coefficients, mod3$coefficients, mod1$coefficients)))
    }
    
    if(epsilon_error < 1){
      
      q12t = pmax(datmat$q12, eps)
      q1t = pmin(datmat$q12 + datmat$q10, 1)
      q2t = pmax(pmin(datmat$q12 + datmat$q02, 1 + q12t - q1t, 1), q12t/q1t)
      
      gammainvhat = (delta*(q1t - q12t) + q12t)*q2t/q12t
      psiinvhat = mean(gammainvhat, na.rm = TRUE)
      
      phihat = gammainvhat*(yj/q2t + (delta*(yi - yi*yj) + yi*yj)/(delta*(q1t - q12t) + q12t) - yi*yj/q12t)
      
      Qnphihat = mean(phihat, na.rm = TRUE)
      
      psiinvmat[paste0("TMLE.", lim)] = psiinvhat
      
      varmat[paste0("TMLE.", lim)] = var(phihat)
    }
  }
  return(result = list(psiinvmat = psiinvmat, varmat = varmat))
}

phihat_risk = function(List1, List2, K = 2, i = 1, j = 2, func = "logit", twolist = TRUE, iter = 50, risk_max = 1, length_alpha = 10, actual, alpha, omega, eps = 0.005){
  
  l = ncol(List1) - K
  
  conforminglists = apply(List1[,1:K], 2, function(col){return(setequal(col, c(0,1)))})
  if(sum(conforminglists) < 2){
    print("Data is not in the required format or lists are degenerate.")
    return(NULL)
  }
  
  if(sum(conforminglists) < K){
    cat("Lists ", which(conforminglists == FALSE), " are not in the required format.\n")
  }
  
  delta_vec = c(1/risk_max, risk_max)
  
  if(l == 0){
    #renaming the columns of List_matrix for ease of use
    colnames(List1) = c(paste("L", 1:K, sep = ''))
    colnames(List2 = colnames(List1))
    q1 = mean(List_matrix[,i])
    q2 = mean(List_matrix[,j])
    q12 = mean(List_matrix[,i]*List_matrix[,j])
    
    for(delta in delta_vec){
      
      psiinv_summary[paste(i, ", ", j, sep = ''),] = (delta*q1 + (1 - delta)*q12)*q2/q12
      var_summary[paste(i, ", ", j, sep = ''),] = q1*q2*(q1*q2 - q12)*(1 - q12)/q12^3/N
    }
    
  }else{
    #renaming the columns of List_matrix for ease of use
    colnames(List1) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List1) - K), sep = ''))
    colnames(List2) = colnames(List1)
    
    yi = List2[,paste("L", i, sep = '')]
    yj = List2[,paste("L", j, sep = '')]
    
    psiinvmat = numeric(0)
    varmat = numeric(0)
    if(mean(List1[,i]*List1[,j]) > eps) {
      
      #colsubset = stringr::str_subset(colnames(psiinv_summary), func)
      if(actual == FALSE){
        qhat = try(get(paste0("qhat_", func))(List1, List2, K, i, j, eps), silent = TRUE)
        if (class(qhat) == "try-error") {
          return(class(qhat))
        }
      }
      
      for(risk_alpha in seq(1, risk_max, length.out = length_alpha)){

        if(actual == TRUE){
          errmat = matrix(rnorm(nrow(List2)*3, 1, omega)/nrow(List2)^alpha, ncol = 3)
          q12 = expit(logit(apply(List2[,-c(1:K)], 1, function(x){return(pi1(x)*pi2_1(x))})) + errmat[,1])
          q1 = apply(List2[,-c(1:K)], 1, function(x){return(pi2(x))})
          q2 = apply(List2[,-c(1:K)], 1, function(x){return(pi1(x))})
          q1 = pmax(q12, expit(logit(q1) + errmat[,2]))
          q2 = pmin(pmax(q12/q1, expit(logit(q2) + errmat[,3])), 1 + q12 - q1)
          
          result = risklim(q1 = q1, q2 = q2, q12 = q12, yi = yi, yj = yj, risk_alpha = risk_alpha)
          psiinvmat = rbind(psiinvmat, c(risk_alpha, result$psiinvmat))
          varmat = rbind(varmat, c(risk_alpha, result$varmat))
        }else{
          q12 = qhat$q12
          q1 = pmin(pmax(q12, qhat$q1), 1)
          q2 = pmax(q12/q1, pmin(qhat$q2, 1 + q12 - q1, 1))
          
          result = risklim(q1 = q1, q2 = q2, q12 = q12, yi = yi, yj = yj, risk_alpha = risk_alpha)
          psiinvmat = rbind(psiinvmat, c(risk_alpha, func, result$psiinvmat))
          varmat = rbind(varmat, c(risk_alpha, func, result$varmat))
        }
      }
    }else{
      cat("Overlap between the lists", i, "and", j, "is less than", eps, "\n")
    }
  }
  colnames(psiinvmat)[1] = "risk_alpha"
  colnames(varmat)[1] = "risk_alpha"
  if(actual == FALSE){
    colnames(psiinvmat)[2] = "model"
    colnames(varmat)[2] = "model"
  }
  return(list(psiinvmat = psiinvmat, varmat = varmat))
}

estim_multi_risk = function(List_matrix, K, l, func = "logit", risk_max = 1, i = i, j = j, length_alpha = 10, actual = FALSE, alpha = 0.25, omega = 1, eps = 0.005){
  
  n = nrow(List_matrix)
  
  if(missing(K)){
    K = ncol(List_matrix) - 1
  }
  
  l = ncol(List_matrix) - K
  List_matrix = na.omit(List_matrix)
  #removing all rows with only 0's
  List_matrix_cov = List_matrix[which(rowSums(List_matrix[,1:K]) > 0),]
  colnames(List_matrix_cov) = c(paste("L", 1:K, sep = ''), paste("x", 1:l, sep = ''))
  #N = number of observed or captured units
  N = nrow(List_matrix_cov)
  
  set1 = sample(N, ceiling(N/2), replace = FALSE)
  List1 = as.data.frame(List_matrix_cov[set1,])
  List2 = as.data.frame(List_matrix_cov[-set1,])
  
  pp = phihat_risk(List1, List2, K = K, func = func, i = i, j = j, risk_max = risk_max, length_alpha = length_alpha, actual = actual, alpha = alpha, omega = omega)
  
  if(pp == "try-error"){
    print("did not run")
    return(pp)
  }
  psimat = 1/pp$psiinvmat
  sigma2mat = pp$varmat
  nmat = N/psimat
  sigma2nmat = N*(sigma2mat + (1 - psimat)/psimat^2)
  psimat[,1] = pp$psiinvmat[,1]
  nmat[,1] = pp$psiinvmat[,1]
  sigma2nmat[,1] = pp$psiinvmat[,1]
  colnames(psimat) = colnames(pp$psiinvmat)[c(1,3,2,5,4,7,6)]
  
  return(list(psimat = as.data.frame(psimat), sigma2mat = as.data.frame(sigma2mat), nmat = as.data.frame(nmat), sigma2nmat = as.data.frame(sigma2nmat), N = N))
}

oddslim <- function(q1, q2, q12, yi, yj, odds_alpha){
  
  psiinvmat = rep(NA, 6)
  names(psiinvmat) = paste(rep(c("PI", "DR", "TMLE"), each = 2), c("l", "u"), sep = '.')
  varmat = psiinvmat
  
  for(lim in c('l', 'u')){
    if(lim == 'u'){
      delta = rep(odds_alpha, length(yi))
      delta = pmax(delta, (1-q1-q2+q12)*q12/(q1-q12)*(q2-q12))
    }else{
      delta = rep(1/odds_alpha, length(yi))
      delta = pmax(delta, (1-q1-q2+q12)*q12/(q1-q12)*(q2-q12))
    }
    gammainvhat = delta*q1*q2/q12 + (1 - delta)*(q1 + q2 - q12)
    psiinvhat = mean(gammainvhat, na.rm = TRUE)
    
    phihat = delta*q1*q2/q12*(yi/q1 + yj/q2 - yi*yj/q12) + (1 - delta)*(yi + yj - yi*yj)
    
    Qnphihat = mean(phihat, na.rm = TRUE)
    
    psiinvhatq = max(Qnphihat, 1)
    
    psiinvmat[paste(c('PI', "DR"), lim, sep = '.')] = c(psiinvhat, psiinvhatq)
    varmat[paste(c('PI', "DR"), lim, sep = '.')] = c(var(phihat), var(phihat))
    
    #print(c(psiinvhat, psiinvhatq))
    
    datmat = as.data.frame(cbind(yi, yj, yi*yj, q1 - q12, q2 - q12, q12, yi*(1 - yj), yj*(1 - yi)))
    datmat[,4:6] = cbind(apply(datmat[,4:6], 2, function(u) {return(pmin(pmax(u, eps), 1 - eps))}))
    colnames(datmat) = c("yi", "yj", "yij", "q10", "q02", "q12", "yi0", "y0j")
    
    epsilon_error = 1
    cnt = 0
    iter = 10; eps = 0.005; twolist = FALSE
    while (abs(epsilon_error) > 0.00001){
      cnt = cnt + 1
      datmat_old = datmat
      if (cnt > iter){break}
      
      ########################### model 1 for q12
      dat1 = cbind(datmat$yij, logit(datmat$q12), delta*(datmat$q10 + datmat$q12)/datmat$q12
                   + delta*(datmat$q02 + datmat$q12)/datmat$q12
                   - delta*(datmat$q10 + datmat$q12)*(datmat$q02 + datmat$q12)/datmat$q12^2 + 1 - delta)
      colnames(dat1) = c("yij", "logitq12", "ratio")
      dat1 = as.data.frame(dat1)
      mod1 = try(glm(yij ~ -1 + offset(logitq12) + ratio
                     , family = binomial(link = logit), data = dat1, na.action = na.omit))
      if (class(mod1) != "try-error"){
        datmat[,"q12"] = predict(mod1, newdata = dat1, type = "response")
        
      }
      summary(cbind(datmat[,4:6], datmat_old[,4:6]))
      datmat$q12 = pmax(pmin(datmat$q12, 1), eps)
      
      ########################### model 2 for q1
      dat2 = cbind(datmat$yi0, logit(datmat$q10), delta*(datmat$q02 + datmat$q12)/datmat$q12)
      colnames(dat2) = c("yi0", "logitq10", "ratio")
      dat2 = as.data.frame(dat2)
      mod2 = try(glm(yi0 ~ -1 + offset(logitq10) + ratio, family = binomial(link = logit), data = dat2, na.action = na.omit))
      if (class(mod2) != "try-error"){
        datmat$q10 = predict(mod2, newdata = dat2, type = "response")
        datmat[,"q10"] = pmin(datmat[,"q10"], 1 - datmat$q12)
      }
      summary(cbind(datmat[,4:6], datmat_old[,4:6]))
      datmat$q10 = pmax(pmin(datmat$q10, 1), eps)
      
      ########################### model 3 for q2
      if (K > 2 | twolist == FALSE){
        dat3 = cbind(datmat$y0j, logit(datmat$q02), delta*(datmat$q10 + datmat$q12)/datmat$q12 + 1 - delta)
        colnames(dat3) = c("y0j", "logitq02", "ratio")
        dat3 = as.data.frame(dat3)
        mod3 = try(glm(y0j ~ -1 + offset(logitq02) + ratio, family = binomial(link = logit), data = dat3, na.action = na.omit))
        if (class(mod3) != "try-error"){
          datmat$q02 = predict(mod3, newdata = dat3, type = "response")
          datmat[,"q02"] = pmin(datmat[,"q02"], 1 - datmat$q10 - datmat$q12)
        }
      }else{
        mod3 = mod2
        datmat[,"q02"] = pmax(0, 1 - datmat$q10 - datmat$q12)
      }
      
      datmat$q02 = pmax(pmin(datmat$q02, 1), eps)
      
      epsilon_error = max(abs(c(mod2$coefficients, mod3$coefficients, mod1$coefficients)))
    }
    
    if(epsilon_error < 1){
      
      q12t = pmax(datmat$q12, eps)
      q1t = pmin(datmat$q12 + datmat$q10, 1)
      q2t = pmax(pmin(datmat$q12 + datmat$q02, 1 + q12t - q1t, 1), q12t/q1t)
      
      gammainvhat = delta*q1t*q2t/q12t + (1 - delta)*(q1t + q2t - q12t)
      psiinvhat = mean(gammainvhat, na.rm = TRUE)
      
      phihat = delta*q1t*q2t/q12t*(yi/q1t + yj/q2t - yi*yj/q12t) + (1 - delta)*(yi + yj - yi*yj)
      
      Qnphihat = mean(phihat, na.rm = TRUE)
      
      psiinvhatq = max(Qnphihat, 1)
      
      psiinvmat[paste0("TMLE.", lim)] = psiinvhat
      
      varmat[paste0("TMLE.", lim)] = var(phihat)
    }
  }
  return(result = list(psiinvmat = psiinvmat, varmat = varmat))
}

phihat_odds = function(List1, List2, K = 2, i = 1, j = 2, func = "logit", nfolds = 2, twolist = TRUE, eps = 0.05, iter = 50, odds_alpha = 1, odds_min = NA, odds_max = NA, length_alpha = length_alpha, actual, alpha, omega){
  
  l = ncol(List1) - K
  
  conforminglists = apply(List1[,1:K], 2, function(col){return(setequal(col, c(0,1)))})
  if(sum(conforminglists) < 2){
    print("Data is not in the required format or lists are degenerate.")
    return(NULL)
  }
  
  if(sum(conforminglists) < K){
    cat("Lists ", which(conforminglists == FALSE), " are not in the required format.\n")
  }
  
  if(l == 0){
    #renaming the columns of List_matrix for ease of use
    colnames(List1) = c(paste("L", 1:K, sep = ''))
    colnames(List2 = colnames(List1))
    q1 = mean(List_matrix[,i])
    q2 = mean(List_matrix[,j])
    q12 = mean(List_matrix[,i]*List_matrix[,j])
    
    for(delta in delta_vec){
      
      psiinv_summary[paste(i, ", ", j, sep = ''),] = delta*q1*q2/q12 + (1 - delta)*(q1 + q2 - q12)
      var_summary[paste(i, ", ", j, sep = ''),] = q1*q2*(q1*q2 - q12)*(1 - q12)/q12^3/N
    }
    
  }else{
    #renaming the columns of List_matrix for ease of use
    colnames(List1) = c(paste("L", 1:K, sep = ''), paste("x", 1:(ncol(List1) - K), sep = ''))
    colnames(List2) = colnames(List1)
    
    yi = List2[,paste("L", i, sep = '')]
    yj = List2[,paste("L", j, sep = '')]
    
    psiinvmat = numeric(0)
    varmat = numeric(0)
    
    if(mean(List1[,i]*List1[,j]) > eps) {
      
      #colsubset = stringr::str_subset(colnames(psiinv_summary), func)
      if(actual == FALSE){
        qhat = try(get(paste0("qhat_", func))(List1, List2, K, i, j, eps), silent = TRUE)
        if (class(qhat) == "try-error") {
          return(class(qhat))
        }
      }

      for(odds_alpha in seq(1, odds_max, length.out = length_alpha)){
        
        if(actual == TRUE){
          errmat = matrix(rnorm(nrow(List2)*3, 1, omega)/nrow(List2)^alpha, ncol = 3)
          q12 = expit(logit(apply(List2[,-c(1:K)], 1, function(x){return(pi1(x)*pi2_1(x))})) + errmat[,1])
          q1 = apply(List2[,-c(1:K)], 1, function(x){return(pi2(x))})
          q2 = apply(List2[,-c(1:K)], 1, function(x){return(pi1(x))})
          q1 = pmax(q12, expit(logit(q1) + errmat[,2]))
          q2 = pmin(pmax(q12/q1, expit(logit(q2) + errmat[,3])), 1 + q12 - q1)
          
          result = oddslim(q1 = q1, q2 = q2, q12 = q12, yi = yi, yj = yj, odds_alpha = odds_alpha)
          psiinvmat = rbind(psiinvmat, c(odds_alpha, result$psiinvmat))
          varmat = rbind(varmat, c(odds_alpha, result$varmat))
        }else{
          q12 = qhat$q12
          q1 = pmin(pmax(q12, qhat$q1), 1)
          q2 = pmax(q12/q1, pmin(qhat$q2, 1 + q12 - q1, 1))
          
          result = oddslim(q1 = q1, q2 = q2, q12 = q12, yi = yi, yj = yj, odds_alpha = odds_alpha)
          psiinvmat = rbind(psiinvmat, c(odds_alpha, func, result$psiinvmat))
          varmat = rbind(varmat, c(odds_alpha, func, result$varmat))
        }
      }
    }else{
      cat("Overlap between the lists", i, "and", j, "is less than", eps, "\n")
    }
    colnames(psiinvmat)[1] = "odds_alpha"
    colnames(varmat)[1] = "odds_alpha"
    if(actual == FALSE){
      colnames(psiinvmat)[2] = "model"
      colnames(varmat)[2] = "model"
    }
  return(list(psiinvmat = psiinvmat, varmat = varmat))
  }
}

estim_multi_odds = function(List_matrix, K, l, func = "logit", odds_max = 1, i = i, j = j, length_alpha = 10, actual = TRUE, alpha = 0.25, omega = 1){
  
  n = nrow(List_matrix)
  
  if(missing(K)){
    K = ncol(List_matrix) - 1
  }
  
  l = ncol(List_matrix) - K
  List_matrix = na.omit(List_matrix)
  #removing all rows with only 0's
  List_matrix_cov = List_matrix[which(rowSums(List_matrix[,1:K]) > 0),]
  colnames(List_matrix_cov) = c(paste("L", 1:K, sep = ''), paste("x", 1:l, sep = ''))
  #N = number of observed or captured units
  N = nrow(List_matrix_cov)
  
  set1 = sample(N, ceiling(N/2), replace = FALSE)
  List1 = as.data.frame(List_matrix_cov[set1,])
  List2 = as.data.frame(List_matrix_cov[-set1,])
  
  pp = phihat_odds(List1, List2, K = K, func = func, i = i, j = j, odds_max = odds_max, length_alpha = length_alpha, actual = actual, alpha = alpha, omega = omega)
  
  if(pp == "try-error"){
    print("did not run")
    return(pp)
  }
  psimat = 1/pp$psiinvmat
  sigma2mat = pp$varmat
  nmat = N/psimat
  sigma2nmat = N*(sigma2mat + (1 - psimat)/psimat)
  psimat[,1] = pp$psiinvmat[,1]
  nmat[,1] = pp$psiinvmat[,1]
  sigma2nmat[,1] = pp$psiinvmat[,1]
  colnames(psimat) = colnames(pp$psiinvmat)[c(1,3,2,5,4,7,6)]
  
  return(list(psimat = as.data.frame(psimat), sigma2mat = as.data.frame(sigma2mat), nmat = as.data.frame(nmat), sigma2nmat = as.data.frame(sigma2nmat), N = N))
}

normq = function(cvrg = 0.95, mu = 0, sigma = 1, D) {
  quant_vec = seq(-3, 3, length.out = 100)
  #seq(from = qnorm(alpha)-D, to = qnorm(1 - alpha), length.out = 100)
  coverage_vec = sapply(quant_vec, function(q){pnorm(D + q) - pnorm(-q)})
  if(length(coverage_vec) > 0){
    return(quant_vec[which.min(abs(coverage_vec - cvrg))])
  } else{
    return(0)
  }
}

confintrvl = function(mus, varmus, cvrg = 0.95){
  
  q = min(normq(cvrg = cvrg, D = diff(range(mus))/max(varmus)))
  return(c(max(0, mus[1] - q*sqrt(varmus[1])), mus[2] + q*sqrt(varmus[2])))
}


#sapply(1:(K - 1), function(i1) {sapply((i1 + 1):K, function(j1) {
#  return(
#    (sapply(colnames(qmat3), function(xx) str_contains(xx, c(as.character(i1), as.character(j1)), logic = "and")))
#  )})})
