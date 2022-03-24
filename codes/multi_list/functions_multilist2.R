library(drpop)
phihat_risk = function(List1, List2, K = 2, i = 1, j = 2, func = "logit", nfolds = 2, twolist = TRUE, eps = 0.05, iter = 50, risk_alpha = 1, risk_min = NA, risk_max = NA, length = length, actual, alpha, omega){

  l = ncol(List1) - K

  conforminglists = apply(List1[,1:K], 2, function(col){return(setequal(col, c(0,1)))})
  if(sum(conforminglists) < 2){
    print("Data is not in the required format or lists are degenerate.")
    return(NULL)
  }

  if(sum(conforminglists) < K){
    cat("Lists ", which(conforminglists == FALSE), " are not in the required format.\n")
  }

  risk_alpha_min = max(min(risk_alpha, 1/risk_alpha), risk_min, na.rm = TRUE)
  risk_alpha_max = min(max(risk_alpha, 1/risk_alpha), risk_max, na.rm = TRUE)
  delta_vec = unique(seq(risk_alpha_min, risk_alpha_max, length.out = length))

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

    psiinv_pi = matrix(0, ncol = length(delta_vec), nrow = 1)
    colnames(psiinv_pi) = delta_vec
    psiinv_bc = psiinv_pi
    psiinv_tmle = psiinv_pi

    phi_bc = matrix(0, ncol = length(delta_vec), nrow = nrow(List2))
    colnames(phi_bc) = colnames(psiinv_pi)
    phi_tmle = phi_bc
    
    sigbc = matrix(0, ncol = length(delta_vec), nrow = 1)
    colnames(sigbc) = colnames(psiinv_pi)

    yi = List2[,paste("L", i, sep = '')]
    yj = List2[,paste("L", j, sep = '')]

    if(mean(List1[,i]*List1[,j]) > eps) {

      #colsubset = stringr::str_subset(colnames(psiinv_summary), func)
      if(actual == FALSE){
        qhat = try(get(paste0("qhat_", func))(List1, List2, K, i, j, eps), silent = TRUE)
        if (class(qhat) == "try-error") {
          return(class(qhat))
        }
      }
      
      for(coli in 1:length(delta_vec)){
        delta = delta_vec[coli]
        if(actual == TRUE){
          errmat = matrix(rnorm(nrow(List2)*3, 1, omega)/nrow(List2)^alpha, ncol = 3)
          q12 = expit(logit(apply(List2[,-c(1:K)], 1, function(x){return(pi3(x)*pi4_1(x))})) + errmat[,1])
          q1 = apply(List2[,-c(1:K)], 1, function(x){return(pi4(x))})
          q2 = apply(List2[,-c(1:K)], 1, function(x){return(pi3(x))})
          q1 = pmax(q12, expit(logit(q1) + errmat[,2]))
          q2 = pmin(pmax(q12/q1, expit(logit(q2) + errmat[,3])), 1 + q12 - q1)
        }else{
          q12 = qhat$q12
          q1 = pmin(pmax(q12, qhat$q1), 1)
          q2 = pmax(q12/q1, pmin(qhat$q2, 1 + q12 - q1, 1))
        }
        gammainvhat = (delta*(q1 - q12) + q12)*q2/q12
        psiinvhat = mean(gammainvhat, na.rm = TRUE)

        phihat = gammainvhat*(yj/q2 + (delta*(yi - yi*yj) + yi*yj)/(delta*(q1 - q12) + q12) - yi*yj/q12)

        Qnphihat = mean(phihat, na.rm = TRUE)

        psiinvhatq = max(Qnphihat, 1)

        psiinv_pi[,coli] = psiinvhat
        psiinv_bc[,coli] = psiinvhatq
#print(c(psiinvhat, psiinvhatq))
        phi_bc[,coli] = phihat
        
        sigbc[1,coli] = mean(gammainvhat*(gammainvhat - 1)*(1/q12 - 1) + gammainvhat*(1-q1-q2+q12)/q12) + var(gammainvhat) +
          (1-delta)*mean((delta*(q1-q12-q2)+q12)*(q1-q12)*q2/q12^2)
        
        datmat = as.data.frame(cbind(yi, yj, yi*yj, q1 - q12, q2 - q12, q12, yi*(1 - yj), yj*(1 - yi)))
        datmat[,4:6] = cbind(apply(datmat[,4:6], 2, function(u) {return(pmin(pmax(u, eps), 1 - eps))}))
        colnames(datmat) = c("yi", "yj", "yij", "q10", "q02", "q12", "yi0", "y0j")

        epsilon_error = 1
        cnt = 0

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

        if(epsilon_error > 1){
          psiinvmat[folds,colsubset][3] = NA
          varmat[folds,colsubset][3] = NA
        }else{
          q12 = pmax(datmat$q12, eps)
          q1 = pmin(datmat$q12 + datmat$q10, 1)
          q2 = pmax(pmin(datmat$q12 + datmat$q02, 1 + q12 - q1, 1), q12/q1)

          gammainvhat = (delta*(q1 - q12) + q12)*q2/q12
          psiinvhat = mean(gammainvhat, na.rm = TRUE)

          phihat = gammainvhat*(yj/q2 + (delta*(yi - yi*yj) + yi*yj)/(delta*(q1 - q12) + q12) - yi*yj/q12)

          Qnphihat = mean(phihat, na.rm = TRUE)

          psiinv_tmle[,coli] = psiinvhat

          phi_tmle[,coli] = phihat
        }
      }
      }else{
        cat("Overlap between the lists", i, "and", j, "is less than", eps, "\n")
    }
  }

  return(list(#psiinvmat = psiinvmat, varmat = varmat,
    psiinv_pi = psiinv_pi, psiinv_bc = psiinv_bc, psiinv_tmle = psiinv_tmle,
    phi_pi = phi_bc, phi_bc = phi_bc, phi_tmle = phi_tmle, sigbc = sigbc
  ))
}

estim_multi_risk = function(List_matrix, K, l, func = "logit", risk_alpha = 1, i = i, j = j, length = 10, actual = FALSE, alpha = 0.25, omega = 0.5, risk_min = NA, risk_max = NA){

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

  pp = phihat_risk(List1, List2, K = K, func = func, i = i, j = j, risk_alpha = risk_alpha, risk_min = risk_min, risk_max = risk_max, length = length, actual = actual, alpha = alpha, omega = omega)

  if(pp == "try-error"){
    print("did not run")
    return(pp)
  }
  psi_pi = 1/pp$psiinv_pi
  psi_bc = 1/pp$psiinv_bc
  psi_tmle = 1/pp$psiinv_tmle

  sigma2_pi = apply(pp$phi_pi, 2, var)
  sigma2_bc = apply(pp$phi_bc, 2, var)
  sigma2_bcalt = pp$sigbc
  #plot(as.numeric(names(sigma2_bc)), sigma2_bc); points(as.numeric(names(sigma2_bc)), sigma2_bcalt, col = "red")#;abline(0,1)
  sigma2_tmle = apply(pp$phi_tmle, 2, var)

  psimat = cbind(t(psi_pi), t(psi_bc), t(psi_tmle))
  rownames(psimat) = colnames(psi_pi)
  colnames(psimat) = c("PI", "BC", "TMLE")

  sigma2mat = cbind(sigma2_pi, sigma2_bc, sigma2_tmle)
  rownames(sigma2mat) = colnames(psi_pi)
  colnames(sigma2mat) = c("PI", "BC", "TMLE")

  nmat = N/psimat
  sigma2n = N*(sigma2mat + (1 - psimat)/psimat)

  return(list(psimat = psimat, sigma2mat = sigma2mat, nmat = nmat, sigma2n = sigma2n,
              N = N))
}

phihat_odds = function(List1, List2, K = 2, i = 1, j = 2, func = "logit", nfolds = 2, twolist = TRUE, eps = 0.05, iter = 50, odds_alpha = 1, odds_min = NA, odds_max = NA, length = length, actual, alpha, omega){
  
  l = ncol(List1) - K
  
  conforminglists = apply(List1[,1:K], 2, function(col){return(setequal(col, c(0,1)))})
  if(sum(conforminglists) < 2){
    print("Data is not in the required format or lists are degenerate.")
    return(NULL)
  }
  
  if(sum(conforminglists) < K){
    cat("Lists ", which(conforminglists == FALSE), " are not in the required format.\n")
  }
  
  odds_alpha_min = max(min(odds_alpha, 1/odds_alpha), odds_min, na.rm = TRUE)
  odds_alpha_max = min(max(odds_alpha, 1/odds_alpha), odds_max, na.rm = TRUE)
  delta_vec = unique(seq(odds_alpha_min, odds_alpha_max, length.out = length))
  
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
    
    psiinv_pi = matrix(0, ncol = length(delta_vec), nrow = 1)
    colnames(psiinv_pi) = delta_vec
    psiinv_bc = psiinv_pi
    psiinv_tmle = psiinv_pi
    
    phi_bc = matrix(0, ncol = length(delta_vec), nrow = nrow(List2))
    colnames(phi_bc) = colnames(psiinv_pi)
    phi_tmle = phi_bc
    
    sigbc = matrix(0, ncol = length(delta_vec), nrow = 1)
    colnames(phi_bc) = colnames(psiinv_pi)
    
    yi = List2[,paste("L", i, sep = '')]
    yj = List2[,paste("L", j, sep = '')]
    
    if(mean(List1[,i]*List1[,j]) > eps) {
      
      #colsubset = stringr::str_subset(colnames(psiinv_summary), func)
      if(actual == FALSE){
        qhat = try(get(paste0("qhat_", func))(List1, List2, K, i, j, eps), silent = TRUE)
        if (class(qhat) == "try-error") {
          return(class(qhat))
        }
      }
      
      for(coli in 1:length(delta_vec)){
        delta = delta_vec[coli]
        if(actual == TRUE){
          errmat = matrix(rnorm(nrow(List2)*3, 1, omega)/nrow(List2)^alpha, ncol = 3)
          q12 = expit(logit(apply(List2[,-c(1:K)], 1, function(x){return(pi3(x)*pi4_1(x))})) + errmat[,1])
          q1 = apply(List2[,-c(1:K)], 1, function(x){return(pi4(x))})
          q2 = apply(List2[,-c(1:K)], 1, function(x){return(pi3(x))})
          q1 = pmax(q12, expit(logit(q1) + errmat[,2]))
          q2 = pmin(pmax(q12/q1, expit(logit(q2) + errmat[,3])), 1 + q12 - q1)
        }else{
          q12 = qhat$q12
          q1 = pmin(pmax(q12, qhat$q1), 1)
          q2 = pmax(q12/q1, pmin(qhat$q2, 1 + q12 - q1, 1))
        }
        gammainvhat = delta*q1*q2/q12 + (1 - delta)*(q1 + q2 - q12)
        psiinvhat = mean(gammainvhat, na.rm = TRUE)
        
        phihat = delta*q1*q2/q12*(yi/q1 + yj/q2 - yi*yj/q12) + (1 - delta)*(yi + yj - yi*yj)
        
        Qnphihat = mean(phihat, na.rm = TRUE)
        
        psiinvhatq = max(Qnphihat, 1)
        
        psiinv_pi[,coli] = psiinvhat
        psiinv_bc[,coli] = psiinvhatq
        #print(c(psiinvhat, psiinvhatq))
        phi_bc[,coli] = phihat
        
        sigbc[1,coli] = mean(gammainvhat*(gammainvhat - 1)*(1/q12 - 1) + gammainvhat*(1-q1-q2+q12)/q12) + var(gammainvhat) -
          delta*(1-delta)*mean((q1-q12 +q2-q12)*(q1-q12)*(q2-q12)/q12^2)
        
        datmat = as.data.frame(cbind(yi, yj, yi*yj, q1 - q12, q2 - q12, q12, yi*(1 - yj), yj*(1 - yi)))
        datmat[,4:6] = cbind(apply(datmat[,4:6], 2, function(u) {return(pmin(pmax(u, eps), 1 - eps))}))
        colnames(datmat) = c("yi", "yj", "yij", "q10", "q02", "q12", "yi0", "y0j")
        
        epsilon_error = 1
        cnt = 0
        
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
          dat2 = cbind(datmat$yi0, logit(datmat$q10), delta*(datmat$q02 + datmat$q12)/datmat$q12 + 1 - delta)
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
        
        if(epsilon_error > 1){
          psiinvmat[folds,colsubset][3] = NA
          varmat[folds,colsubset][3] = NA
        }else{
          q12 = pmax(datmat$q12, eps)
          q1 = pmin(datmat$q12 + datmat$q10, 1)
          q2 = pmax(pmin(datmat$q12 + datmat$q02, 1 + q12 - q1, 1), q12/q1)
          
          gammainvhat = delta*q1*q2/q12 + (1 - delta)*(q1 + q2 - q12)
          psiinvhat = mean(gammainvhat, na.rm = TRUE)
          
          phihat = delta*q1*q2/q12*(yi/q1 + yj/q2 - yi*yj/q12) + (1 - delta)*(yi + yj - yi*yj)
                                    
          Qnphihat = mean(phihat, na.rm = TRUE)
          
          psiinv_tmle[,coli] = psiinvhat
          
          phi_tmle[,coli] = phihat
        }
      }
    }else{
      cat("Overlap between the lists", i, "and", j, "is less than", eps, "\n")
    }
  }
  
  return(list(#psiinvmat = psiinvmat, varmat = varmat,
    psiinv_pi = psiinv_pi, psiinv_bc = psiinv_bc, psiinv_tmle = psiinv_tmle,
    phi_pi = phi_bc, phi_bc = phi_bc, phi_tmle = phi_tmle, sigbc = sigbc
  ))
}

estim_multi_odds = function(List_matrix, K, l, func = "logit", odds_alpha = 1, i = i, j = j, length = 10, actual = FALSE, alpha = 0.25, omega = 0.5, odds_min = NA, odds_max = NA){
  
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
  
  pp = phihat_odds(List1, List2, K = K, func = func, i = i, j = j, odds_alpha = odds_alpha, odds_min = odds_min, odds_max = odds_max, length = length, actual = actual, alpha = alpha, omega = omega)
  
  if(pp == "try-error"){
    print("did not run")
    return(pp)
  }
  psi_pi = 1/pp$psiinv_pi
  psi_bc = 1/pp$psiinv_bc
  psi_tmle = 1/pp$psiinv_tmle
  
  sigma2_pi = apply(pp$phi_pi, 2, var)
  sigma2_bc = apply(pp$phi_bc, 2, var)
  sigma2_bcalt = pp$sigbc
  #plot(as.numeric(names(sigma2_bc)), sigma2_bc); points(as.numeric(names(sigma2_bc)), sigma2_bcalt, col = "red")#;abline(0,1)
  sigma2_tmle = apply(pp$phi_tmle, 2, var)
  
  psimat = cbind(t(psi_pi), t(psi_bc), t(psi_tmle))
  rownames(psimat) = colnames(psi_pi)
  colnames(psimat) = c("PI", "BC", "TMLE")
  
  sigma2mat = cbind(sigma2_pi, sigma2_bc, sigma2_tmle)
  rownames(sigma2mat) = colnames(psi_pi)
  colnames(sigma2mat) = c("PI", "BC", "TMLE")
  
  nmat = N/psimat
  sigma2n = N*(sigma2mat + (1 - psimat)/psimat)
  
  return(list(psimat = psimat, sigma2mat = sigma2mat, nmat = nmat, sigma2n = sigma2n,
              N = N))
}

normq = function(alpha, mu = 0, sigma = 1, D) {
  quant_vec = seq(-3, 3, length.out = 100)
  #seq(from = qnorm(alpha)-D, to = qnorm(1 - alpha), length.out = 100)
  coverage_vec = sapply(quant_vec, function(q){pnorm(D + q) - pnorm(-q)})
  if(length(coverage_vec) > 0){
    return(quant_vec[which.min(abs(coverage_vec - alpha))])
  } else{
    return(0)
  }
}

confintvrl = function(thetas, covtheta, alpha, N){
  
  if(class(thetas) == "numeric") {
    q = min(normq(alpha = alpha, D = diff(range(thetas))/max(diag(covtheta))))
    return(c(max(0, thetas[1] - q*sqrt(covtheta[1,1])), thetas[2] + q*sqrt(covtheta[2,2])))
  }else{
    ci_matrix = numeric(0)
    for(l in 1:length(covtheta)){
      q = min(normq(alpha = alpha, D = diff(range(thetas[[l]]))/sqrt(max(diag(covtheta[[l]])))))
      ci_matrix = rbind(ci_matrix, c(max(0,thetas[[l]][1] - q*sqrt(covtheta[[l]][1,1])), thetas[[l]][2] + q*sqrt(covtheta[[l]][2,2])))
    }
    0
    return(ci_matrix)
  }
}


#sapply(1:(K - 1), function(i1) {sapply((i1 + 1):K, function(j1) {
#  return(
#    (sapply(colnames(qmat3), function(xx) str_contains(xx, c(as.character(i1), as.character(j1)), logic = "and")))
#  )})})
