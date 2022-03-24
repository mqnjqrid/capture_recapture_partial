expit = function(x) {
  exp(x)/(1 + exp(x))
}
logit = function(x) {
  log(x/(1 - x))
}

n0 = 5000

l = 3
risk = 1.2
##   K = 1 case
#alp = -3#0.42
#-2 #0.72 ## parameter to control capture probability
#-1.5 #0.84
#-1 #0.8
#0#0.96
##   K = 3 case
alp = #-5#0.28
  #-4.3
  #-4#0.5
-3.5#0.57
#-2 #0.86 ## parameter to control capture probability
#-1.5 #0.84
#-1 #0.8
#0#0.96
pi1 = function(x) {
  expit (alp +0.6 + sum(c(0.3#
  )*x))
}
if(risk >= 1){
  pi2_0 = function(x) {
  expit (alp + 0.6 + sum(c(0.5#
  )*x))/risk
}
  pi2_1 = function(x) {
  expit (alp + 0.6 + sum(c(0.5#
  )*x))
}
} else {
  pi2_0 = function(x) {
    expit (alp + 0.6 + sum(c(0.5#
    )*x))
  }
  pi2_1 = function(x) {
    expit (alp + 0.6 + sum(c(0.5#
    )*x))*risk
  }
}
pi2 = function(x){
  pi2_0(x)*(1 - pi1(x)) + pi2_1(x)*pi1(x)
}
dat_p_K = function(n, K = 2, l){
  x = matrix(rnorm(n*l, 2, 1), nrow = n, ncol = l)#
#  y1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - pi1(xi), pi1(xi)))}))
#  y2 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - pi2(xi), pi2(xi)))}))
  p1 = unlist(apply(x, 1, function(xi) {pi1(xi)}))
  p2 = unlist(apply(x, 1, function(xi) {pi2(xi)}))
  y1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - pi1(xi), pi1(xi)))}))
  ## list 2 depends on list 1 to violate independence assumption
  y2_0 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - pi2_0(xi), pi2_0(xi)))}))
  y2_1 = unlist(apply(x, 1, function(xi) {sample(c(0, 1), 1, replace = TRUE, prob = c( 1 - pi2_1(xi), pi2_1(xi)))}))
  y2 = y2_0 * (1 - y1) + y2_1 * y1

  xp = do.call("cbind", lapply(1:ncol(x),
                               function(li){
                                 if(li%%4 == 1){
                                   return(exp(x[,li]/2))
                                 }else if(li%%4 == 2){
                                   return(x[,li]/(1 + exp(x[,li -1])) + 10)
                                 }else if(li%%4 == 3){
                                   return((x[,li]*x[,li-2]/25 + 0.6)^3)
                                 }else{
                                   return((x[,li -2] + x[,li] + 20)^2)
                                 }
                               }))
  List_matrix = cbind(y1, y2, x)
  List_matrix_xstar = cbind(y1, y2, xp )# 5*(sin(5*x)))
  return(list(List_matrix = List_matrix, List_matrix_xstar = List_matrix_xstar, p1 = p1, p2 = p2))
}
psi0 = 1 - mean(apply(matrix(rnorm(n0*l, 2, 1), ncol = l), 1, function(x){return((1 - pi1(x))*(1 - pi2_0(x)))}))
psi0

xmat = matrix(rnorm(n0*l, 2, 1), nrow = n0, ncol = l)
riskvec = apply(xmat, 1, function(x){return(pi1(x)*pi2_1(x)/((1 - pi1(x))*pi2_0(x))*(1 - pi1(x))/pi1(x))})
#hist(riskvec)


# List_matrix = dat_p_K(n0, l = l)$List_matrix
# psi = 1 - mean((1-List_matrix[,1])*(1-List_matrix[,2]))
# risk = (mean(List_matrix[List_matrix[,1] == 1,2])/mean(List_matrix[List_matrix[,1] == 0,2]))
# psi;risk
# 
# risk = sum(List_matrix[,1] == 1 & List_matrix[,2] == 1)/sum(List_matrix[,1] == 0 & List_matrix[,2] == 1)*
#   sum(List_matrix[,1] == 0)/sum(List_matrix[,1] == 1)   
# risk
# 
# xmat = matrix(rnorm(n0*l, 2, 1), nrow = n0, ncol = l)
# riskvec = apply(xmat, 1, function(x){return(pi1(x)*pi2(x)/((1 - pi1(x))*pi2_0(x))*(1 - pi2(x))/pi2(x))})