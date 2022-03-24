
if(FALSE){
  List_matrix = dat_p_K(5000, K = 2, l = 3)$List_matrix
  K = 2; l = 3; odds_max = 2; i = 1; j = 2; length_alpha = 3; alpha = 0.25; omega = 1
  actual = TRUE; iter = 10; eps = 0.005; twolist = FALSE
  pp = phihat_odds(List1, List2, K = K, func = func, i = i, j = j, odds_max = odds_max, length_alpha = length_alpha, actual = actual, alpha = alpha, omega = omega)
  
  estim_multi_odds(List_matrix, K, l, odds_max = 2, i = i, j = j, length_alpha = length_alpha, actual = actual, alpha = alpha, omega = omega)
}

library(reshape2)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(glue)
library(colorspace)
library(beepr)
library(tidyr)

n0 = 5000
K = 2
l = 3
eps = 0.005

source("C:/Users/manja/Dropbox/capture_recapture_partial/codes/multi_list/simulation_multi_list_odds.R")
#source("C:/Users/manja/Dropbox/capture_recapture/codes/multi_list/functions_multi_psi_theta_new_notation.R")
source("C:/Users/manja/Dropbox/capture_recapture_partial/codes/multi_list/functions_multilist3.R")

func = "rangerlogit"
length_alpha = 10
odds_max = 2
actual = TRUE

#n_vec= c(1000, 5000, 10000, 15000, 20000)
#n_vec= c(25000, 30000, 35000, 40000, 45000)
n_vec = 5000#c(1:3)*1000
psi0

cvrgprob_vec = c(0.95)

# data table to store the estimated values
pestim = numeric(0)

for(n0 in n_vec){print(n0)
  for(s in 1:100){print(s)
    #print(s)
    ################## generating data and calculating psi, n estimates
    datap = dat_p_K(n0, K, l)
    if(actual == TRUE){
      List_matrix = datap$List_matrix
    }else{
      List_matrix = datap$List_matrix_xstar
    }
    #       List_matrix = List_matrix[colSums(List_matrix[,1:K])>0,]
    estim = estim_multi_odds(List_matrix, K, l, actual = actual, odds_max = odds_max, i = 1, j = 2, alpha = 0.25, omega = 1, length_alpha = length_alpha)
    
    p1 = melt(estim$psimat, value.name = "psi", id.vars = "odds_alpha")
    p2 = melt(estim$sigma2mat, value.name = "sigma2", id.vars = "odds_alpha")
    p3 = melt(estim$nmat, value.name = "n", id.vars = "odds_alpha")
    p4 = melt(estim$sigma2nmat, value.name = "sigma2n", id.vars = "odds_alpha")
    
    pall = merge(p1, merge(p2, merge(p3, p4, by = c("odds_alpha", "variable")), by = c("odds_alpha", "variable")), by = c("odds_alpha", "variable"))
    
    pall$n0 = n0
    pall$N = estim$N
    
    pestim = rbind(pestim, pall)
  }
  beep(sound = 10)
  print(summary(pall))
}
colnames(pestim)[2] = "model"

dtl = dat_p_K(n = 10000, K = 2, l = l)
dt = dtl$List_matrix
p1 = dtl$p1[dt[,1] + dt[,2] > 0]
p2 = dtl$p2[dt[,1] + dt[,2] > 0]
dt = dt[dt[,1] + dt[,2] > 0,]

gamm = 1 - (1 - p1)*(1 - p2)
q1 = p1/gamm
q2 = p2/gamm
q12 = p1*p2/gamm
psi0l = data.frame(odds_alpha = sort(unique(pestim$odds_alpha)), psi0 = 1)
psi0l$psi0 = sapply(psi0l$odds_alpha, 
                    function(delta) {
                      delta1 = pmax(1/delta, (1-q1-q2+q12)*q12/(q1-q12)*(q2-q12))
                      1/mean((delta1*(q1-q12) + q12)*q2/q12)})

psi0u = data.frame(odds_alpha = sort(unique(pestim$odds_alpha)), psi0 = 1)
psi0u$psi0 = sapply(psi0u$odds_alpha, 
                    function(delta) {
                      delta1 = pmax(1/delta, (1-q1-q2+q12)*q12/(q1-q12)*(q2-q12))
                      1/mean((delta1*(q1-q12) + q12)*q2/q12)})

psi0l$psiinv0 = 1/psi0u$psi0
psi0u$psiinv0 = 1/psi0l$psi0
psi0l$bound = 'l'
psi0u$bound = 'u'
psi0mat = rbind(psi0l, psi0u)

pclim = pestim
pclim = separate(pclim, col = "model", into = c("method", "bound"), sep = "\\.")

pclim = merge(pclim, psi0mat, by = c("odds_alpha", "bound"))
pclim$varpsi = pclim$sigma2/pclim$N*pclim$psi^4
pclim$biaspsi = abs(pclim$psi - pclim$psi0)
pclim$biasn = abs(pclim$n - pclim$N*pclim$psiinv0)
pclim$msepsi = ((pclim$psi - pclim$psi0)^2)
pclim$msen = ((pclim$n - pclim$N*pclim$psiinv0)^2)

pclim = merge(pclim %>% filter(bound == 'l'), pclim %>% filter(bound == 'u'), suffixes = c('.l', '.u'), by = c("odds_alpha", 'n0', 'N', "method")) %>% subset(select = -c(bound.l, bound.u))
pclim$psi0 = psi0

coveragemat = do.call("rbind", lapply(1:nrow(pclim), function(i){
  c(confintrvl(mus = c(pclim$psi.l[i], pclim$psi.u[i]), varmus = c(pclim$varpsi.l[i], pclim$varpsi.u[i])),
    confintrvl(mus = c(pclim$n.l[i], pclim$n.u[i]), varmus = c(pclim$sigma2n.l[i], pclim$sigma2n.u[i])))
}))
colnames(coveragemat) = c("psici.l", "psici.u", "nci.l", "nci.u")
pclim = cbind.data.frame(pclim, coveragemat)

#$pclim$psici.l = pclim$psi.l - 1.95*sqrt(pclim$varpsi.l)
#pclim$psici.u = pclim$psi.u + 1.95*sqrt(pclim$varpsi.u)
#pclim$nci.l = pclim$n.l - 1.95*sqrt(pclim$sigma2n.l)
#pclim$nci.u = pclim$n.u + 1.95*sqrt(pclim$sigma2n.u)

pclim$cvrgpsi = (psi0 < pclim$psici.u)*(psi0 > pclim$psici.l)
pclim$cvrgn = (pclim$n0 < pclim$nci.u)*(pclim$n0 > pclim$nci.l)
pclim = aggregate(cbind(psi.l, varpsi.l, n.l, sigma2n.l, psici.l, nci.l, psi.u, biaspsi.l, msepsi.l, biasn.l, msen.l,
                        varpsi.u, n.u, sigma2n.u, psici.u, nci.u, biaspsi.u, msepsi.u, biasn.u, msen.u, cvrgpsi, cvrgn)~odds_alpha + method + n0, data = pclim, mean)
pclim$rmsepsi.l = sqrt(pclim$msepsi.l)
pclim$rmsen.l = sqrt(pclim$msen.l)
pclim$rmsepsi.u = sqrt(pclim$msepsi.u)
pclim$rmsen.u = sqrt(pclim$msen.u)

#load('C:/Users/manja/Dropbox/capture_recapture_partial/codes/multi_list/dataodds_K2_l3_psi044.Rdata')

color1 = "red"
color2 = "#E69F00"
color3 = "dodgerblue1" #  "#56B4E9"
lsize = 0.5
tsize = 12
options(scipen = 5)
gbasic = ggplot(pclim[pclim$method != "PI",], aes(x = odds_alpha, color = method, fill = method)) +
  scale_color_manual(values = c("PI" = color1, "DR" = color2, "TMLE" = color3)) +
  scale_fill_manual(values = c("PI" = color1, "DR" = color2, "TMLE" = color3)) +
  ylab(NULL) +
  facet_wrap(~n0) +
  labs(x = "odds_ratio") +
  #geom_vline(xintercept = odds) +
  theme_bw() +
  theme(text = element_text(size = tsize)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))

g1 = gbasic + geom_line(aes(y = biaspsi.l), size = lsize, linetype = "dotdash") +
  geom_line(aes(y = biaspsi.u), size = lsize) +
  geom_text(aes(x = 1.15, y = 0.15, label = " solid = upper\n broken = lower"), color = "brown") +
  ggtitle(substitute(paste("Bias for ", psi)))
g2 = gbasic + geom_line(aes(y = rmsepsi.l), size = lsize, linetype = "dotdash") +
  geom_line(aes(y = rmsepsi.u), size = lsize) +
  ggtitle(substitute(paste("RMSE for ", psi)))
g3 = gbasic + geom_line(aes(y = psi.l), size = lsize) +
  geom_line(aes(y = psi.u), size = lsize) +
  geom_line(aes(y = psici.l), size = lsize, linetype = "dashed") +
  geom_line(aes(y = psici.u), size = lsize, linetype = "dashed") +
  ggtitle(substitute(paste("CI for ", psi)))
g4 = gbasic + geom_line(aes(y = cvrgpsi), size = lsize) +
  ggtitle(substitute(paste("Coverage for ", psi)))

g5 = gbasic + geom_line(aes(y = biasn.l), size = lsize, linetype = "dotdash") +
  geom_line(aes(y = biasn.u), size = lsize) +
  ggtitle(paste("Bias for n"))
g6 = gbasic + geom_line(aes(y = rmsen.l), size = lsize, linetype = "dotdash") +
  geom_line(aes(y = rmsen.u), size = lsize) +
  ggtitle(paste("RMSE for n"))
g7 = gbasic + geom_line(aes(y = n.l), size = lsize) +
  geom_line(aes(y = n.u), size = lsize) +
  geom_line(aes(y = nci.l), size = lsize, linetype = "dashed") +
  geom_line(aes(y = nci.u), size = lsize, linetype = "dashed") +
  ggtitle(paste("CI for n"))
g8 = gbasic + geom_line(aes(y = cvrgn), size = lsize) +
  ggtitle(paste("Coverage for n"))


ggall = ggarrange(g1, g2, g3, g4, g5, g6, g7, g8, nrow = 2, ncol = 4, common.legend = TRUE, legend = "bottom")

ggall = annotate_figure(ggall, top = text_grob(paste("Summary for capture probability", round(psi0, 2)), size = 15))
ggall
ggall

pdf(paste0("C:/Users/manja/Dropbox/capture_recapture_partial/codes/images/lineplot_odds_K", K, '_l', l, '_psi0', round(psi0*100), '.pdf'), height = 6.5, width = 15, onefile = FALSE)
#ggarrange(g11, g5, nrow = 1, common.legend = TRUE, legend = "bottom")
ggall
dev.off()

save(pestim, psi0mat, psi0, n0, pclim, file = paste0("C:/Users/manja/Dropbox/capture_recapture_partial/codes/multi_list/dataodds_K", K, '_l', l, '_psi0', round(psi0*100), '.Rdata'))

