library(reshape2)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(glue)
library(colorspace)
library(beepr)

n0 = 5000
K = 2
l = 3

source("C:/Users/manja/Dropbox/capture_recapture_partial/codes/multi_list/simulation_multi_list_odds.R")
#source("C:/Users/manja/Dropbox/capture_recapture/codes/multi_list/functions_multi_psi_theta_new_notation.R")
source("C:/Users/manja/Dropbox/capture_recapture_partial/codes/multi_list/functions_multilist2.R")

func = "rangerlogit"
length = 20
odds_alpha = odds + 0.25
odds_min = odds - 0.25
odds_max = odds + 0.25
actual = FALSE

#n_vec= c(1000, 5000, 10000, 15000, 20000)
#n_vec= c(25000, 30000, 35000, 40000, 45000)
n_vec = 10000#c(1:3)*1000
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
    estim = estim_multi_odds(List_matrix, K, l, actual = actual, func = func, odds_alpha = odds_alpha, i = 2, j = 1, odds_min = odds_min, odds_max = odds_max, alpha = 0.25, omega = 1, length = length)
    
    p1 = melt(estim$psimat, value.name = "psi")
    p2 = melt(estim$sigma2mat, value.name = "sigma2")
    p3 = melt(estim$nmat, value.name = "n")
    p4 = melt(estim$sigma2n, value.name = "sigma2n")
    
    pall = merge(p1, merge(p2, merge(p3, p4, by = c("Var1", "Var2")), by = c("Var1", "Var2")), by = c("Var1", "Var2"))
    colnames(pall)[1:2] = c("delta", "model")
    
    pall$n0 = n0
    pall$N = estim$N
    
    pestim = rbind(pestim, pall)
  }
  beep(sound = 10)
  print(summary(pall))
}

#save(pestim, file = paste0("C://Users/manja/Dropbox/capture_recapture/codes/multi_list/multilist_psi0", round(psi0*100), "_K", K, "_l", l, "_", func, ".Rdata"))

if(FALSE){
  beep(sound = 8)
  psi = matrix(NA, nrow = length(beta_vec), ncol = 3)
  psi[,1] = beta_vec
  colnames(psi) = c("beta", "psi0_l", "psi0_u")
  xmat = matrix(rnorm(n0*l, 2, 1), ncol = l)
  psi_odds = function(delta){
    pmax(pmin(1/mean(apply(xmat, 1, function(x) {
      gammax = 1 - (1 - pi1(x))*(1 - pi2(x))*(1 - pi3(x))*(1 - pi4_0(x))
      q1 = pi3(x)
      q2 = pi3(x)*pi4_1(x) + (1 - pi3(x))*pi4_0(x)
      q12 = pi3(x)*pi4_1(x)
      return((delta*(q1-q12) + q12)*q2/q12)
    })), 1), 0)
  }
  psi[,"psi0_u"] = psi_odds(1/odds_alpha)
  psi[,"psi0_l"] = psi_odds(odds_alpha)
  
  pestim = merge(pestim, psi, by = "beta", all.x = TRUE)
}
#K = 2; l = 3; func = "logit"; psi0 = 0.35
#load(paste0("C://Users/manja/Dropbox/capture_recapture/codes/multi_list/multilist_psi0", round(psi0*100), "_K", K, "_l", l, "_", func, ".Rdata"))

pc = pestim
pc$psi0 = psi0
pc$biaspsi = abs(pc$psi - pc$psi0)
pc$biasn = abs(pc$n - pc$N/pc$psi0)
pc$msepsi = ((pc$psi - pc$psi0)^2)
pc$msen = ((pc$n - pc$N/pc$psi0)^2)
pc$varpsi = pc$sigma2/pc$N*pc$psi^4
pc$cvrgpsi = abs(pc$biaspsi) < 1.95*sqrt(pc$varpsi)
pc$cvrgn = abs(pc$biasn) < 1.95*sqrt(pc$sigma2n)
pc$cvrgpsi = (psi0 < pc$psi + 1.95*sqrt(pc$varpsi))*(psi0 > pc$psi - 1.95*sqrt(pc$varpsi))
pc$cvrgn = (pc$n0 < pc$n + 1.95*sqrt(pc$sigma2n))*(pc$n0 > pc$n - 1.95*sqrt(pc$sigma2n))
pc = aggregate(cbind(psi, varpsi, n, sigma2n, biaspsi, biasn, msepsi, msen, cvrgpsi, cvrgn, psi0)~delta + model + n0, data = pc, mean)
pc$rmsepsi = sqrt(pc$msepsi)
pc$rmsen = sqrt(pc$msen)

color1 = "red"
color2 = "#E69F00"
color3 = "dodgerblue1" #  "#56B4E9"
lsize = 0.5
tsize = 12
options(scipen = 5)
gbasic = ggplot(pc, aes(x = delta, color = model, fill = model)) +
  scale_color_manual(values = c("PI" = color1, "BC" = color2, "TMLE" = color3)) +
  scale_fill_manual(values = c("PI" = color1, "BC" = color2, "TMLE" = color3)) +
  ylab(NULL) +
  facet_wrap(~n0) +
  labs(x = "odds_ratio") +
  geom_vline(xintercept = odds) +
  theme_bw() +
  theme(text = element_text(size = tsize)) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))

g1 = gbasic + geom_line(aes(y = biasn), size = lsize) +
  ggtitle("Bias of n")

g2 = gbasic + geom_line(aes(y = cvrgn), size = lsize) +
  ggtitle("Coverage of n")

g3 = gbasic + geom_line(aes(y = rmsepsi), size = lsize) +
  ggtitle("RMSE of n")

g4 = gbasic +
  geom_line(aes(y = n + 1.95*sqrt(sigma2n)), size = lsize) +
  geom_line(aes(y = n - 1.95*sqrt(sigma2n)), size = lsize) +
  geom_line(aes(y = n0), linetype = "dashed", color = "black") +
  ggtitle('95% CI of n')

g5 = gbasic + geom_line(aes(y = biaspsi), size = lsize) +
  ggtitle("Bias of psi")

g6 = gbasic + geom_line(aes(y = cvrgpsi), size = lsize) +
  ggtitle("Coverage of psi")

g7 = gbasic + geom_line(aes(y = rmsen), size = lsize) +
  ggtitle("RMSE of psi")

g8 = gbasic +
  #  geom_ribbon(data = .%>% filter(variable == "max"), aes(ymin = n0, ymax = n + nq*nsd), alpha = 0.3, show.legend = FALSE) +
  #  geom_ribbon(data = .%>% filter(variable == "min"), aes(ymin = n0, ymax = n - nq*nsd), alpha = 0.3, show.legend = FALSE) +
  geom_line(aes(y = psi + 1.95*sqrt(varpsi)), size = lsize) +
  geom_line(aes(y = psi - 1.95*sqrt(varpsi)), size = lsize) +
  geom_line(aes(y = psi0), linetype = "dashed", color = "black") +
  #ggtitle(glue("95% CI of psi"))
  labs(title = substitute(paste('95% CI of ', psi, ', true ', psi, ' = ', var), list(var = round(psi0, 2))))

g9 = ggplot(pc, aes(x = model, fill = model)) +
  scale_fill_manual(values = c("PI" = color1, "BC" = color2, "TMLE" = color3)) +
  geom_bar(aes(y = cvrgpsi), position = "dodge", stat = "summary", fun = "median") +
  ggtitle("Coverage of psi")
ggall = ggarrange(g1, g2, g3, g4, g5, g6, g7, g8, nrow = 2, ncol = 4, common.legend = TRUE, legend = "bottom")

ggall = annotate_figure(ggall, top = text_grob(paste("Summary for capture probability", round(psi0, 2)), size = 15))
ggall
ggall

g9

pdf(paste0("C:/Users/manja/Dropbox/capture_recapture_partial/codes/images/summary_odds_K", K, '_l', l, '_psi0', round(psi0*100), '.pdf'), height = 6.5, width = 15, onefile = FALSE)
#ggarrange(g11, g5, nrow = 1, common.legend = TRUE, legend = "bottom")
ggall
dev.off()

save(pestim, pc, file = paste0("C:/Users/manja/Dropbox/capture_recapture_partial/codes/multi_list/dataodds_K", K, '_l', l, '_psi0', round(psi0*100), '.Rdata'))

