#### Explore the skewness and symmetricity of parameter estimation for NB model
#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(glmGamPoi)
library(fitdistrplus)
library(dplyr)
# library(Seurat)
library(moments)

## Replicate 100 times
nr_cell = c(10, 25, 50, 100, 500, 1000, 5000)
mu_true = c(0.001,0.01,0.05, seq(0.1,0.5,0.1),seq(1,3,0.5), 5, 10, 30)
rpl = c()
rep = as.numeric(args[1])

res = data.frame(mu = c(0.001,0.01,0.05, seq(0.1,0.5,0.1),seq(1,3,0.5), 5, 10, 30), mu_est = NA, var_est = NA, theta_moe = NA, prop_moe_NA = NA, theta_mle = NA, 
                 sw_moe = NA, sw_mle = NA, sk_moe = NA, sk_mle = NA)

for(k in 1:length(nr_cell)){

res = data.frame(mu = c(0.001,0.01,0.05, seq(0.1,0.5,0.1),seq(1,3,0.5), 5, 10, 30), nr_cell = NA, mu_est = NA, var_est = NA, theta_moe = NA, prop_moe_NA = NA, theta_mle = NA,
                 sw_moe = NA, sw_mle = NA, sk_moe = NA, sk_mle = NA)

for(j in 1:length(mu_true)){

# Simulate a matrix with 1000 individuals x 500 cells per ind on the same gene with Poisson distribution
m = 1000
# n = 5000 - as.numeric(args[1]) * 50 + 50
n = nr_cell[k]
# mu = 0.2
mu = res$mu[j]
size = 1 # inverse of the dispersion
data = matrix(NA, nrow = m, ncol = n)

for(i in 1:m){
  # data[i,] = rpois(n, mu)
  data[i,] = rnbinom(n = n, size = size, mu = mu)
}

fit = glmGamPoi::glm_gp(as.matrix(data), size_factors = FALSE, verbose = T)
# This is the original estimate from CR-MLE. No further adjustment or shrinkage is added.
theta_mle = fit$overdispersions

m_est = rowMeans(fit$Mu, na.rm=T)

v_est = apply(data, 1, function(x)var(x, na.rm=T))

theta_moe = (v_est-m_est)/(m_est^2)

sw_moe = tryCatch(shapiro.test(theta_moe), error = function(e) {a=list();a["p.value"]=NA;return(a)} )
sw_mle = tryCatch(shapiro.test(theta_mle), error = function(e) {a=list();a["p.value"]=NA;return(a)} )

res[j,] = c(res$mu[j], n, mean(m_est,na.rm=T), mean(v_est,na.rm=T), mean(theta_moe,na.rm=T), prop_moe_NA = sum(is.na(theta_moe)) / m, mean(theta_mle,na.rm=T), sw_moe$p.value, sw_mle$p.value, 
            skewness(theta_moe,na.rm=T), skewness(theta_mle,na.rm=T))

# print(j)
}

rpl = rbind(rpl, res)
print(k)
}

rpl$replicate = rep


write.table(rpl,paste0("theta_est_inflation_results_theta1_rep",rep,".txt"), row.names = F, col.names = T, quote = F)




