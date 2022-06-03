#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(sydneyPaternity)

mat_gno <- read.csv("MaternalGenotypes.csv")[1:23,2:12] #check indices
mat_gno <- lapply(apply(mat_gno, 2, strsplit, split="/"), function(x) { x <- table(unlist(x)); x <- x[names(x)!="0"]; x <- x/sum(x)} )
mat_gno <- mat_gno[lapply(mat_gno, length) > 1]

SEED=as.numeric(args[1])
ERR=as.numeric(args[2])
set.seed(SEED)
DF <- data.frame()

for(k in 0:5){

##dummy allele frequencies
#n_all <- 5
#n_loc <- 4
#all_freq <- lapply(1:n_loc, function(x) rep(1/n_all, n_all))

n_loc <- length(mat_gno)
all_freq <- mat_gno
dropout <- rep(ERR, n_loc)
mistype <- rep(ERR, n_loc)

sim <- simulate_sibling_group(c(20-k, rep(1,k)), all_freq, dropout, mistype, 0.0)
pheno <- abind::abind(sim$observed_maternal_genotypes, sim$observed_offspring_genotypes, along=2)
pheno <- array(as.numeric(pheno), dim(pheno))

fit <- sample_paternity_and_error_rates_from_joint_posterior(pheno, mother=1, update_allele_frequencies=FALSE, number_of_mcmc_samples=1100)

burnin <- 101
n_father <- apply(fit$paternity[,burnin:ncol(fit$paternity)], 2, function(x) length(unique(x)))
post_prob <- sapply(1:20, function(i) sum(n_father==i))/length(n_father)
mean_dropout <- mean(fit$dropout_rate[1,burnin:ncol(fit$dropout_rate)])
mean_mistyping <- mean(fit$mistyping_rate[1,burnin:ncol(fit$mistyping_rate)])

df <- data.frame(seed=SEED, mistyping=ERR, dropout=ERR, truth=k+1, n_father=1:20, post_prob=post_prob, post_dropout=mean_dropout, post_mistyping=mean_mistyping)
DF <- rbind(DF, df)
}

save(DF, file=paste0("power_", SEED, ".", ERR*100, ".RData"))

