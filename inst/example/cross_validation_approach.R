#this script is for proof of concept and some targeted testing of inferring multiple maternity via cross validation

library(paternityDP)
library(abind)

#for reproducible results
set.seed(29800)

#set realistic allele frequencies here
allele_freqs <- rep(list(rep(1/5,5)),5)

#----------- 3 mothers x 3 fathers, monogamy, have phenotypes for all --------------#
num_fathers <- 1 #per mother
num_mothers <- 2
dropout_rate <- 0.01
mistyping_rate <- 0.01
proportion_missing_data <- 0.0
number_of_offspring <- 6
num_loci <- length(allele_freqs)

data <- list()
colonies <- list()
for (mother in 1:num_mothers)
{
colonies[[mother]] <- simulate_sibling_group(number_of_offspring=number_of_offspring,
                              proportion_of_sperm_per_father=rep(1/num_fathers,num_fathers),
                              allele_frequencies_per_msat=allele_freqs,
                              rate_of_allelic_dropout_per_locus=rep(dropout_rate,num_loci),
                              rate_of_allelic_mistyping_per_locus=rep(mistyping_rate,num_loci),
                              probability_of_missing_data=proportion_missing_data)
data[[mother]] <- array("0", c(2,0,num_loci))
data[[mother]] <- abind(data[[mother]], colonies[[mother]]$observed_maternal_genotypes, along=2)
data[[mother]] <- abind(data[[mother]], abind(colonies[[mother]]$observed_paternal_genotypes, colonies[[mother]]$observed_paternal_genotypes, along=1), along=2)
data[[mother]] <- abind(data[[mother]], colonies[[mother]]$observed_offspring_genotypes, along=2)
dimnames(data[[mother]])[[2]] <- paste0("colony", mother, "_", dimnames(data[[mother]])[[2]])
}

list_of_genotype_arrays_to_txt(data, filename = "example_data.txt")

phenotypes <- genotype_array_from_txt("example_data.txt")

#proof of c
oof <- paternityDP:::sample_parentage_and_error_rates_from_joint_posterior(
  phenotypes, 
  grep("mother", dimnames(phenotypes)[[2]]),
  grep("father", dimnames(phenotypes)[[2]]),
  c(3:4),
  number_of_mcmc_samples = 1000,
  thinning_interval = 1
  )

#check overfitting
phenotypes2 <- paternityDP:::add_unsampled_parents_to_phenotype_array(phenotypes, 2)
oof2 <- paternityDP:::sample_parentage_and_error_rates_from_joint_posterior(
  phenotypes2,                                                                                
  grep("mother", dimnames(phenotypes2)[[2]]),
  grep("father", dimnames(phenotypes2)[[2]]),
  c(3,4),
  number_of_mcmc_samples = 1000,
  thinning_interval = 1
  )

#check underfitting
drop_mom <- grep("mother", dimnames(phenotypes)[[2]])[1]
drop_dad <- grep("father", dimnames(phenotypes)[[2]])[1]
phenotypes3 <- phenotypes[,-c(drop_mom,drop_dad),]
oof3 <- paternityDP:::sample_parentage_and_error_rates_from_joint_posterior(
  phenotypes3,                                                                                
  grep("mother", dimnames(phenotypes3)[[2]]),
  grep("father", dimnames(phenotypes3)[[2]]),
  c(1,2),
  number_of_mcmc_samples = 1000,
  thinning_interval = 1
  )


phenotypes4 <- phenotypes
phenotypes4 <- phenotypes4[,-grep("mother", dimnames(phenotypes4)[[2]])[-1],]
phenotypes4 <- phenotypes4[,-grep("father", dimnames(phenotypes4)[[2]]),]
woo <- cross_validation_for_number_of_parents(phenotypes4, 1, 5, 5)
plot_cross_validation_scores(woo)

