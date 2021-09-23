library(sydneyPaternity)

#for reproducible results
set.seed(1)

#set realistic allele frequencies here
allele_freqs <- list(
                     "btms0066"=c("133"=0.203, "136"=0.318, "139"=0.066, "142"=0.096,"145"=0.051,"148"=0.267),
                     "b124"=c("247"=0.103, "251"=0.045, "253"=0.115,"255"=0.115,"259"=0.245, "261"=0.061, "263"=0.058, "265"=0.097, "271"=0.030, "275"=0.021, "277"=0.064, "279"=0.009,"287"=0.036),
                     "btern01"= c("128"=0.004, "130"=0.458, "132"=0.152,"134"=0.097,"136"=0.099, "138"=0.103, "140"=0.046, "142"=0.025, "144"=0.017),
                     "bt28" = c("178"=1.00),
                     "btms0062" = c("261"=0.022, "263"=0.022, "267"=0.012,"269"=0.100,"271"=0.028, "275"=0.062, "277"=0.068, "279"=0.046, "281"=0.006, "283"=0.167, "285"=0.056, "287"=0.068,"289"=0.149,"291"=0.100,"295"=0.30,"299"=0.018,"303"=0.040,"305"=0.004),
                     "btms0073" = c("120"=0.053,"123"=0.907,"126"=0.036,"129"=0.004),
                     "bt10"=c("147"=0.009,"149"=0.084,"151"=0.142,"153"=0.202,"155"=0.142,"157"=0.090,"159"=0.145,"161"=0.151,"165"=0.012,"167"=0.024),
                     "bl11"=c("130"=0.132,"132"=0.355,"134"=0.026,"136"=0.263,"138"=0.132,"140"=0.066,"142"=0.026),
                     "bt30"=c("187"=0.100,"190"=0.900),
                     "b96"=c("238"=0.464,"240"=0.332,"244"=0.045,"246"=0.043,"250"=0.115),
                     "btms0081"=c("307"=0.487,"310"=0.513))

# set simulation parameters here, all combinations of values will be considered
grid_of_simulation_parameters <- 
  expand.grid(replicate=1:3, #probably want more like 1:100 for actual analysis
              num_fathers=1:10, #range of fathers in colony
              error_rates_in_simulation=0.01, #used for simulating data
              error_rates_in_estimation=0.01, #used for estimating from simulated data
              proportion_missing_data=0.1, #use something realistic
              number_of_offspring=20, 
              use_true_allele_freqs=FALSE) #if FALSE use uniform allele freqs for estimation

#run simulations
simulation_results <- data.frame()
for(i in 1:nrow(grid_of_true_values)){

num_fathers <- grid_of_simulation_parameters$num_fathers[i]
replicate <- grid_of_simulation_parameters$replicate[i]
error_rates_in_simulation <- grid_of_simulation_parameters$error_rates_in_simulation[i]
error_rates_in_estimation <- grid_of_simulation_parameters$error_rates_in_estimation[i]
proportion_missing_data <- grid_of_simulation_parameters$proportion_missing_data[i]
number_of_offspring <- grid_of_simulation_parameters$number_of_offspring[i]
use_true_allele_freqs <- grid_of_simulation_parameters$use_true_allele_freqs[i]

# 0. Simulate data
num_loci <- length(allele_freqs)
out <- simulate_sibling_group(number_of_offspring=number_of_offspring,
                              proportion_of_sperm_per_father=rep(1/num_fathers,num_fathers),
                              allele_frequencies_per_msat=allele_freqs,
                              rate_of_allelic_dropout_per_locus=rep(error_rates_in_simulation,num_loci),
                              rate_of_allelic_mistyping_per_locus=rep(error_rates_in_simulation,num_loci),
                              probability_of_missing_data=proportion_missing_data)

# 1. Here we output simulation to a file. Not necessary b/c simulation is already an R object, but this
#    allows us to write the rest of the example like you are working with your actual data.
#    The first row in the data is the maternal genotype, this is important later on.
list_of_genotype_arrays_to_txt(list(out$observed_maternal_genotype, out$observed_offspring_genotypes), 
                               filename = "example_data.txt")

# 2. Load the genotype data we just saved into an array
my_genotype_data <- genotype_array_from_txt("example_data.txt")

# 3. Recode alleles to integers (uses positions in "allele_freqs" as integer codes)
my_genotype_data_recoded <- recode_genotype_array_to_integers(my_genotype_data, allele_freqs)

# 4. Pull out maternal genotype (see note above: if not first row in txt, you'll need to modify. Will want to split up colonies too. Ask me later)
mother_index <- 1
maternal_genotype <- my_genotype_data_recoded[,mother_index,]
offspring_genotypes <- my_genotype_data_recoded[,-mother_index,]

# 5. Find MLE of paternity vector
starting_paternity <- rep(1, dim(offspring_genotypes)[2]) #starting state for optimization, everybody has same father
starting_allele_freqs <- if(use_true_allele_freqs) allele_freqs else lapply(allele_freqs, function(x) rep(1/length(x), length(x)))
dropout_rates <- rep(error_rates_in_estimation, length(allele_freqs)) #for now just set to true values (used in simulations)
mistyping_rates <- rep(error_rates_in_estimation, length(allele_freqs)) #for now just set to true values (used in simulations)
estimated_paternity <- optimize_paternity_given_error_rates(starting_paternity,offspring_genotypes,maternal_genotype,starting_allele_freqs,dropout_rates,mistyping_rates)

# 6. Compare estimated and true paternities
mismatch_proportion <- sum(paternity_vector_to_adjacency_matrix(estimated_paternity$paternity) != paternity_vector_to_adjacency_matrix(out$offspring_paternity))/length(paternity_vector_to_adjacency_matrix(estimated_paternity$paternity))
diff_in_fathers <- length(unique(estimated_paternity$paternity)) - length(unique(out$offspring_paternity))

simulation_results <- rbind(simulation_results,
data.frame(num_fathers=num_fathers, replicate=replicate, "mismatch_prop"=mismatch_proportion, "diff_in_fathers"=diff_in_fathers))

}

library(ggplot2)
#plot results, error vs number of fathers
ggplot(simulation_results) + geom_boxplot(aes(x=num_fathers, y=diff_in_fathers, group=num_fathers)) +
  theme_bw() + xlab("True number of fathers") + ylab("Error in number of fathers (estimated - true)")

###----don't include in power analysis
# 7. Get measure of uncertainty (not necessary for power analysis)
# this is super important for your actual data
set.seed(1)
mcmc_paternity <- 
  sample_paternity_given_error_rates(starting_paternity,
                                     offspring_genotypes,
                                     maternal_genotype,
                                     uniform_allele_freqs,
                                     dropout_rates,
                                     mistyping_rates,
                                     max_iter = 100)
mcmc_number_of_fathers <- apply(mcmc_paternity$paternity,2,function(x) length(unique(x)))

library(ggplot2)
ggplot(data.frame(x=mcmc_number_of_fathers)) + theme_bw() + 
  geom_histogram(aes(x=mcmc_number_of_fathers, y=..count../length(mcmc_number_of_fathers)), binwidth=1, fill="white", color="black") +
  geom_vline(xintercept=length(unique(estimated_paternity$paternity)), color="red") +
  xlim(0,20) + xlab("Estimated number of fathers for colony") + ylab("Posterior probability")
