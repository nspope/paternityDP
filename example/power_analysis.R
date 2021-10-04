library(sydneyPaternity)

#for reproducible results
set.seed(29800)

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
  expand.grid(replicate=1:10, #probably want more like 1:100 for actual analysis
              num_fathers=1:4, #range of fathers in colony
              error_rates_in_simulation=0.05, #used for simulating data
              error_rates_in_estimation=NA, #used for estimating from simulated data
              proportion_missing_data=0.1, #use something realistic
              number_of_offspring=20, 
              use_true_allele_freqs=FALSE) #if FALSE use uniform allele freqs for estimation

#run simulations
simulation_results <- data.frame()
simulation_results_raw <- list()
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
my_genotype_data_raw <- as.numeric(my_genotype_data) #TODO
my_genotype_data <- array(my_genotype_data_raw, dim(my_genotype_data))

# 3. Jointly estimate paternity and error rates
set.seed(1)
mcmc_paternity <- 
  sample_paternity_and_error_rates_from_joint_posterior(my_genotype_data, 
                                                        mother = 1,
                                                        number_of_mcmc_samples=1000)

simulation_results_raw[[i]] <- mcmc_paternity

posterior_probs <- table(mcmc_paternity$number_of_fathers)/1000
if (num_fathers %in% names(posterior_probs))
{
  posterior_prob_of_truth <- posterior_probs[names(posterior_probs)==num_fathers]
} else {
  posterior_prob_of_truth <- 0
}
if ("1" %in% names(posterior_probs))
{
  posterior_prob_of_monandry <- posterior_probs[names(posterior_probs)==1]
} else {
  posterior_prob_of_monandry <- 0
}

simulation_results <- rbind(simulation_results,
 data.frame(num_fathers=num_fathers, 
            replicate=replicate, 
            posterior_prob_of_truth=posterior_prob_of_truth,
            posterior_prob_of_monandry=posterior_prob_of_monandry,
            true_dropout_error=error_rates_in_simulation,
            true_mistyping_error=error_rates_in_simulation,
            est_dropout_error=mean(mcmc_paternity$dropout_rate),
            est_mistyping_error=mean(mcmc_paternity$mistyping_rate)))
}

save(simulation_results_raw, simulation_results, file="simulation_results.RData")


