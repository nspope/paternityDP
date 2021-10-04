library(sydneyPaternity)

# 1. Load the genotype data we just saved into an array
my_genotype_data <- genotype_array_from_txt("colony2.txt")
my_genotype_data #look and make sure names are OK

# 2. Jointly estimate paternity and error rates
set.seed(1)
mcmc_paternity <- 
  sample_paternity_and_error_rates_from_joint_posterior(my_genotype_data,
                                                        mother = 1,
                                                        number_of_mcmc_samples=1000)

# 3. Visualize joint posterior distribution of error rates/paternity
plot_posterior(mcmc_paternity)

# 4. Visualize dropout and mistyping errors
plot_genotyping_errors(mcmc_paternity, my_genotype_data)
