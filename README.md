```{r}
# 0. Simulate data
set.seed(1)
num_fathers <- 4
num_loci <- 4
num_allel <- 2
allele_freqs <- lapply(1:num_loci, function(x) rep(1/num_allel, num_allel))

out <- simulate_sibling_group(number_of_offspring=20,
                              proportion_of_sperm_per_father=rep(1/num_fathers,num_fathers),
                              allele_frequencies_per_msat=allele_freqs,
                              rate_of_allelic_dropout_per_locus=rep(0.1,num_loci),
                              rate_of_allelic_mistyping_per_locus=rep(0.1,num_loci),
                              probability_of_missing_data=0.1)

# 0. Here we output simulation to a file. Not necessary b/c simulation is already an R object, but this
#    allows us to write the rest of the example like you are working with your actual data.
#    The first row in the data is the maternal genotype, this is important later on.
list_of_genotype_arrays_to_txt(list(out$observed_maternal_genotype, out$observed_offspring_genotypes), 
                               filename = "example_data.txt")

# 1. Load the genotype data we just saved into an array
my_genotype_data <- genotype_array_from_txt("example_data.txt")

# 2. Calculate allele frequencies from genotype array
my_allele_freqs <- allele_frequencies_from_genotype_array(my_genotype_data)

# 3. Recode alleles to integers (uses positions in "allele_freqs" as integer codes)
my_genotype_data_recoded <- recode_genotype_array_to_integers(my_genotype_data, my_allele_freqs)

# 4. Pull out maternal genotype (see note above: if not first row in txt, you'll need to modify. Will want to split up colonies too. Ask me later)
mother_index <- 1
maternal_genotype <- my_genotype_data_recoded[,mother_index,]
offspring_genotypes <- my_genotype_data_recoded[,-mother_index,]

# 5. Find MLE of paternity vector
starting_paternity <- rep(1, dim(offspring_genotypes)[2]) #starting state for optimization, everybody has same father
dropout_rates <- rep(0.1, length(my_allele_freqs)) #for now just set to true values (used in simulations)
mistyping_rates <- rep(0.1, length(my_allele_freqs)) #for now just set to true values (used in simulations)
estimated_paternity <- optimize_paternity_given_error_rates(starting_paternity,offspring_genotypes,maternal_genotype,my_allele_freqs,dropout_rates,mistyping_rates)

# 6. Compare estimated and true paternities
cbind(estimated_paternity$paternity, out$offspring_paternity)
paternity_vector_to_adjacency_matrix(estimated_paternity$paternity)
paternity_vector_to_adjacency_matrix(estimated_paternity$paternity) == paternity_vector_to_adjacency_matrix(out$offspring_paternity)
c(length(unique(estimated_paternity$paternity)), length(unique(out$offspring_paternity)))

# 7. Get measure of uncertainty (not necessary for power analysis)
set.seed(1)
estimated_paternity_samples <- sample_paternity_given_error_rates(starting_paternity,offspring_genotypes,maternal_genotype,my_allele_freqs,dropout_rates,mistyping_rates)
number_of_fathers_samples <- apply(estimated_paternity_samples$paternity,2,function(x) length(unique(x)))
plot(as.numeric(names(table(number_of_fathers_samples))),
     unname(table(number_of_fathers_samples))/length(number_of_fathers_samples), xlim=c(0,20), ylim=c(0,1), type="h")
```
