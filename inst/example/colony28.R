library(sydneyPaternity)

#---1. Load the genotype data into an array
genotype_data <- genotype_array_from_txt("colony28_genotypes.txt")
genotype_data <- remove_loci_with_excessive_missingness(genotype_data, 0.5) #drop loci w/ >50% missing loci
genotype_data <- remove_samples_with_excessive_missingness(genotype_data, 0.5, always_keep = "Qu") #drop samples w/ >50% missing, retain samples with "Qu" in name
genotype_data <- remove_monomorphic_loci(genotype_data) #drop loci w/ 1 allele
the_queen <- grep("Qu", dimnames(genotype_data)[[2]]) #find queen index (sample name with "Qu" in it)

#---2. Use cross-validation to check for multiple maternity (this can take quite awhile to run)
multimaternity_fit <- sample_parentage_with_multiple_chains(genotype_data, mother=the_queen, number_of_chains=3, thinning_interval=1, burn_in=0, number_of_mcmc_samples=1000)
plot_trace(multimaternity_fit) #check MCMC convergence
plot_parentage(multimaternity_fit) #we can divvy up "mixed" colonies based on this

#---3. Remove floaters
maternity <- plot_parentage(multimaternity_fit, return_consensus_maternity=TRUE) #PLEASE check that this makes sense wrt the parentage plot. if there's a lot of uncertainty regarding parentage, this prediction will be worthless
floater_prob <- plot_parentage(multimaternity_fit, return_extramaternity_probabilities=TRUE) #PLEASE check that this makes sense wrt the parentage plot.
cbind(maternity, floater_prob) #take a look, does this make sense? choose a threshold for line below
bees_to_remove <- rownames(floater_prob)[floater_prob >= 0.5] #make sure threshold makes sense
genotype_data_clean <- genotype_data[,!(colnames(genotype_data) %in% bees_to_remove),] #drop "bees_to_remove"

#---3. Refit model assuming that there's a single maternity (after perhaps removing some bees)
multipaternity_fit <- sample_paternity_and_error_rates_from_joint_posterior(genotype_data_clean, mother=grep("Qu",colnames(genotype_data_clean)))
posterior_prob_number_of_fathers <- table(multipaternity_fit$number_of_fathers)/length(multipaternity_fit$number_of_fathers)

# you can access the paternity information here
multipaternity_fit$paternity
colnames(genotype_data_clean) #correspond to rows in the above ^^^

# make a csv or table that looks like this:
#
# colony number_of_fathers posterior_prob
#     28                 1          0.905 
#     28                 2          0.090 
#     28                 3          0.005 

