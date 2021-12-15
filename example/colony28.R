library(sydneyPaternity)

#---1. Load the genotype data into an array
genotype_data <- genotype_array_from_txt("colony28_genotypes.txt")
genotype_data <- remove_loci_with_excessive_missingness(genotype_data, 0.5) #drop loci w/ >50% missing loci
genotype_data <- remove_samples_with_excessive_missingness(genotype_data, 0.5) #drop samples w/ >50% missing
genotype_data <- remove_monomorphic_loci(genotype_data) #drop loci w/ 1 allele
the_queen <- grep("Qu", dimnames(genotype_data)[[2]]) #find queen index (sample name with "Qu" in it)

#---2. Use cross-validation to check for multiple maternity (this can take quite awhile to run)
cv_scores <- cross_validation_for_number_of_parents(genotype_data, the_queen, max_fathers=5, max_mothers=5, nburnin=1000)
plot_cross_validation_scores(cv_scores) #shows predictive accuracy across different #'s of parents

#look at MCMC traces for CV scores; these should look more or less straight, and not have a step-like pattern
#if they do have a step-like pattern then its likely the CV scores haven't converged yet (==are biased)
#in this case, double the burnin samples in "cv_scores <- ..." above and check again, 
#if still not good double again, etc
plot_cross_validation_traces(cv_scores)
#this one^^^looks fine, an example of what you don't want it to look like is at the bottom of the script

#it took awhile so save scores
save(cv_scores, file="colony28_cv.RData")

#visualise parentage
plot_cross_validation_parentage(cv_scores)

#---here is an example of what you DON'T want to see in the trace plots
cv_scores_bad <- cross_validation_for_number_of_parents(genotype_data, the_queen, max_fathers=5, max_mothers=5, nburnin=0, nthin=1)
plot_cross_validation_traces(cv_scores_bad) 
#^^^notice step-like pattern, means MCMC sampler doesn't have sufficient burnin
