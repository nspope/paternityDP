library(paternityDP)

genotype_data <- genotype_array_from_txt("colony2_genotypes.txt")

#drop loci w/ >50% missing loci
genotype_data <- remove_loci_with_excessive_missingness(genotype_data, 0.5) 

#drop samples w/ >50% missing, retain samples with "Qu" in name
genotype_data <- remove_samples_with_excessive_missingness(genotype_data, 0.5, always_keep = "Qu") 

#drop loci w/ 1 allele
genotype_data <- remove_monomorphic_loci(genotype_data)

#find queen index (sample name with "Qu" in it)
the_queen <- grep("Qu", colnames(genotype_data))

#fit model
model_fit <- sample_paternity_and_error_rates_from_joint_posterior(
  genotype_data, 
  mother=the_queen,
  number_of_mcmc_samples=1000,
  global_genotyping_error_rates=FALSE, #use same rate across loci?
  add_unsampled_allele=FALSE #include an unobserved allele per locus?
)

#outputs
model_fit$paternity #columns are MCMC samples, rows are individuals
model_fit$dropout_error #columns are MCMC samples, rows are loci
model_fit$mistyping_error #columns are MCMC samples, rows are loci

#posterior probabilities for number of paternities
table(model_fit$number_of_fathers) / length(model_fit$number_of_fathers)

