simulate_sibling_group <- function(number_of_offspring, 
                                   proportion_of_sperm_per_father, 
                                   allele_frequencies_per_msat,
                                   rate_of_allelic_dropout_per_locus,
                                   rate_of_allelic_mistyping_per_locus,
                                   probability_of_missing_data)
{
  number_of_loci <- length(allele_frequencies_per_msat)
  number_of_msats <- length(allele_frequencies_per_msat)
  number_of_alleles_per_msat <- sapply(allele_frequencies_per_msat, length)
  number_of_fathers <- length(proportion_of_sperm_per_father)
  number_of_mothers <- 1 #for now, restrict to single mother per sib-group
  maternal_ploidy <- 2
  paternal_ploidy <- 1
  offspring_ploidy <- 2

  # reparameterize error rates to probabilities 
  probability_of_allelic_dropout_per_locus <- rate_of_allelic_dropout_per_locus
  probability_of_allelic_mistyping_per_locus <- rate_of_allelic_mistyping_per_locus

  # check inputs
  stopifnot(number_of_offspring >= 1)
  stopifnot(number_of_fathers >= 1)
  stopifnot(length(probability_of_allelic_dropout_per_locus) == number_of_loci)
  stopifnot(length(probability_of_allelic_mistyping_per_locus) == number_of_loci)

  # allocate output
  true_maternal_genotypes <- array(NA, dim=c(maternal_ploidy, number_of_mothers, number_of_loci), 
                                   dimnames=list(paste0("allele",1:maternal_ploidy),paste0("mother",1:number_of_mothers),paste0("locus",1:number_of_loci)))
  observed_maternal_genotypes <- array(NA, dim=c(maternal_ploidy, number_of_mothers, number_of_loci),
                                       dimnames=list(paste0("allele",1:maternal_ploidy),paste0("mother",1:number_of_mothers),paste0("locus",1:number_of_loci)))
  true_paternal_genotypes <- array(NA, dim=c(paternal_ploidy, number_of_fathers, number_of_loci),
                                   dimnames=list(paste0("allele",1:paternal_ploidy),paste0("father",1:number_of_fathers),paste0("locus",1:number_of_loci)))
  observed_paternal_genotypes <- array(NA, dim=c(paternal_ploidy, number_of_fathers, number_of_loci),
                                       dimnames=list(paste0("allele",1:paternal_ploidy),paste0("father",1:number_of_fathers),paste0("locus",1:number_of_loci)))
  true_offspring_genotypes <- array(NA, dim=c(offspring_ploidy, number_of_offspring, number_of_loci),
                                    dimnames=list(paste0("allele",1:offspring_ploidy),paste0("offspring",1:number_of_offspring),paste0("locus",1:number_of_loci)))
  observed_offspring_genotypes <- array(NA, dim=c(offspring_ploidy, number_of_offspring, number_of_loci),
                                        dimnames=list(paste0("allele",1:offspring_ploidy),paste0("offspring",1:number_of_offspring),paste0("locus",1:number_of_loci)))
  offspring_paternity <- rep(NA, number_of_offspring)
  offspring_maternity <- rep(1, number_of_offspring) #for now, restrict to single mother per sib-group
  names(offspring_paternity) <- names(offspring_maternity) <- paste0("offspring",1:number_of_offspring)

  # simulate paternity and maternity of offspring
  offspring_paternity[] <- sample(1:number_of_fathers,
                                  size=number_of_offspring,
                                  prob=proportion_of_sperm_per_father,
                                  replace=TRUE)
  offspring_maternity[] <- 1 #for now all sibs have same mother

  for (locus in 1:number_of_loci)
  {
    # if there aren't 'labels' associated with alleles, use integers
    if (is.null(names(allele_frequencies_per_msat[[locus]]))) names(allele_frequencies_per_msat[[locus]]) <- as.character(1:number_of_alleles_per_msat[locus])

    # simulate true parental genotypes
    for (father in 1:number_of_fathers)
    {
      true_paternal_genotypes[,father,locus] <- sample(names(allele_frequencies_per_msat[[locus]]),
                                                       size=paternal_ploidy, 
                                                       prob=allele_frequencies_per_msat[[locus]],
                                                       replace=TRUE)
    }
    for (mother in 1:number_of_mothers)
    {
      true_maternal_genotypes[,mother,locus] <- sample(names(allele_frequencies_per_msat[[locus]]),
                                                       size=maternal_ploidy, 
                                                       prob=allele_frequencies_per_msat[[locus]],
                                                       replace=TRUE)
    }

    # simulate true offspring genotypes, conditional on paternity/maternity and parental genotypes
    for (sib in 1:number_of_offspring)
    {
      true_offspring_genotypes[1,sib,locus] <- sample(true_maternal_genotypes[,offspring_maternity[sib],locus], size=1)
      true_offspring_genotypes[2,sib,locus] <- sample(true_paternal_genotypes[,offspring_paternity[sib],locus], size=1)
    }

    # in practice genotyping is prone to error. For microsats, this is modelled by
    # "allele dropout", where only one allele amplifies (observed genotype always appears homozygous)
    # "mistyping", where a given allele gets mistakenly genotyped as another
    # "missing data", where neither allele amplifies
    simulate_genotyping_errors <- function(genotype, 
                                           possible_alleles,
                                           probability_of_allelic_dropout,
                                           probability_of_allelic_mistyping,
                                           probability_of_missing_data)
    { 
      # this function replicates the error model used by COLONY, see (Wang 2004 Genetics) and (Wang 2018 Methods Ecology Evolution)
      # (the parameterization of error rates is a bit different to make it easier to read)
      if (runif(1) < 2*probability_of_allelic_dropout) 
      { 
        genotype[] <- sample(genotype, size=1) 
      }
      for (allele in 1:length(genotype)) 
      { 
        if (length(possible_alleles) > 1)
        {
          if (runif(1) < probability_of_allelic_mistyping) 
          { 
            genotype[allele] <- sample(possible_alleles[possible_alleles != genotype[allele]], size=1) 
          } 
        }
      }
      if (runif(1) < probability_of_missing_data) 
      { 
        genotype[] <- NA 
      }
      return(genotype)
    }
    for (father in 1:number_of_fathers) 
    { 
      observed_paternal_genotypes[,father,locus] <- 
        simulate_genotyping_errors(genotype = true_paternal_genotypes[,father,locus], 
                                   possible_alleles = names(allele_frequencies_per_msat[[locus]]), 
                                   probability_of_allelic_dropout = probability_of_allelic_dropout_per_locus[locus], 
                                   probability_of_allelic_mistyping = probability_of_allelic_mistyping_per_locus[locus], 
                                   probability_of_missing_data = probability_of_missing_data)
    }
    for (mother in 1:number_of_mothers)
    { 
      observed_maternal_genotypes[,mother,locus] <- 
        simulate_genotyping_errors(genotype = true_maternal_genotypes[,mother,locus], 
                                   possible_alleles = names(allele_frequencies_per_msat[[locus]]), 
                                   probability_of_allelic_dropout = probability_of_allelic_dropout_per_locus[locus], 
                                   probability_of_allelic_mistyping = probability_of_allelic_mistyping_per_locus[locus], 
                                   probability_of_missing_data = probability_of_missing_data)
    }
    for (sib in 1:number_of_offspring) 
    { 
      observed_offspring_genotypes[,sib,locus] <- 
        simulate_genotyping_errors(genotype = true_offspring_genotypes[,sib,locus], 
                                   possible_alleles = names(allele_frequencies_per_msat[[locus]]), 
                                   probability_of_allelic_dropout = probability_of_allelic_dropout_per_locus[locus], 
                                   probability_of_allelic_mistyping = probability_of_allelic_mistyping_per_locus[locus], 
                                   probability_of_missing_data = probability_of_missing_data)
    }
  }

  list("true_offspring_genotypes"=true_offspring_genotypes,
       "observed_offspring_genotypes"=observed_offspring_genotypes,
       "true_maternal_genotypes"=true_maternal_genotypes,
       "observed_maternal_genotypes"=observed_maternal_genotypes,
       "true_paternal_genotypes"=true_paternal_genotypes,
       "observed_paternal_genotypes"=observed_paternal_genotypes,
       "offspring_paternity"=offspring_paternity,
       "offspring_maternity"=offspring_maternity,
       "observed_number_of_fathers"=length(unique(offspring_paternity))
       )
}


perturb_paternity_vector <- function(paternity, number_of_changes = 1)
{
  unique_fathers <- unique(paternity)
  paternity <- match(paternity, unique_fathers)
  who_to_shift <- sample(1:length(paternity), number_of_changes)
  paternity[who_to_shift] <- sample(1:(length(unique_fathers)+1), number_of_changes, replace=TRUE)
  unique_fathers <- unique(paternity)
  paternity <- match(paternity, unique_fathers)-1
  paternity
}

genotype_array_to_txt <- function(genotype_array, filename)
{
  ploidy <- dim(genotype_array)[1]
  number_of_samples <- dim(genotype_array)[2]
  number_of_loci <- dim(genotype_array)[3]
  sink(filename)
  try({
    cat(paste0("locus",1:number_of_loci),"\n",sep="\t")
    for (i in 1:number_of_samples)
    {
      for (j in 1:number_of_loci)
      {
        cat(paste(genotype_array[,i,j],collapse="/"),"\t",sep="")
      }
      cat("\n")
    }
  })
  sink()
}

list_of_genotype_arrays_to_txt <- function(genotype_list, filename)
{
  number_of_loci <- unique(sapply(genotype_list, function(x) dim(x)[3]))  
  stopifnot(length(number_of_loci) == 1)
  sink(filename)
  try({
    cat(paste0("locus",1:number_of_loci),"\n",sep="\t")
    for(k in 1:length(genotype_list)){
      genotype_array <- genotype_list[[k]]
      ploidy <- dim(genotype_array)[1]
      number_of_samples <- dim(genotype_array)[2]
      for (i in 1:number_of_samples)
      {
        for (j in 1:number_of_loci)
        {
          cat(paste(genotype_array[,i,j],collapse="/"),"\t",sep="")
        }
        cat("\n")
      }
    }
  })
  sink()
}

genotype_array_from_txt <- function(filename, missing_data = "NA")
{
  tb <- read.table(filename, header=TRUE)
  number_of_samples <- nrow(tb)
  number_of_loci <- ncol(tb)
  ploidy <- length(strsplit(as.character(tb[1,1]), split="/")[[1]])
  out <- array(NA, c(ploidy,number_of_samples,number_of_loci))
  for(i in 1:number_of_samples)
  {
    for(j in 1:number_of_loci)
    {
      out[,i,j] <- strsplit(as.character(tb[i,j]),split="/")[[1]]
    }
  }
  out[out==missing_data] <- NA
  out
}

allele_frequencies_from_genotype_array <- function(genotypes)
{
  number_of_loci <- dim(genotypes)[3]
  frequencies <- list()
  for(i in 1:number_of_loci)
  {
    alleles <- unique(c(genotypes[,,i]))
    alleles <- alleles[!is.na(alleles)]
    frequencies[[i]] <- table(genotypes[,,i])/sum(!is.na(genotypes[,,i]))
  }
  frequencies
}

recode_genotype_array_to_integers <- function(genotypes, allele_frequencies = NULL)
{
  if (is.null(allele_frequencies)) allele_frequencies <- allele_frequencies_from_genotype_array(genotypes)
  number_of_loci <- dim(genotypes)[3]
  stopifnot(number_of_loci == length(allele_frequencies))
  recoded_genotypes <- array(NA, dim(genotypes))
  for(i in 1:number_of_loci)
  {
    alleles <- unique(c(genotypes[,,i]))
    alleles <- alleles[!is.na(alleles)]
    stopifnot(all(alleles %in% names(allele_frequencies[[i]])))
    recoded_genotypes[,,i] <- match(genotypes[,,i], alleles)-1
  }
  recoded_genotypes
}

paternity_vector_to_adjacency_matrix <- function(pat)
{
  1.*(as.matrix(dist(pat))==0)
}
