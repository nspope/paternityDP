simulate_sibling_group <- function(number_of_offspring_per_father, 
                                   allele_frequencies_per_msat,
                                   rate_of_allelic_dropout_per_locus,
                                   rate_of_allelic_mistyping_per_locus,
                                   probability_of_missing_data)
{
  number_of_loci <- length(allele_frequencies_per_msat)
  number_of_msats <- length(allele_frequencies_per_msat)
  number_of_alleles_per_msat <- sapply(allele_frequencies_per_msat, length)
  number_of_fathers <- length(number_of_offspring_per_father)
  number_of_mothers <- 1 #for now, restrict to single mother per sib-group
  maternal_ploidy <- 2
  paternal_ploidy <- 1
  offspring_ploidy <- 2

  # reparameterize error rates to probabilities 
  probability_of_allelic_dropout_per_locus <- rate_of_allelic_dropout_per_locus
  probability_of_allelic_mistyping_per_locus <- rate_of_allelic_mistyping_per_locus

  # check inputs
  number_of_offspring <- sum(number_of_offspring_per_father)
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
  offspring_paternity[] <- rep(1:number_of_fathers, number_of_offspring_per_father)
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
      if (runif(1) < probability_of_missing_data) 
      { 
        genotype[] <- NA 
        return(genotype)
      }
      dropouts <- 0
      mistypes <- 0
      if (length(unique(genotype)) > 1)
      {
        if (runif(1) < 2*probability_of_allelic_dropout) 
        { 
          genotype[] <- sample(genotype, size=1) 
          dropouts <- dropouts+1
        }
      }
      for (allele in 1:length(genotype)) 
      { 
        if (length(possible_alleles) > 1)
        {
          if (runif(1) < probability_of_allelic_mistyping) 
          { 
            genotype[allele] <- sample(possible_alleles[possible_alleles != genotype[allele]], size=1) 
            mistypes <- mistypes+1
          } 
        }
      }
      attr(genotype, "errors") <- c(dropouts,mistypes)
      return(genotype)
    }

    errors <- c(0,0)
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
      observed_maternal_genotypes[,mother,locus] <- draw <-
        simulate_genotyping_errors(genotype = true_maternal_genotypes[,mother,locus], 
                                   possible_alleles = names(allele_frequencies_per_msat[[locus]]), 
                                   probability_of_allelic_dropout = probability_of_allelic_dropout_per_locus[locus], 
                                   probability_of_allelic_mistyping = probability_of_allelic_mistyping_per_locus[locus], 
                                   probability_of_missing_data = probability_of_missing_data)
      if (!any(is.na(draw))) errors <- errors + attr(draw, "errors")
    }
    for (sib in 1:number_of_offspring) 
    { 
      observed_offspring_genotypes[,sib,locus] <- draw <-
        simulate_genotyping_errors(genotype = true_offspring_genotypes[,sib,locus], 
                                   possible_alleles = names(allele_frequencies_per_msat[[locus]]), 
                                   probability_of_allelic_dropout = probability_of_allelic_dropout_per_locus[locus], 
                                   probability_of_allelic_mistyping = probability_of_allelic_mistyping_per_locus[locus], 
                                   probability_of_missing_data = probability_of_missing_data)
      if (!any(is.na(draw))) errors <- errors + attr(draw, "errors")
    }
  }
  names(errors) <- c("dropouts", "mistypes")

  list("true_offspring_genotypes"=true_offspring_genotypes,
       "observed_offspring_genotypes"=observed_offspring_genotypes,
       "true_maternal_genotypes"=true_maternal_genotypes,
       "observed_maternal_genotypes"=observed_maternal_genotypes,
       "true_paternal_genotypes"=true_paternal_genotypes,
       "observed_paternal_genotypes"=observed_paternal_genotypes,
       "offspring_paternity"=offspring_paternity,
       "offspring_maternity"=offspring_maternity,
       "observed_number_of_fathers"=length(unique(offspring_paternity)),
       "errors"=errors
       )
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
        cat(dimnames(genotype_array)[[2]][i],"\t",sep="")
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
  out_raw <- as.numeric(out)
  out <- array(out_raw, dim(out))
  dimnames(out) <- list(paste0("allele",1:2), rownames(tb), colnames(tb))
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

plot_genotyping_errors <- function(mcmc_paternity, my_genotype_data) 
{
  library(ggplot2)

  mistyping_errors <- mcmc_paternity$mistyping_errors
  rownames(mistyping_errors) <- dimnames(my_genotype_data)[[2]]
  colnames(mistyping_errors) <- dimnames(my_genotype_data)[[3]]
  mistyping_errors <- as.data.frame.table(mistyping_errors)
  colnames(mistyping_errors) <- c("Locus", "Sample", "Errors")
  mistyping_errors$type <- "Mistyping"

  dropout_errors <- mcmc_paternity$dropout_errors
  rownames(dropout_errors) <- dimnames(my_genotype_data)[[2]]
  colnames(dropout_errors) <- dimnames(my_genotype_data)[[3]]
  dropout_errors <- as.data.frame.table(dropout_errors)
  colnames(dropout_errors) <- c("Locus", "Sample", "Errors")
  dropout_errors$type <- "Dropout"

  ggplot(rbind(mistyping_errors,dropout_errors)) +
    geom_tile(aes(x=Sample, y=Locus, fill=Errors)) +
    theme_bw() + 
    theme(text=element_text(family="Garamond"), 
          axis.text.x=element_text(angle=90),
          panel.grid=element_blank(), legend.position="bottom") +
    xlab("Locus") + ylab("Sample") + 
    guides(fill=guide_colorbar("E[# Errors]", direction="horizontal", title.position="top")) +
    scale_fill_gradient(low="white", high="firebrick") +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    facet_grid(~type)
}

plot_posterior_number_of_fathers <- function(mcmc_paternity)
{
  library(ggplot2)
  mcmc_number_of_fathers <- mcmc_paternity$number_of_fathers
  ggplot(data.frame(x=mcmc_number_of_fathers)) + theme_bw() + theme(panel.grid=element_blank()) +
    geom_histogram(aes(x=mcmc_number_of_fathers, y=..count../length(mcmc_number_of_fathers)), binwidth=1, fill="gray90", color="black") +
    xlim(0,20) + xlab("Estimated number of fathers for colony") + ylab("Posterior probability")
}

plot_posterior <- function(mcmc_paternity)
{
  library(ggplot2)
  ggplot(data.frame(n_fathers=mcmc_paternity$number_of_fathers, mistyping_rate=mcmc_paternity$mistyping_rate[1,])) +
    geom_bin2d(aes(x=mistyping_rate,y=n_fathers,fill=..count../1000), binwidth=c(0.005,1)) +
    theme_bw() + 
    theme(text=element_text(family="Garamond"), panel.grid=element_blank(), legend.position="bottom") +
    xlab("Proportion of erroneous phenotypes") +
    ylab("Number of fathers in sample") +
    guides(fill=guide_colorbar("Posterior probability", direction="horizontal", title.position="top")) +
    scale_x_continuous(expand=c(0,0), limits=c(0,0.3)) +
    scale_y_continuous(expand=c(0,0), limits=c(-1,15),breaks=seq(0.5,14.5,1),labels=1:15)
}

add_unsampled_to_phenotype_array <- function(phenotypes, fathers = 0, mothers = 0, offspring = 0)
{
  new_phenotypes <- array(0, dim=c(dim(phenotypes)[1], dim(phenotypes)[2]+fathers+mothers+offspring, dim(phenotypes)[3]))
  sample_names <- dimnames(phenotypes)[[2]]
  if (fathers>0) sample_names <- c(sample_names, paste0("add_father", 1:fathers))
  if (mothers>0) sample_names <- c(sample_names, paste0("add_mother", 1:mothers))
  if (offspring>0) sample_names <- c(sample_names, paste0("add_offspring", 1:offspring))
  dimnames(new_phenotypes) <- list(dimnames(phenotypes)[[1]], sample_names, dimnames(phenotypes)[[3]])
  for (i in 1:dim(phenotypes)[2]) new_phenotypes[,i,] <- phenotypes[,i,]
  new_phenotypes
}

cross_validation_for_number_of_parents <- 
  function(phenotypes, 
           sampled_mothers=1, #indices
           max_fathers=1, 
           max_mothers=0,
           nmcmc=500, nburnin=500, nthin=20)
{
  if (max_mothers < length(sampled_mothers)) stop("max_mothers < length(sampled_mothers)")
  if (max_mothers < 1 | max_fathers < 1) stop("max_mothers < 1 | max_fathers < 1")

  loo_cv <- data.frame()
  raw_lpd <- list()
  nmothers <- max_mothers
  nfathers <- max_fathers
  start <- if (length(sampled_mothers) == 0) 1 else length(sampled_mothers)
  for(m in start:nmothers)
  {
    for(f in 1:nfathers)
    {
      cat(m, " ", f, "\n", sep="")
      raw_lpd[[paste0(m,"_",f)]] <- c()
      new_phenotypes <- add_unsampled_to_phenotype_array(phenotypes, mothers=m-length(sampled_mothers), fathers=f)
      offspring <- 1:dim(phenotypes)[2]
      if (length(sampled_mothers)>0) offspring <- offspring[-sampled_mothers]
      for (i in offspring)
      {
        fit <- sample_parentage_and_error_rates_from_joint_posterior(
                 new_phenotypes,
                 c(sampled_mothers, grep("add_mother", dimnames(new_phenotypes)[[2]])),
                 grep("add_father", dimnames(new_phenotypes)[[2]]),
                 c(i),
                 number_of_mcmc_samples = nmcmc,
                 burn_in_samples = nburnin,
                 thinning_interval = nthin)
        raw_lpd[[paste0(m,"_",f)]] <- rbind(raw_lpd[[paste0(m,"_",f)]], fit$holdout_deviance[i,])
      }
      elpd <- rowMeans(raw_lpd[[paste0(m,"_",f)]])
      loo_cv <- rbind(loo_cv, data.frame(elpd=elpd, sample=offspring, fathers=f, mothers=m))
    }
  }
  average_scores <- aggregate(loo_cv$elpd, by=list(fathers=loo_cv$fathers, mothers=loo_cv$mothers), mean)
  best_fathers <- average_scores$fathers[which.max(average_scores$x)]
  best_mothers <- average_scores$mothers[which.max(average_scores$x)]
  new_phenotypes <- add_unsampled_to_phenotype_array(phenotypes, mothers=best_mothers-length(sampled_mothers), fathers=best_fathers, offspring=1) 
  fit <- sydneyPaternity:::sample_parentage_and_error_rates_from_joint_posterior(
           new_phenotypes,
           c(sampled_mothers, grep("add_mother", dimnames(new_phenotypes)[[2]])),
           grep("add_father", dimnames(new_phenotypes)[[2]]),
           grep("add_offspring", dimnames(new_phenotypes)[[2]]),
           number_of_mcmc_samples = nmcmc,
           burn_in_samples = nburnin,
           thinning_interval = nthin)
  list(loo_cv=loo_cv, raw_lpd=raw_lpd, top_model=fit, dim_names=dimnames(new_phenotypes))
}

plot_cross_validation_scores <- function(cv_output)
{
  library(ggplot2)
  cv_output$loo_cv$fathers <- factor(cv_output$loo_cv$fathers)
  cv_output$loo_cv$mothers <- factor(cv_output$loo_cv$mothers)
  ggplot(cv_output$loo_cv) + 
    stat_summary(geom="pointrange", aes(x=mothers, y=elpd, group=fathers, color=fathers), position=position_dodge(width=0.3)) +
    theme_bw() + theme(panel.grid=element_blank(), text=element_text(family="Garamond"), strip.background=element_blank()) + 
    ylab("Average leave-one-out cross-validation score\n(higher is better)") + 
    xlab("# mothers in model") + labs(color="# fathers")
}

plot_cross_validation_traces <- function(cv_scores)
{
  library(ggplot2)
  id <- strsplit(names(cv_scores$raw_lpd), split="_")
  out <- c()
  for(i in 1:length(cv_scores$raw_lpd))
  {
    tmp <- data.frame(elpd=colMeans(cv_scores$raw_lpd[[i]]))
    tmp$rep <- 1:nrow(tmp)
    tmp$mothers <- paste0("mothers=",id[[i]][1])
    tmp$fathers <- paste0("fathers=",id[[i]][2])
    out <- rbind(out, tmp)
  }
  ggplot(out, aes(x=rep, y=elpd)) + geom_point(alpha=0.1, shape=19) + facet_wrap(fathers~mothers) + 
    theme_bw() + theme(panel.grid=element_blank(), text=element_text(family="Garamond"), strip.background=element_blank()) + 
    geom_smooth(method="loess", alpha=0.5, size=0.5, formula=y~x) +
    ylab("Average predictive score") + xlab("MCMC iteration")
}

plot_cross_validation_parentage <- function(cv_output)
{
  library(ggplot2)
  model <- cv_output$top_model
  offspring <- which(!is.na(cv_output$top_model$maternity[,1]))
  names(offspring) <- cv_output$dim_names[[2]][offspring]

  maternity_matrix <- matrix(0, length(offspring), length(offspring))
  for(i in 1:ncol(model$maternity))
  {
    maternity_matrix <- maternity_matrix + paternity_vector_to_adjacency_matrix(model$maternity[offspring,i])
  }
  maternity_matrix <- maternity_matrix / ncol(model$maternity)
  rownames(maternity_matrix) <- colnames(maternity_matrix) <- names(offspring)

  paternity_matrix <- matrix(0, length(offspring), length(offspring))
  for(i in 1:ncol(model$paternity))
  {
    paternity_matrix <- paternity_matrix + paternity_vector_to_adjacency_matrix(model$paternity[offspring,i])
  }
  paternity_matrix <- paternity_matrix / ncol(model$paternity)
  rownames(paternity_matrix) <- colnames(paternity_matrix) <- names(offspring)

  ord <- colnames(paternity_matrix)[hclust(1-as.dist(paternity_matrix))$order]

  df_paternity <- as.data.frame(which(!is.na(paternity_matrix), arr.ind=TRUE))
  df_paternity$probability <- c(paternity_matrix)
  df_paternity$row <- colnames(paternity_matrix)[df_paternity$row]
  df_paternity$col <- colnames(paternity_matrix)[df_paternity$col]
  df_paternity$type <- 'Share father'
  df_maternity <- as.data.frame(which(!is.na(maternity_matrix), arr.ind=TRUE))
  df_maternity$probability <- c(maternity_matrix)
  df_maternity$row <- colnames(maternity_matrix)[df_maternity$row]
  df_maternity$col <- colnames(maternity_matrix)[df_maternity$col]
  df_maternity$type <- 'Share mother'

  df <- rbind(df_maternity, df_paternity)
  df$row <- factor(df$row, levels=ord)
  df$col <- factor(df$col, levels=ord)

  ggplot(df) +
    geom_tile(aes(y=row, x=col, fill=probability), color="black") +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradient(low="white", high="dodgerblue") +
    theme_bw() + facet_grid(~type) + 
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), axis.title=element_blank(), strip.background=element_blank(), text=element_text(family="Garamond"))
}

remove_loci_with_excessive_missingness <- function(genotypes, threshold)
{
  prop_missing <- apply(genotypes, 3, function(x) sum(is.na(x) | x == 0)/length(c(x)))
  cat("Dropped ", sum(prop_missing > threshold), " loci\n", sep="")
  genotypes[,,prop_missing <= threshold]
}

remove_samples_with_excessive_missingness <- function(genotypes, threshold, always_keep = NA)
{
  prop_missing <- apply(genotypes, 2, function(x) sum(is.na(x) | x == 0)/length(c(x)))
  if (!is.na(always_keep)) retain <- grepl(always_keep, dimnames(genotypes)[[2]]) else retain <- rep(FALSE, length(prop_missing))
  cat("Dropped ", sum(prop_missing > threshold), " samples\n", sep="")
  genotypes[,prop_missing <= threshold | retain,]
}

remove_monomorphic_loci <- function(genotypes)
{
  num_alleles <- apply(genotypes, 3, function(x) length(unique(x[x!=0 & !is.na(x)])))
  cat("Dropped ", sum(num_alleles < 2), " loci\n", sep="")
  genotypes[,,num_alleles > 1]
}

fit_dp_model <- #doesnt work delete at some point
  function(phenotypes, 
           sampled_mothers=1, #indices
           nmcmc=500, nburnin=500, nthin=20, sample_from_prior=FALSE, concentration = 1, random_initialization = FALSE)
{
  models <- list()
  offspring <- 1:dim(phenotypes)[2]
  if (length(sampled_mothers)>0) offspring <- offspring[-sampled_mothers]
  new_phenotypes <- add_unsampled_to_phenotype_array(phenotypes, mothers=length(offspring), fathers=length(offspring), offspring=0)
  fathers <- grep("add_father", dimnames(new_phenotypes)[[2]])
  mothers <- c(sampled_mothers, grep("add_mother", dimnames(new_phenotypes)[[2]]))
  fit <- sample_parentage_and_error_rates_from_joint_posterior_alt(
                 new_phenotypes,
                 mothers = mothers,
                 fathers = fathers,
                 concentration = concentration,
                 number_of_mcmc_samples = nmcmc,
                 burn_in_samples = nburnin,
                 thinning_interval = nthin,
                 sample_from_prior = sample_from_prior,
                 random_initialization = random_initialization)
  fit$sample_names <- dimnames(new_phenotypes)[[2]]
  fit$imputed_genotypes <- get_imputed_genotypes(fit)
  dimnames(fit$imputed_genotypes) <- dimnames(new_phenotypes)
  fit$paternal_genotypes <- fit$imputed_genotypes[,fathers,]
  fit$maternal_genotypes <- fit$imputed_genotypes[,mothers,]
  fit
}

get_imputed_genotypes <- function(model) #useful for previous versions
{
  genotypes <- array(NA, c(2,nrow(model$paternity),length(model$genotypes)))
  for(locus in 1:length(model$genotypes)){
  allele_names <- model$allele_lengths[[locus]]
  genotypes[,,locus] <- apply(model$genotypes[[locus]], 3, function(x) allele_names[which(x==max(x),arr.ind=TRUE)[1,]])
  }
  genotypes
}

#--------------------------------------------------- simulation w/ fixed distance

simulate_sibling_group_with_fixed_parental_genotypes <- 
  function(maternal_genotypes,
           paternal_genotypes,
           number_of_offspring, 
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
  stopifnot(dim(paternal_genotypes)[1] == 1 & 
            dim(paternal_genotypes)[2] == number_of_fathers &
            dim(paternal_genotypes)[3] == number_of_loci) 
  stopifnot(dim(maternal_genotypes)[1] == 2 & 
            dim(maternal_genotypes)[2] == number_of_mothers &
            dim(maternal_genotypes)[3] == number_of_loci) 
            
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
      true_paternal_genotypes[,father,locus] <- paternal_genotypes[,father,locus] 
    }
    for (mother in 1:number_of_mothers)
    {
      true_maternal_genotypes[,mother,locus] <- maternal_genotypes[,mother,locus]
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
      if (runif(1) < probability_of_missing_data) 
      { 
        genotype[] <- NA 
        return(genotype)
      }
      dropouts <- 0
      mistypes <- 0
      if (length(unique(genotype)) > 1)
      {
        if (runif(1) < 2*probability_of_allelic_dropout) 
        { 
          genotype[] <- sample(genotype, size=1) 
          dropouts <- dropouts+1
        }
      }
      for (allele in 1:length(genotype)) 
      { 
        if (length(possible_alleles) > 1)
        {
          if (runif(1) < probability_of_allelic_mistyping) 
          { 
            genotype[allele] <- sample(possible_alleles[possible_alleles != genotype[allele]], size=1) 
            mistypes <- mistypes+1
          } 
        }
      }
      attr(genotype, "errors") <- c(dropouts,mistypes)
      return(genotype)
    }

    errors <- c(0,0)
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
      observed_maternal_genotypes[,mother,locus] <- draw <-
        simulate_genotyping_errors(genotype = true_maternal_genotypes[,mother,locus], 
                                   possible_alleles = names(allele_frequencies_per_msat[[locus]]), 
                                   probability_of_allelic_dropout = probability_of_allelic_dropout_per_locus[locus], 
                                   probability_of_allelic_mistyping = probability_of_allelic_mistyping_per_locus[locus], 
                                   probability_of_missing_data = probability_of_missing_data)
      if (!any(is.na(draw))) errors <- errors + attr(draw, "errors")
    }
    for (sib in 1:number_of_offspring) 
    { 
      observed_offspring_genotypes[,sib,locus] <- draw <-
        simulate_genotyping_errors(genotype = true_offspring_genotypes[,sib,locus], 
                                   possible_alleles = names(allele_frequencies_per_msat[[locus]]), 
                                   probability_of_allelic_dropout = probability_of_allelic_dropout_per_locus[locus], 
                                   probability_of_allelic_mistyping = probability_of_allelic_mistyping_per_locus[locus], 
                                   probability_of_missing_data = probability_of_missing_data)
      if (!any(is.na(draw))) errors <- errors + attr(draw, "errors")
    }
  }
  names(errors) <- c("dropouts", "mistypes")

  list("true_offspring_genotypes"=true_offspring_genotypes,
       "observed_offspring_genotypes"=observed_offspring_genotypes,
       "true_maternal_genotypes"=true_maternal_genotypes,
       "observed_maternal_genotypes"=observed_maternal_genotypes,
       "true_paternal_genotypes"=true_paternal_genotypes,
       "observed_paternal_genotypes"=observed_paternal_genotypes,
       "offspring_paternity"=offspring_paternity,
       "offspring_maternity"=offspring_maternity,
       "observed_number_of_fathers"=length(unique(offspring_paternity)),
       "errors"=errors
       )
}

simulate_mixed_colony <- function(number_of_offspring, new_mother_distance, new_father_distance, new_parents_proportion, allele_frequencies_per_msat, rate_of_allelic_dropout, rate_of_allelic_mistyping, probability_of_missing_data)
{
  num_loci <- length(allele_frequencies_per_msat)

  # simulate core colony
  offspring_in_core_colony <- floor(number_of_offspring * new_parents_proportion)
  core_colony <- simulate_sibling_group (number_of_offspring_per_father = c(offspring_in_core_colony), 
                          allele_frequencies_per_msat = allele_frequencies_per_msat,
                          rate_of_allelic_dropout_per_locus = rep(rate_of_allelic_dropout, num_loci),
                          rate_of_allelic_mistyping_per_locus = rep(rate_of_allelic_mistyping, num_loci),
                          probability_of_missing_data = probability_of_missing_data)

  # simulate extra parents genotypes
  maternal_genotype <- apply(core_colony$true_maternal_genotypes, c(2,3), sort)
  extramaternal_genotype <- maternal_genotype
  maxit <- 1000
  if (new_mother_distance > 0){
    diffs <- sample(1:length(maternal_genotype), new_mother_distance)
    todo <- TRUE
    it <- 1
    while (todo)
    {
      dummy <- array(NA, dim(maternal_genotype))
      dummy[diffs] <- "1"
      location <- which(!is.na(dummy), arr.ind=TRUE)
      for (i in 1:nrow(location))
      {
        l <- location[i,]
        extramaternal_genotype[l[1],l[2],l[3]] <- sample(names(allele_frequencies_per_msat[[l[3]]]), 1)
      }
      extramaternal_genotype <- apply(extramaternal_genotype, c(2,3), sort)
      if (sum(extramaternal_genotype != maternal_genotype) == new_mother_distance) todo <- FALSE
      it <- it + 1
      if (it > maxit) stop("rejection sampling failed, too few unique alleles for requested distance?")
    }
  }

  # simulate extra father
  paternal_genotype <- core_colony$true_paternal_genotypes
  extrapaternal_genotype <- paternal_genotype
  maxit <- 1000
  if (new_father_distance > 0){
    diffs <- sample(1:length(paternal_genotype), new_father_distance)
    todo <- TRUE
    it <- 1
    while (todo)
    {
      dummy <- array(NA, dim(paternal_genotype))
      dummy[diffs] <- "1"
      location <- which(!is.na(dummy), arr.ind=TRUE)
      for (i in 1:nrow(location))
      {
        l <- location[i,]
        extrapaternal_genotype[l[1],l[2],l[3]] <- sample(names(allele_frequencies_per_msat[[l[3]]]), 1)
      }
      if (sum(extrapaternal_genotype != paternal_genotype) == new_father_distance) todo <- FALSE
      it <- it + 1
      if (it > maxit) stop("rejection sampling failed, too few unique alleles for requested distance?")
    }
  }

  # simulate extra sib group
  offspring_in_extra_colony <- number_of_offspring - offspring_in_core_colony
  extra_colony <- simulate_sibling_group_with_fixed_parental_genotypes (
           maternal_genotypes = extramaternal_genotype,
           paternal_genotypes = extrapaternal_genotype,
           number_of_offspring = offspring_in_extra_colony, 
           proportion_of_sperm_per_father = 1, 
           allele_frequencies_per_msat = allele_frequencies_per_msat,
           rate_of_allelic_dropout_per_locus = rep(rate_of_allelic_dropout,num_loci),
           rate_of_allelic_mistyping_per_locus = rep(rate_of_allelic_mistyping,num_loci),
           probability_of_missing_data = probability_of_missing_data)

  # combine two groups; queen genotype is omitted in extra colony
  library(abind)
  core_colony <- abind(core_colony$observed_maternal_genotypes, core_colony$observed_offspring_genotypes, along=2)
  extra_colony <- extra_colony$observed_offspring_genotypes
  dimnames(core_colony)[[2]] <- paste0("core_", dimnames(core_colony)[[2]])
  dimnames(extra_colony)[[2]] <- paste0("extra_", dimnames(extra_colony)[[2]])
  combined <- abind(core_colony, extra_colony, along=2)
  combined <- array(as.numeric(combined[]), dim(combined), dimnames=dimnames(combined))
  maternal_genotype <- array(as.numeric(maternal_genotype[]), dim(maternal_genotype), dimnames=dimnames(maternal_genotype))
  paternal_genotype <- array(as.numeric(paternal_genotype[]), dim(paternal_genotype), dimnames=dimnames(paternal_genotype))
  extramaternal_genotype <- array(as.numeric(extramaternal_genotype[]), dim(extramaternal_genotype), dimnames=dimnames(extramaternal_genotype))
  extrapaternal_genotype <- array(as.numeric(extrapaternal_genotype[]), dim(extrapaternal_genotype), dimnames=dimnames(extrapaternal_genotype))
  attr(combined, "maternal_genotype") <- maternal_genotype[,,]
  attr(combined, "paternal_genotype") <- paternal_genotype[,,]
  attr(combined, "extramaternal_genotype") <- extramaternal_genotype[,,]
  attr(combined, "extrapaternal_genotype") <- extrapaternal_genotype[,,]
  combined
}

sample_parentage_with_multiple_chains <- function(phenotypes, mother = 1, number_of_chains = 3, burn_in = 0, number_of_mcmc_samples = 1000, thinning_interval = 1, maternity = NA, lambda_mother = 0., lambda_father = 0., alpha = 1., mistyping_rate = 0.01, dropout_rate = 0.01, update_error_rates = TRUE, update_allele_frequencies = TRUE)
{
  if (all(is.na(maternity))) maternity <- rep(0, ncol(phenotypes))
  fits <- lapply(1:number_of_chains, function(i) 
    sample_parentage_and_error_rates(phenotypes,maternity=maternity,mother=mother,
                                     burn_in=burn_in,number_of_mcmc_samples=number_of_mcmc_samples,
                                     concentration = alpha,
                                     lambda_mother = lambda_mother,
                                     lambda_father = lambda_father,
                                     starting_mistyping_rate = mistyping_rate,
                                     starting_dropout_rate = dropout_rate,
                                     update_error_rates = update_error_rates,
                                     update_allele_frequencies = update_allele_frequencies,
                                     thinning_interval=thinning_interval))
  attr(fits, "phenotypes") <- phenotypes
  fits
}

plot_trace <- function(list_of_models)
{
  library(ggplot2)
  df <- Reduce(rbind, lapply(1:length(list_of_models), function(i) data.frame(chain=i, iter=1:length(list_of_models[[i]]$deviance), loglik=-0.5*list_of_models[[i]]$deviance)))
  ggplot(df) + geom_point(aes(x=iter, y=loglik, color=factor(chain))) + theme_bw() + xlab("MCMC iteration") + ylab("Log probability of offspring phenotypes") + theme(legend.position="none")
}

plot_parentage <- function(list_of_models, return_consensus_maternity = FALSE, return_extramaternity_probabilities = FALSE)
{
  library(ggplot2)
  genotype_data <- attr(list_of_models, "phenotypes")
  offspring <- which(!is.na(list_of_models[[1]]$maternity[,1]))
  names(offspring) <- dimnames(genotype_data)[[2]][offspring]

  extramaternity_prob <- matrix(0, length(offspring), 1)
  maternity_matrix <- matrix(0, length(offspring), length(offspring))
  number_of_mothers <- c()
  for(j in 1:length(list_of_models))
  {
    for(i in 1:ncol(list_of_models[[j]]$maternity))
    {
      maternity_vector <- list_of_models[[j]]$maternity[offspring,i]
      maternity_matrix <- maternity_matrix + paternity_vector_to_adjacency_matrix(list_of_models[[j]]$maternity[offspring,i])
      extramaternity_prob <- extramaternity_prob + 1*(list_of_models[[j]]$maternity[offspring,i]>0)
      number_of_mothers <- c(number_of_mothers, length(unique(maternity_vector[!is.na(maternity_vector)])))
    }
  }
  extramaternity_prob <- extramaternity_prob / max(maternity_matrix)
  maternity_matrix <- maternity_matrix / max(maternity_matrix)
  rownames(maternity_matrix) <- colnames(maternity_matrix) <- names(offspring)
  rownames(extramaternity_prob) <- names(offspring)
  colnames(extramaternity_prob) <- c("Extramaternity")
  if (return_extramaternity_probabilities) return(extramaternity_prob)

  number_of_mothers <- table(number_of_mothers)
  best_number_of_mothers <- as.numeric(names(number_of_mothers[which.max(number_of_mothers)]))
  consensus_maternity <- matrix(cutree(hclust(as.dist(1-maternity_matrix)), best_number_of_mothers), length(offspring), 1)
  rownames(consensus_maternity) <- colnames(maternity_matrix); colnames(consensus_maternity) <- "Maternity"
  if (return_consensus_maternity) return(consensus_maternity)

  paternity_matrix <- matrix(0, length(offspring), length(offspring))
  for(j in 1:length(list_of_models))
  {
    for(i in 1:ncol(list_of_models[[j]]$paternity))
    {
      paternity_matrix <- paternity_matrix + paternity_vector_to_adjacency_matrix(list_of_models[[j]]$paternity[offspring,i])
    }
  }
  paternity_matrix <- paternity_matrix / max(paternity_matrix)
  rownames(paternity_matrix) <- colnames(paternity_matrix) <- names(offspring)

  ord <- colnames(paternity_matrix)[hclust(1-as.dist(paternity_matrix))$order]

  df_paternity <- as.data.frame(which(!is.na(paternity_matrix), arr.ind=TRUE))
  df_paternity$probability <- c(paternity_matrix)
  df_paternity$row <- rownames(paternity_matrix)[df_paternity$row]
  df_paternity$col <- colnames(paternity_matrix)[df_paternity$col]
  df_paternity$type <- 'Share father'
  df_maternity <- as.data.frame(which(!is.na(maternity_matrix), arr.ind=TRUE))
  df_maternity$probability <- c(maternity_matrix)
  df_maternity$row <- rownames(maternity_matrix)[df_maternity$row]
  df_maternity$col <- colnames(maternity_matrix)[df_maternity$col]
  df_maternity$type <- 'Share mother'
  df_extramaternity <- as.data.frame(which(!is.na(extramaternity_prob), arr.ind=TRUE))
  df_extramaternity$probability <- c(extramaternity_prob)
  df_extramaternity$row <- rownames(extramaternity_prob)[df_extramaternity$row]
  df_extramaternity$col <- colnames(extramaternity_prob)[df_extramaternity$col]
  df_extramaternity$type <- 'Floater'

  df <- rbind(df_maternity, df_paternity, df_extramaternity)
  df$row <- factor(df$row, levels=c(ord, "Extramaternity"))
  df$col <- factor(df$col, levels=c(ord, "Extramaternity"))

  ggplot(df) +
    geom_tile(aes(y=row, x=col, fill=probability), color="black") +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradient(low="white", high="dodgerblue", limits=c(0,1)) +
    theme_bw() + facet_grid(~type, scales="free_x", space="free") + 
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), axis.title=element_blank(), strip.background=element_blank(), text=element_text(family="Garamond"))
}

group_by_ibs <- function(phenotypes, mother, return_matrix=FALSE) #starting values for maternity vector
{
  shares <- matrix(0, dim(phenotypes)[2], dim(phenotypes)[2])
  counts <- matrix(0, dim(phenotypes)[2], dim(phenotypes)[2])
  for (i in 1:dim(phenotypes)[3])
    for (j1 in 1:dim(phenotypes)[2])
      for (j2 in 1:dim(phenotypes)[2])
      {
        a <- phenotypes[,j1,i]; b <- phenotypes[,j2,i]
        if (all(!is.na(c(a,b))))
        {
          counts[j1,j2] <- counts[j1,j2] + 2
          if (a[1] %in% b)
          {
            d <- which(b==a[1])
            b <- b[-d]
            shares[j1,j2] <- shares[j1,j2] + 1
          }
          if (a[2] %in% b)
          {
            shares[j1,j2] <- shares[j1,j2] + 1
          }
        }
      }
  ibs_matrix <- shares/counts
  if (return_matrix) return(ibs_matrix)
  init_maternity <- cutree(hclust(as.dist(1-ibs_matrix)),2)
  init_maternity[init_maternity == init_maternity[mother]] <- 0
  init_maternity2 <- init_maternity
  unq <- unique(init_maternity)
  for(i in 1:length(unq)) init_maternity2[init_maternity==unq[i]] <- i
  init_maternity2-1
  #make sure to have whatsit recode to contiguous integers
}
