#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <vector>

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]

arma::umat recode_to_contiguous_integers (arma::umat input)
{
  arma::uvec uniq = arma::unique(input);
  for(unsigned i=0; i<uniq.n_elem; ++i) input.replace(uniq.at(i), i);
  return input;
}

arma::uword sample (const arma::vec& pvec) 
{
  arma::uword K = pvec.n_elem;
  arma::uvec opts = arma::linspace<arma::uvec>(0L, K - 1L, K);
  return arma::conv_to<arma::uword>::from(
      Rcpp::RcppArmadillo::sample(opts, 1L, true, pvec)
  );
}

// [[Rcpp::export]]
double genotyping_error_model 
 (const arma::uvec& phenotype,
  const unsigned& genotype0, 
  const unsigned& genotype1, 
  const unsigned& number_of_alleles,
  const double& dropout_rate, 
  const double& mistyping_rate)
{
  // from Eqs 1 & 2 in Wang 2004 Genetics
  const double e1 = dropout_rate;
  const double e2 = mistyping_rate/double(number_of_alleles-1);
  const double E2 = mistyping_rate;

  const arma::uvec genotype = {genotype0, genotype1};
  const bool phenotype_is_homozygous = phenotype[0] == phenotype[1];
  const bool genotype_is_homozygous = genotype[0] == genotype[1];

  if (number_of_alleles == 1) return 1.; //monomorphic loci

  if (genotype_is_homozygous) 
  {
    if (phenotype_is_homozygous && phenotype[0] == genotype[0])
    {
      return std::pow(1.-E2, 2);
    } else if ((phenotype[0] == genotype[0] && phenotype[1] != genotype[0]) || 
               (phenotype[0] != genotype[0] && phenotype[1] == genotype[0]) ){
      return 2.*e2*(1-E2);
    } else if (phenotype[0] != genotype[0] && phenotype[1] != genotype[0]) {
      return (2.-int(phenotype_is_homozygous))*std::pow(e2, 2);
    }  
  } else {
    if ( (phenotype[0] == genotype[0] && phenotype[1] == genotype[1]) ||
         (phenotype[1] == genotype[0] && phenotype[0] == genotype[1]) )
    {
      return std::pow(1.-E2, 2) + std::pow(e2, 2) - 2.*e1*std::pow(1.-E2-e2, 2);
    } else if (phenotype_is_homozygous && (phenotype[0] == genotype[0] || phenotype[0] == genotype[1])) {
      return e2*(1.-E2) + e1*std::pow(1.-E2-e2, 2);
    } else if (phenotype[0] != genotype[0] && phenotype[0] != genotype[1] && 
               phenotype[1] != genotype[0] && phenotype[1] != genotype[1] ){
      return (2.-int(phenotype_is_homozygous))*std::pow(e2, 2);
    } else {
      // Wang 2018 has (1-E2-e2) and Wang 2004 has (1-E2+e2), latter is correct
      return e2*(1.-E2+e2);
    }
  }
  return 0.;
}

// [[Rcpp::export]]
int genotyping_error_model_class
 (const arma::uvec& phenotype,
  const unsigned& genotype0, 
  const unsigned& genotype1) 
{
  // from Eqs 1 & 2 in Wang 2004 Genetics
  const arma::uvec genotype = {genotype0, genotype1};
  const bool phenotype_is_homozygous = phenotype[0] == phenotype[1];
  const bool genotype_is_homozygous = genotype[0] == genotype[1];

  if (genotype_is_homozygous) 
  {
    if (phenotype_is_homozygous && phenotype[0] == genotype[0])
    {
      return 1;
    } else if ((phenotype[0] == genotype[0] && phenotype[1] != genotype[0]) || 
               (phenotype[0] != genotype[0] && phenotype[1] == genotype[0]) ){
      return 2;
    } else if (phenotype[0] != genotype[0] && phenotype[1] != genotype[0]) {
      return 3;
    }  
  } else {
    if ( (phenotype[0] == genotype[0] && phenotype[1] == genotype[1]) ||
         (phenotype[1] == genotype[0] && phenotype[0] == genotype[1]) )
    {
      return 4;
    } else if (phenotype_is_homozygous && (phenotype[0] == genotype[0] || phenotype[0] == genotype[1])) {
      return 5;
    } else if (phenotype[0] != genotype[0] && phenotype[0] != genotype[1] && 
               phenotype[1] != genotype[0] && phenotype[1] != genotype[1] ){
      return 6;
    } else {
      return 7;
    }
  }
  return 0;
}

// [[Rcpp::export]]
arma::uvec simulate_genotyping_errors
 (const arma::uvec& phenotype,
  const unsigned& genotype0, 
  const unsigned& genotype1, 
  const unsigned& number_of_alleles,
  const double& dropout_rate, 
  const double& mistyping_rate)
{
  const double e1 = dropout_rate;
  const double e2 = mistyping_rate/double(number_of_alleles-1);
  const double E2 = mistyping_rate;

  const arma::uvec genotype = {genotype0, genotype1};
  const bool phenotype_is_homozygous = phenotype[0] == phenotype[1];
  const bool genotype_is_homozygous = genotype[0] == genotype[1];

  if (number_of_alleles == 1) return arma::uvec({0,0}); //monomorphic loci

  if (genotype_is_homozygous) 
  {
    if (phenotype_is_homozygous && phenotype[0] == genotype[0])
    {
      return arma::uvec({0, 0});
    } else if ((phenotype[0] == genotype[0] && phenotype[1] != genotype[0]) || 
               (phenotype[0] != genotype[0] && phenotype[1] == genotype[0]) ){
      return arma::uvec({0, 1});
    } else if (phenotype[0] != genotype[0] && phenotype[1] != genotype[0]) {
      return arma::uvec({0, 2});
    }  
  } else {
    if ( (phenotype[0] == genotype[0] && phenotype[1] == genotype[1]) ||
         (phenotype[1] == genotype[0] && phenotype[0] == genotype[1]) )
    {
      // either: no errors occurred OR 
      //         there was no dropout error but two typing errors
      //         there was a dropout error that got fixed by single typing error
      // prop.table(c((1-2*e1)*(1-E2)^2, (1-2*e1)*e2*e2, 4*e1*e2*(1-E2)))
      // sum(c((1-2*e1)*(1-E2)^2, (1-2*e1)*e2*e2, 4*e1*e2*(1-E2)))
      // (1-E2)^2 + e2^2 - 2*e1*(1-E2-e2)^2
      // [0,0] [0,2] [1,1]
      const arma::vec probs = {(1-2*e1)*(1-E2)*(1-E2), (1-2*e1)*e2*e2, 4*e1*e2*(1-E2)};
      const arma::umat counts = {{0,0},{0,2},{1,1}};
      return counts.row(sample(probs)).t();
    } else if (phenotype_is_homozygous && (phenotype[0] == genotype[0] || phenotype[0] == genotype[1])) {
      // one class 2 error OR there was a dropout error and no typing errors
      // prop.table(c(e1*(1-E2)^2, (1-2*e1)*e2*(1-E2), e1*e2*e2))
      // sum(c(e1*(1-E2)^2, (1-2*e1)*e2*(1-E2), e1*e2*e2))
      // e2*(1-E2) + e1*(1-E2-e2)^2
      // [0 1] [1 0] [1 2]
      const arma::vec probs = {e1*(1-E2)*(1-E2), (1-2*e1)*e2*(1-E2), e1*e2*e2};
      const arma::umat counts = {{1,0},{0,1},{1,2}};
      return counts.row(sample(probs)).t();
    } else if (phenotype[0] != genotype[0] && phenotype[0] != genotype[1] && 
               phenotype[1] != genotype[0] && phenotype[1] != genotype[1] ){
      // two class 2 errors occurred, regardless of whether class 1 error occurs
      // prop.table(c( (1-2*e1)*e2*e2 , 2*e1*e2*e2 ))
      // [0 2] [1 2]
      const arma::vec probs = {(1-2*e1)*e2*e2, 2*e1*e2*e2};
      const arma::umat counts = {{0,2},{1,2}};
      return counts.row(sample(probs)).t();
    } else {
      // "otherwise" ... there's one match but phenotype is heterozygous?
      // 1. could have: sequencing error at one, no sequencing error at other
      // 2. sequencing error gets fixed by second sequencing error
      // dropout could happen in either case
      // prop.table(c( (1-2*e1)*e2*(1-E2), (1-2*e1)*e2*e2, 2*e1*e2*(1-E2), 2*e1*e2*e2 ))
      // [0 1] [0 2] [1 1] [1 2]
      const arma::vec probs = {(1-2*e1)*e2*(1-E2), (1-2*e1)*e2*e2, 2*e1*e2*(1-E2), 2*e1*e2*e2};
      const arma::umat counts = {{0,1},{0,2},{1,1},{1,2}};
      return counts.row(sample(probs)).t();
    }
  }
  return arma::uvec{{0,0}};
}

arma::umat sample_genotyping_errors_given_paternity
 (const arma::uvec& paternity,
  const arma::umat& offspring_phenotypes, 
  const arma::uvec& maternal_phenotype, 
  const arma::vec& allele_frequencies, 
  const double& dropout_rate, 
  const double& mistyping_rate)
{
  // simulate from conditional posterior of error events given phenotypes and paternity
  
  const unsigned number_of_alleles = allele_frequencies.n_elem;
  const arma::uvec fathers = arma::unique(paternity);
  const arma::vec allele_frequencies_normalized = allele_frequencies / arma::accu(allele_frequencies);

  if (offspring_phenotypes.n_rows != 2) Rcpp::stop("offspring phenotypes must have 2 rows");
  if (offspring_phenotypes.n_cols != paternity.n_elem) Rcpp::stop("offspring phenotypes must have column for each individual");
  if (maternal_phenotype.n_elem != 2) Rcpp::stop("maternal phenotype must have 2 elements");
  if (offspring_phenotypes.max() >= number_of_alleles) Rcpp::stop("offspring allele out of range");
  if (maternal_phenotype.max() >= number_of_alleles) Rcpp::stop("maternal allele out of range");
  if (arma::any(allele_frequencies_normalized < 0.)) Rcpp::stop("negative allele frequencies");
  if (dropout_rate <= 0. || mistyping_rate <= 0.) Rcpp::stop("negative genotyping error rates");

  // alternatively pass in as mutable argument
  arma::uvec maternal_genotype (2);
  arma::umat offspring_genotypes (2, paternity.n_elem);
  arma::uvec paternal_genotypes (fathers.n_elem);
  maternal_genotype.fill(arma::datum::nan);
  offspring_genotypes.fill(arma::datum::nan);
  paternal_genotypes.fill(arma::datum::nan);

  // simulate maternal genotype; there are choose(k,2)+k possible genotypes
  // marginalize over paternal & offspring genotypes to calculate posterior genotype probabilities
  unsigned number_of_genotypes = number_of_alleles*(number_of_alleles + 1)/2;
  arma::vec maternal_genotype_posterior (number_of_genotypes);
  arma::umat possible_maternal_genotypes (2, number_of_genotypes);
  unsigned genotype = 0;
  for (unsigned w=1; w<=number_of_alleles; ++w) // first maternal allele 
  { 
    for (unsigned v=w; v<=number_of_alleles; ++v) // second maternal allele
    {
      double maternal_genotype_probability = 
        (2.-int(w==v)) * allele_frequencies_normalized[w-1] * allele_frequencies_normalized[v-1]; //hwe prior
      double log_halfsib_likelihood = log(maternal_genotype_probability); 
      if (arma::prod(maternal_phenotype))
      {
        double maternal_phenotype_probability = 
          genotyping_error_model(maternal_phenotype, w, v, number_of_alleles, dropout_rate, mistyping_rate);
        log_halfsib_likelihood += log(maternal_phenotype_probability);
      }
      for (auto father : fathers)
      {
        double fullsib_likelihood = 0.;
        double running_maximum = -arma::datum::inf;
        arma::uvec offspring_from_father = arma::find(paternity == father);
        for (unsigned u=1; u<=number_of_alleles; ++u) // paternal allele
        {
          double paternal_genotype_probability = 
            allele_frequencies_normalized[u-1]; //hwe prior
          double log_fullsib_likelihood = log(paternal_genotype_probability);
          for (auto offspring : offspring_from_father)
          {
            arma::uvec offspring_phenotype = offspring_phenotypes.col(offspring);
            if (arma::prod(offspring_phenotype)) { 
              double offspring_phenotype_probability = // Mendelian segregation probs * phenotype probabilities
                0.5 * genotyping_error_model(offspring_phenotype, w, u, number_of_alleles, dropout_rate, mistyping_rate) + 
                0.5 * genotyping_error_model(offspring_phenotype, v, u, number_of_alleles, dropout_rate, mistyping_rate); 
              log_fullsib_likelihood += log(offspring_phenotype_probability);
            }
          }
          if (log_fullsib_likelihood <= running_maximum) //underflow protection
          {
            fullsib_likelihood += exp(log_fullsib_likelihood - running_maximum);
          } else {
            fullsib_likelihood *= exp(running_maximum - log_fullsib_likelihood);
            fullsib_likelihood += 1.;
            running_maximum = log_fullsib_likelihood;
          }
        }
        log_halfsib_likelihood += log(fullsib_likelihood) + running_maximum;
      }
      maternal_genotype_posterior.at(genotype) = log_halfsib_likelihood;
      possible_maternal_genotypes.col(genotype) = arma::uvec({w, v});
      genotype++;
    }
  }
  maternal_genotype_posterior -= maternal_genotype_posterior.max();
  maternal_genotype = 
    possible_maternal_genotypes.col(sample(arma::exp(maternal_genotype_posterior)));

  // simulate paternal genotype; there are k possible genotypes
  // paternal genotypes are conditionally independent with fixed maternal genotype
  // calculate posterior genotype probabilities by marginalizing over offspring genotypes
  arma::vec paternal_genotype_posterior (number_of_alleles);
  for (auto father : fathers)
  {
    arma::uvec offspring_from_father = arma::find(paternity == father);
    for (unsigned u=1; u<=number_of_alleles; ++u) // paternal allele
    {
      double paternal_genotype_probability = 
          allele_frequencies_normalized[u-1]; //hwe prior
      double log_fullsib_likelihood = log(paternal_genotype_probability);
      for (auto offspring : offspring_from_father)
      {
        arma::uvec offspring_phenotype = offspring_phenotypes.col(offspring);
        if (arma::prod(offspring_phenotype)) { 
          double offspring_phenotype_probability = // Mendelian segregation probs * phenotype probabilities
            0.5 * genotyping_error_model(offspring_phenotype, maternal_genotype[0], u, number_of_alleles, dropout_rate, mistyping_rate) + 
            0.5 * genotyping_error_model(offspring_phenotype, maternal_genotype[1], u, number_of_alleles, dropout_rate, mistyping_rate); 
          log_fullsib_likelihood += log(offspring_phenotype_probability);
        }
      }
      paternal_genotype_posterior[u-1] = log_fullsib_likelihood;
    }
    paternal_genotype_posterior -= paternal_genotype_posterior.max();
    paternal_genotypes.at(father) = sample(arma::exp(paternal_genotype_posterior));
  }

  // simulate offspring genotypes
  for (unsigned sib=0; sib<paternity.n_elem; ++sib)
  {
    arma::uvec offspring_phenotype = offspring_phenotypes.col(sib);
    arma::vec offspring_genotype_posterior (2);
    arma::umat possible_offspring_genotypes (2, 2);
    for (unsigned i=0; i<2; ++i)
    {
      offspring_genotype_posterior[i] = log(0.5); //Mendelian segregation
      if (arma::prod(offspring_phenotype)) { 
        offspring_genotype_posterior[i] += 
          log(genotyping_error_model(offspring_phenotype, maternal_genotype[i], 
                paternal_genotypes.at(paternity.at(sib)), number_of_alleles, dropout_rate, mistyping_rate));
      }
      possible_offspring_genotypes.col(i) = 
        arma::uvec({maternal_genotype[i], paternal_genotypes.at(paternity.at(sib))});
    }
    offspring_genotype_posterior -= offspring_genotype_posterior.max();
    offspring_genotypes.col(sib) = 
      possible_offspring_genotypes.col(sample(arma::exp(offspring_genotype_posterior)));
  }

  //simulate numbers of errors given genotypes and phenotypes
  //then new error rates via beta-binomial conjugacy, w/ beta(1,1) prior
  unsigned sampled_phenotypes = 0;
  unsigned sampled_heterozygotes = 0;
  arma::uvec counts_of_errors = {0,0}; //dropouts, mistypes
  //arma::umat error_counts (2, paternity.n_elem + 1, arma::fill::zeros); //could store this to see where errors are likely
  if (arma::prod(maternal_phenotype))
  {
    sampled_phenotypes++;
    if (maternal_genotype[0] != maternal_genotype[1]) sampled_heterozygotes++;
    //error_counts.col(paternity.n_elem + 1) = ...;
    counts_of_errors += 
      simulate_genotyping_errors(maternal_phenotype, maternal_genotype[0], 
          maternal_genotype[1], number_of_alleles, dropout_rate, mistyping_rate);
  }
  for (unsigned sib=0; sib<paternity.n_elem; ++sib)
  {
    arma::uvec offspring_phenotype = offspring_phenotypes.col(sib);
    if (arma::prod(offspring_phenotype)) { 
      sampled_phenotypes++;
      if (offspring_genotypes.at(0,sib) != offspring_genotypes.at(1,sib)) sampled_heterozygotes++;
      //error_counts.col(sib) = ...;
      counts_of_errors += 
        simulate_genotyping_errors(offspring_phenotype, offspring_genotypes.at(0,sib), 
          offspring_genotypes.at(1,sib), number_of_alleles, dropout_rate, mistyping_rate);
    }
  }

  ////now doing the rate sampling outside
  //double new_dropout_error = // a single dropout can occur when the genotype is heterozygous
  //  0.5 * R::rbeta(1. + counts_of_errors[0], 1. + sampled_heterozygotes - counts_of_errors[0]);
  //double new_mistyping_error = // both alleles can be mistyped per genotype regardless of zygosity
  //  R::rbeta(1. + counts_of_errors[1], 1. + 2.*sampled_phenotypes - counts_of_errors[1]);

  return arma::umat({{counts_of_errors[0], sampled_heterozygotes - counts_of_errors[0]},
                     {counts_of_errors[1], 2*sampled_phenotypes - counts_of_errors[1]}});
}

// [[Rcpp::export]]
Rcpp::List sample_error_rates_given_paternity
 (arma::uvec paternity, 
  arma::ucube offspring_phenotypes, 
  arma::umat maternal_phenotype,
  std::vector<arma::vec> allele_frequencies,
  arma::vec dropout_rate,
  arma::vec mistyping_rate,
  const unsigned max_iter = 1000,
  const unsigned global_genotyping_error_rates = false)
{
  //TODO make this generate AF prior, initialize stuff as for larger sampling routine
  
  paternity = recode_to_contiguous_integers(paternity);
  const unsigned number_of_loci = allele_frequencies.size();

  if (maternal_phenotype.n_cols != number_of_loci) Rcpp::stop("must have maternal phenotypes for each locus");
  if (offspring_phenotypes.n_slices != number_of_loci) Rcpp::stop("must have offspring phenotypes for each locus");
  if (dropout_rate.n_elem != number_of_loci) Rcpp::stop("must have dropout rates for each locus");
  if (mistyping_rate.n_elem != number_of_loci) Rcpp::stop("must have mistyping rates for each locus");

  arma::mat dropout_rate_samples (number_of_loci, max_iter);
  arma::mat mistyping_rate_samples (number_of_loci, max_iter);

  for (unsigned iter=0; iter<max_iter; ++iter)
  {
    arma::umat global_error_counts = arma::zeros<arma::umat>(2,2);
    for (unsigned locus=0; locus<number_of_loci; ++locus)
    {
      arma::umat error_counts =
        sample_genotyping_errors_given_paternity(paternity, offspring_phenotypes.slice(locus), 
            maternal_phenotype.col(locus), allele_frequencies[locus], dropout_rate[locus], mistyping_rate[locus]);
      if (global_genotyping_error_rates)
      {
        global_error_counts += error_counts;
      } else {
        dropout_rate[locus] = 0.5 * R::rbeta(1. + error_counts.at(0,0), 1. + error_counts.at(0,1));
        mistyping_rate[locus] = R::rbeta(1. + error_counts.at(1,0), 1. + error_counts.at(1,1));
      }
    }
    if (global_genotyping_error_rates)
    {
        dropout_rate.fill(0.5 * R::rbeta(1. + global_error_counts.at(0,0), 1. + global_error_counts.at(0,1)));
        mistyping_rate.fill(R::rbeta(1. + global_error_counts.at(1,0), 1. + global_error_counts.at(1,1)));
    }
    dropout_rate_samples.col(iter) = dropout_rate;
    mistyping_rate_samples.col(iter) = mistyping_rate;
    if (iter % 100 == 0) Rcpp::Rcout << "[" << iter << "]" << std::endl;
  }
  return Rcpp::List::create(
      Rcpp::_["dropout_rate"] = dropout_rate_samples,
      Rcpp::_["mistyping_rate"] = mistyping_rate_samples);
}

double paternity_loglikelihood_by_locus 
 (const arma::uvec& paternity,
  const arma::umat& offspring_phenotypes, 
  const arma::uvec& maternal_phenotype, 
  const arma::vec& allele_frequencies, 
  const double& dropout_rate, 
  const double& mistyping_rate)
{
  // likelihood of offspring paternity given offspring phenotypes, maternal phenotype, haplodiploidy
  // modified from Eqs 3 & 4 in Wang 2004 Genetics
  const unsigned number_of_alleles = allele_frequencies.n_elem;
  const unsigned number_of_genotypes = number_of_alleles*(number_of_alleles+1)/2;
  const arma::uvec fathers = arma::unique(paternity);
  const arma::vec allele_frequencies_normalized = allele_frequencies / arma::accu(allele_frequencies);

  if (offspring_phenotypes.n_rows != 2) Rcpp::stop("offspring phenotypes must have 2 rows");
  if (offspring_phenotypes.n_cols != paternity.n_elem) Rcpp::stop("offspring phenotypes must have column for each individual");
  if (maternal_phenotype.n_elem != 2) Rcpp::stop("maternal phenotype must have 2 elements");
  if (offspring_phenotypes.max() >= number_of_alleles) Rcpp::stop("offspring allele out of range");
  if (maternal_phenotype.max() >= number_of_alleles) Rcpp::stop("maternal allele out of range");
  if (arma::any(allele_frequencies_normalized < 0.)) Rcpp::stop("negative allele frequencies");
  if (dropout_rate <= 0. || mistyping_rate <= 0.) Rcpp::stop("negative genotyping error rates");

  double halfsib_likelihood = 0.;
  double halfsib_running_maximum = -arma::datum::inf;
  for (unsigned w=1; w<=number_of_alleles; ++w) // first maternal allele 
  { 
    for (unsigned v=w; v<=number_of_alleles; ++v) // second maternal allele
    {
      double maternal_genotype_probability = 
        (2.-int(w==v)) * allele_frequencies_normalized[w-1] * allele_frequencies_normalized[v-1]; //hwe prior
      double log_halfsib_likelihood = log(maternal_genotype_probability); 
      if (arma::prod(maternal_phenotype))
      {
        double maternal_phenotype_probability = 
          genotyping_error_model(maternal_phenotype, w, v, number_of_alleles, dropout_rate, mistyping_rate);
        log_halfsib_likelihood += log(maternal_phenotype_probability);
      }
      for (auto father : fathers)
      {
        double fullsib_likelihood = 0.;
        double fullsib_running_maximum = -arma::datum::inf;
        arma::uvec offspring_from_father = arma::find(paternity == father);
        for (unsigned u=1; u<=number_of_alleles; ++u) // paternal allele
        {
          double paternal_genotype_probability = 
            allele_frequencies_normalized[u-1]; //hwe prior
          double log_fullsib_likelihood = log(paternal_genotype_probability);
          for (auto offspring : offspring_from_father)
          {
            arma::uvec offspring_phenotype = offspring_phenotypes.col(offspring);
            if (arma::prod(offspring_phenotype)) { 
              double offspring_phenotype_probability = // Mendelian segregation probs * phenotype probabilities
                0.5 * genotyping_error_model(offspring_phenotype, w, u, number_of_alleles, dropout_rate, mistyping_rate) + 
                0.5 * genotyping_error_model(offspring_phenotype, v, u, number_of_alleles, dropout_rate, mistyping_rate); 
              log_fullsib_likelihood += log(offspring_phenotype_probability);
            } 
          }
          if (log_fullsib_likelihood <= fullsib_running_maximum) //underflow protection
          {
            fullsib_likelihood += exp(log_fullsib_likelihood - fullsib_running_maximum);
          } else {
            fullsib_likelihood *= exp(fullsib_running_maximum - log_fullsib_likelihood);
            fullsib_likelihood += 1.;
            fullsib_running_maximum = log_fullsib_likelihood;
          }
        }
        log_halfsib_likelihood += log(fullsib_likelihood) + fullsib_running_maximum;
      }
      if (log_halfsib_likelihood <= halfsib_running_maximum) //underflow protection
      {
        halfsib_likelihood += exp(log_halfsib_likelihood - halfsib_running_maximum);
      } else {
        halfsib_likelihood *= exp(halfsib_running_maximum - log_halfsib_likelihood);
        halfsib_likelihood += 1.;
        halfsib_running_maximum = log_halfsib_likelihood;
      }
    }
  }
  return log(halfsib_likelihood) + halfsib_running_maximum;
}

// [[Rcpp::export]]
arma::imat missing_data_problem (arma::imat input)
{
  input.t().print("input");
  arma::ivec input_unique = arma::unique(input);
  input_unique.t().print("unique");
  return input;
}

// [[Rcpp::export]]
double paternity_loglikelihood 
 (arma::uvec paternity, 
  arma::ucube offspring_phenotypes, 
  arma::umat maternal_phenotype,
  std::vector<arma::vec> allele_frequencies,
  arma::vec dropout_rate,
  arma::vec mistyping_rate)
{
  // check number of loci match
  const unsigned number_of_loci = allele_frequencies.size();
  if (maternal_phenotype.n_cols != number_of_loci) Rcpp::stop("must have maternal phenotypes for each locus");
  if (offspring_phenotypes.n_slices != number_of_loci) Rcpp::stop("must have offspring phenotypes for each locus");
  if (dropout_rate.n_elem != number_of_loci) Rcpp::stop("must have dropout rates for each locus");
  if (mistyping_rate.n_elem != number_of_loci) Rcpp::stop("must have mistyping rates for each locus");

  double log_likelihood = 0.;
  for (unsigned locus=0; locus<number_of_loci; ++locus)
  {
    log_likelihood += 
      paternity_loglikelihood_by_locus(paternity, offspring_phenotypes.slice(locus), maternal_phenotype.col(locus), allele_frequencies[locus], dropout_rate[locus], mistyping_rate[locus]);
  }
  return log_likelihood;
}

// [[Rcpp::export]]
Rcpp::List optimize_paternity_given_error_rates
 (arma::uvec paternity, 
  arma::ucube offspring_phenotypes, 
  arma::umat maternal_phenotype,
  std::vector<arma::vec> allele_frequencies,
  arma::vec dropout_rate,
  arma::vec mistyping_rate)
{
  const unsigned max_iter = 1000;
  const double convergence_tolerance = 1e-8;
  unsigned iter;
  double current_loglik = -arma::datum::inf;
  double old_loglik = -arma::datum::inf;
  for (iter=0; iter<=max_iter; ++iter)
  {
    paternity = recode_to_contiguous_integers(paternity);
    old_loglik = current_loglik;
    for (unsigned sib=0; sib<paternity.n_elem; ++sib)
    {
      unsigned current_number_of_fathers = paternity.max() + 1;
      arma::vec log_likelihood (current_number_of_fathers + 1);
      for (unsigned father=0; father<=current_number_of_fathers; ++father)
      {
        paternity[sib] = father;
        log_likelihood[father] = paternity_loglikelihood(paternity, offspring_phenotypes, maternal_phenotype, allele_frequencies, dropout_rate, mistyping_rate);
      }
      paternity[sib] = log_likelihood.index_max();
      paternity = recode_to_contiguous_integers(paternity);
      current_loglik = log_likelihood.max();
      //std::cout << iter << " " << sib << " " << log_likelihood.index_max() << " " << log_likelihood.max() << std::endl;
      //log_likelihood.t().print("loglik");
      //paternity.t().print("paternity");
    }
    Rcpp::Rcout << "[" << iter << "] " << "loglik: " << current_loglik << ", delta: " << current_loglik - old_loglik << std::endl;
    if (current_loglik - old_loglik < convergence_tolerance) break;
  }
  return Rcpp::List::create(
      Rcpp::_["paternity"] = paternity,
      Rcpp::_["loglikelihood"] = current_loglik,
      Rcpp::_["iterations"] = iter,
      Rcpp::_["converged"] = iter < max_iter
      );
}

std::vector<arma::vec> collapse_alleles_and_generate_genotype_prior
 (arma::ucube& phenotypes, 
  const bool add_unsampled_allele = true)
{
  // recode alleles into 1-based contiguous integers and return allele frequency prior
  
  const unsigned num_loci = phenotypes.n_slices;
  std::vector<arma::vec> allele_frequencies;

  for (unsigned locus=0; locus<num_loci; ++locus)
  {
    arma::uvec alleles = arma::nonzeros(phenotypes.slice(locus));
    for (int i=1; i<=alleles.n_elem; ++i)
    {
      phenotypes.slice(locus).replace(alleles.at(i-1), i);
    }
    unsigned number_of_alleles = alleles.n_elem + int(add_unsampled_allele);
    arma::vec frequencies (number_of_alleles); 
    frequencies.fill(1./double(number_of_alleles));
    allele_frequencies.push_back(frequencies); 
  }
  return allele_frequencies;
}

// [[Rcpp::export]]
Rcpp::List collapse_alleles_and_generate_prior_wrapper 
 (arma::ucube phenotypes, 
  const unsigned mother = 1,
  const bool add_unsampled_allele = true)
{
  //TODO right now missing data are zeros but that will not work b/c i use alleles in zero-based indexing
  //     i think i could make the genotype stuff 1-based. this would be loops for(v=; for(w=; for(u=; then change allele_frequencies[u]
  //EVERYTHING IS BROKEN -- NO I THINK I FIXED IT
  std::vector<arma::vec> out = collapse_alleles_and_generate_genotype_prior(phenotypes, add_unsampled_allele);
  arma::umat maternal_phenotype = phenotypes.tube(arma::span::all, arma::span(mother-1));
  arma::ucube offspring_phenotypes = phenotypes;
  offspring_phenotypes.shed_col(mother-1);
  return Rcpp::List::create(
      Rcpp::_["frequencies"] = out,
      Rcpp::_["offspring"] = offspring_phenotypes,
      Rcpp::_["maternal"] = maternal_phenotype);
}

// [[Rcpp::export]]
Rcpp::List sample_paternity_and_error_rates_from_joint_posterior
 (arma::ucube phenotypes, 
  const unsigned mother = 1,
  const unsigned number_of_mcmc_samples = 1000,
  const bool global_genotyping_error_rates = true)
{
  // samples from posterior distribution of full sib groups with Dirichlet process prior,
  // using algorithm 8 from Neal 2000 JCGS with m = 1
  
  if (mother > phenotypes.n_cols || mother < 1) Rcpp::stop("1-based index of mother out of range");

  // negative allele codes are treated as missing
  phenotypes.elem(arma::find(phenotypes == 0)).fill(arma::datum::nan);

  // TODO pass in phenotypes and index of mother, modify "collapse_alleles_and_generate_genotype_prior" to take single input
  const unsigned max_iter = number_of_mcmc_samples;
  const unsigned num_loci = phenotypes.n_slices;
  const unsigned num_offspring = phenotypes.n_cols - 1;

  // priors (hardcoded for now)
  const double alpha = 1.; //dirichlet process concentration parameter
  const arma::vec dropout_rate_prior = {{1.,1.}}; //beta(number of dropout homozygotes, number of heterozygotes)
  const arma::vec mistyping_rate_prior = {{1.,1.}}; //beta(number of mistypes, number of correct calls)
  std::vector<arma::vec> allele_frequencies = 
    collapse_alleles_and_generate_genotype_prior(phenotypes, true); //creates uniform frequency prior

  // split maternal, offspring phenotypes
  arma::umat maternal_phenotype = phenotypes.tube(arma::span::all, arma::span(mother-1));
  arma::ucube offspring_phenotypes = phenotypes; offspring_phenotypes.shed_col(mother-1);

  // initialize (could draw from prior instead)
  arma::uvec paternity = arma::ones<arma::uvec>(num_offspring);
  arma::vec dropout_rate (num_loci); dropout_rate.fill(0.05);
  arma::vec mistyping_rate (num_loci); mistyping_rate.fill(0.05);

  // storage
  arma::umat paternity_samples (num_offspring, max_iter);
  arma::mat dropout_rate_samples (num_loci, max_iter);
  arma::mat mistyping_rate_samples (num_loci, max_iter);
  arma::vec deviance_samples (max_iter);

  double deviance = 0.;
  paternity = recode_to_contiguous_integers(paternity);
  for (unsigned iter=0; iter<max_iter; ++iter)
  {
    // update paternity vector
    for (unsigned sib=0; sib<num_offspring; ++sib)
    {
      // tally size of sib groups
      unsigned current_number_of_fathers = paternity.max() + 1;
      arma::uvec offspring_per_father (current_number_of_fathers + 1, arma::fill::zeros);
      for (auto i : paternity) offspring_per_father[i]++;
      bool sib_is_not_singleton = offspring_per_father[paternity[sib]] > 1;
      offspring_per_father[paternity[sib]]--; 

      // conditional paternity probabilities
      arma::vec log_likelihood (current_number_of_fathers + unsigned(sib_is_not_singleton));
      for (unsigned father=0; father<log_likelihood.n_elem; ++father)
      {
        paternity[sib] = father;
        log_likelihood[father] = paternity_loglikelihood(paternity, offspring_phenotypes, maternal_phenotype, allele_frequencies, dropout_rate, mistyping_rate);
        double log_prior = offspring_per_father[father] > 0 ? log(double(offspring_per_father[father])) : log(alpha);
        log_likelihood[father] += -log(double(paternity.n_elem) - 1. + alpha) + log_prior;
      }

      // sample new father
      paternity[sib] = sample(arma::exp(log_likelihood - log_likelihood.max()));
      deviance = -2 * log_likelihood[paternity[sib]];
      paternity = recode_to_contiguous_integers(paternity); 
      //why recode? indices will increase, if pre-existing singleton is moved to a father with a higher index
    }

    // update error rates via data augmentation
    arma::umat global_error_counts = arma::zeros<arma::umat>(2,2);
    for (unsigned locus=0; locus<num_loci; ++locus)
    {
      arma::umat error_counts =
        sample_genotyping_errors_given_paternity(paternity, offspring_phenotypes.slice(locus), 
            maternal_phenotype.col(locus), allele_frequencies[locus], dropout_rate[locus], mistyping_rate[locus]);
      if (global_genotyping_error_rates)
      {
        global_error_counts += error_counts;
      } else {
        dropout_rate[locus] = 0.5 * R::rbeta(1. + error_counts.at(0,0), 1. + error_counts.at(0,1));
        mistyping_rate[locus] = R::rbeta(1. + error_counts.at(1,0), 1. + error_counts.at(1,1));
      }
    }
    if (global_genotyping_error_rates)
    {
        dropout_rate.fill(0.5 * R::rbeta(1. + global_error_counts.at(0,0), 1. + global_error_counts.at(0,1)));
        mistyping_rate.fill(R::rbeta(1. + global_error_counts.at(1,0), 1. + global_error_counts.at(1,1)));
    }

    // store state
    paternity_samples.col(iter) = paternity;
    dropout_rate_samples.col(iter) = dropout_rate;
    mistyping_rate_samples.col(iter) = mistyping_rate;
    deviance_samples.at(iter) = deviance;
    if (iter % 100 == 0) Rcpp::Rcout << "[" << iter << "] " << "deviance: " << deviance << std::endl;
  }
  return Rcpp::List::create(
    Rcpp::_["paternity"] = paternity_samples,
    Rcpp::_["dropout_rate"] = dropout_rate_samples,
    Rcpp::_["mistyping_rate"] = mistyping_rate_samples,
    Rcpp::_["deviance"] = deviance_samples);
}

// optimize_error_rates_given_paternity ==> do this in R with nelder mead or something

