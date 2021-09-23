#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <vector>

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]

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
  const double e2 = mistyping_rate;
  const double E1 = dropout_rate/(1 + dropout_rate);
  const double E2 = mistyping_rate/double(number_of_alleles-1);
  const arma::uvec genotype = {genotype0, genotype1};
  const bool phenotype_is_homozygous = phenotype[0] == phenotype[1];
  const bool genotype_is_homozygous = genotype[0] == genotype[1];
  if (number_of_alleles == 1) return 1.; //monomorphic loci
  if (genotype_is_homozygous) 
  {
    if (phenotype_is_homozygous && phenotype[0] == genotype[0])
    {
      return std::pow(1.-E2, 2);
    } else if ((phenotype[0] == genotype[0] && phenotype[1] != genotype[0]) || (phenotype[0] != genotype[0] && phenotype[1] == genotype[0])) {
      return 2.*e2*(1-E2);
    } else if (phenotype[0] != genotype[0] && phenotype[1] != genotype[0]) {
      return (2.-int(phenotype_is_homozygous))*std::pow(e2, 2);
    }  
  } else {
    if (phenotype[0] == genotype[0] && phenotype[1] == genotype[1])
    {
      return std::pow(1.-E2, 2) + std::pow(e2, 2) + 2.*e1*std::pow(1.-E2-e2, 2);
    } else if (phenotype_is_homozygous && (phenotype[0] == genotype[0] || phenotype[0] == genotype[1])) {
      return e2*(1.-E2) + e1*std::pow(1.-E2-e2, 2);
    } else if (phenotype[0] != genotype[0] && phenotype[0] != genotype[1] && phenotype[1] != genotype[0] && phenotype[1] != genotype[1]) {
      return (2.-int(phenotype_is_homozygous))*std::pow(e2, 2);
    } else {
      return e2*(1.-E2+e2);
    }
  }
  return 0.;
}

// [[Rcpp::export]]
double paternity_loglikelihood_by_locus 
 (const arma::uvec& paternity,
  const arma::umat& offspring_phenotypes, 
  const arma::uvec& maternal_phenotype, 
  const arma::vec& allele_frequencies, 
  const double& dropout_rate, 
  const double& mistyping_rate)
{
  // likelihood of offspring paternity give offspring phenotypes, maternal phenotype, haplodiploidy
  // modified from Eqs 3 & 4 in Wang 2004 Genetics
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
  double halfsib_likelihood = 0.;
  for (unsigned w=0; w<number_of_alleles; ++w) // first maternal allele 
  { 
    for (unsigned v=w; v<number_of_alleles; ++v) // second maternal allele
    {
      double maternal_genotype_probability = 
        (2.-int(w==v)) * allele_frequencies_normalized[w] * allele_frequencies_normalized[v];
      double log_halfsib_likelihood = log(maternal_genotype_probability); 
      if (maternal_phenotype.is_finite())
      {
        double maternal_phenotype_probability = 
          genotyping_error_model(maternal_phenotype, w, v, number_of_alleles, dropout_rate, mistyping_rate);
        log_halfsib_likelihood += log(maternal_phenotype_probability);
      }
      for (auto father : fathers)
      {
        double fullsib_likelihood = 0.;
        arma::uvec offspring_from_father = arma::find(paternity == father);
        for (unsigned u=0; u<number_of_alleles; ++u) // paternal allele
        {
          double paternal_genotype_probability = allele_frequencies_normalized[u];
          double log_fullsib_likelihood = log(paternal_genotype_probability);
          for (auto offspring : offspring_from_father)
          {
            arma::uvec offspring_phenotype = offspring_phenotypes.col(offspring);
            if (offspring_phenotype.is_finite()) { 
              double offspring_phenotype_probability = // Mendelian segregation probs * phenotype probabilities
                0.5 * genotyping_error_model(offspring_phenotype, w, u, number_of_alleles, dropout_rate, mistyping_rate) + 
                0.5 * genotyping_error_model(offspring_phenotype, v, u, number_of_alleles, dropout_rate, mistyping_rate); 
              log_fullsib_likelihood += log(offspring_phenotype_probability);
            }
          }
          fullsib_likelihood += exp(log_fullsib_likelihood);
        }
        log_halfsib_likelihood += log(fullsib_likelihood);
      }
      halfsib_likelihood += exp(log_halfsib_likelihood);
    }
  }
  return log(halfsib_likelihood);
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
arma::uvec recode_to_contiguous_integers (arma::uvec input)
{
  arma::uvec uniq = arma::unique(input);
  for(unsigned i=0; i<uniq.n_elem; ++i) input.replace(uniq.at(i), i);
  return input;
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

// [[Rcpp::export]]
Rcpp::List sample_paternity_given_error_rates
 (arma::uvec paternity, 
  arma::ucube offspring_phenotypes, 
  arma::umat maternal_phenotype,
  std::vector<arma::vec> allele_frequencies,
  arma::vec dropout_rate,
  arma::vec mistyping_rate,
  const unsigned max_iter)
{
  // samples from posterior distribution with Dirichlet process prior using algorithm 8 from Neal 2000 JCGS with m = 1
  const unsigned max_fathers = paternity.n_elem;
  const double alpha = 1.;
  double deviance;
  arma::umat paternity_samples (paternity.n_elem, max_iter+1);
  arma::vec deviance_samples (max_iter+1);
  paternity = recode_to_contiguous_integers(paternity);
  for (unsigned iter=0; iter<=max_iter; ++iter)
  {
    for (unsigned sib=0; sib<paternity.n_elem; ++sib)
    {
      // tally size of sib groups
      unsigned current_number_of_fathers = paternity.max() + 1;
      arma::uvec offspring_per_father (current_number_of_fathers + 1, arma::fill::zeros);
      for (auto i : paternity) offspring_per_father[i]++;
      bool sib_is_not_singleton = offspring_per_father[paternity[sib]] > 1;
      offspring_per_father[paternity[sib]]--; //if was singleton, we should have 0 before the "buffer zero" at the end, and will never hit buffer
      // prior and likelihood
      arma::vec log_likelihood (current_number_of_fathers + unsigned(sib_is_not_singleton));
      for (unsigned father=0; father<log_likelihood.n_elem; ++father)
      {
        paternity[sib] = father;
        log_likelihood[father] = paternity_loglikelihood(paternity, offspring_phenotypes, maternal_phenotype, allele_frequencies, dropout_rate, mistyping_rate);
        double log_prior = offspring_per_father[father] > 0 ? log(double(offspring_per_father[father])) : log(alpha);
        log_likelihood[father] += -log(double(paternity.n_elem) - 1. + alpha) + log_prior;
      }
      // sample new paternity
      Rcpp::IntegerVector possible_fathers (log_likelihood.n_elem);
      Rcpp::NumericVector paternity_probability (log_likelihood.n_elem);
      for (unsigned father=0; father<log_likelihood.n_elem; ++father)
      {
        possible_fathers[father] = father;
        paternity_probability[father] = exp(log_likelihood[father]);
      }
      Rcpp::IntegerVector draw = Rcpp::sample(possible_fathers, 1, false, paternity_probability);
      paternity[sib] = draw(0);
      deviance = -2 * log_likelihood[draw(0)];
      paternity = recode_to_contiguous_integers(paternity); //why? without this it crawls up if preexsiting singleton is rejected
    }
    // store state
    paternity_samples.col(iter) = paternity;
    deviance_samples.at(iter) = deviance;
    //Rcpp::Rcout << "[" << iter << "] " << "loglik: " << deviance << std::endl;
  }
  return Rcpp::List::create(
      Rcpp::_["paternity"] = paternity_samples,
      Rcpp::_["deviance"] = deviance_samples
      );
}

// optimize_error_rates_given_paternity
