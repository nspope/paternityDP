#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <vector>
#include <tuple>

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
double log_ascending_factorial (const double x, const unsigned r)
{
  return ::Rf_lgammafn(x + double(r)) - ::Rf_lgammafn(x);
}

// [[Rcpp::export]]
double log_descending_factorial (const double x, const unsigned r)
{
  return ::Rf_lgammafn(x + 1.) - ::Rf_lgammafn(x - double(r) + 1.);
}

// [[Rcpp::export]]
double log_uniform_MFM_prior (const unsigned n, const unsigned t, const double gamma, const unsigned max_number_of_components)
{
  // MFM coefficients assuming uniform prior; from Eq 3.2 in Miller & Harrison 2015
  // see https://github.com/jwmi/BayesianMixtures.jl/blob/master/src/MFM.jl for efficient implementation
  double out = 0.;
  double running_maximum = -arma::datum::inf;
  for (unsigned k=t; k<=max_number_of_components; ++k)
  {
    double incr = log_descending_factorial(double(k), t) -
      log_ascending_factorial(gamma * double(k), n) - 
      log(double(max_number_of_components));
    if (incr <= running_maximum) //underflow protection
    {
      out += exp(incr - running_maximum);
    } else {
      out *= exp(running_maximum - incr);
      out += 1.;
      running_maximum = incr;
    }
  }
  return log(out) + running_maximum;
}

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

std::vector<arma::vec> collapse_alleles_and_generate_genotype_prior
 (arma::ucube& phenotypes, 
  const bool add_unsampled_allele = false)
{
  // recode alleles into 1-based contiguous integers and return allele frequency prior
  
  const unsigned num_loci = phenotypes.n_slices;
  std::vector<arma::vec> allele_frequencies;

  for (unsigned locus=0; locus<num_loci; ++locus)
  {
    arma::uvec alleles = arma::unique(arma::nonzeros(phenotypes.slice(locus)));
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

std::vector<arma::uvec> unique_alleles
 (const arma::ucube& phenotypes, 
  const bool add_unsampled_allele = false)
{
  const unsigned num_loci = phenotypes.n_slices;
  std::vector<arma::uvec> allele_names;
  for (unsigned locus=0; locus<num_loci; ++locus)
  {
    arma::uvec alleles = arma::unique(arma::nonzeros(phenotypes.slice(locus)));
    if (add_unsampled_allele) alleles = arma::join_vert(alleles, arma::uvec({999}));
    allele_names.push_back(alleles);
  }
  return allele_names;
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

// ---------------------------------------------------------------------------- //

std::vector<arma::uvec> sample_genotyping_errors_and_allele_counts_given_paternity
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
  if (offspring_phenotypes.max() > number_of_alleles) Rcpp::stop("offspring allele out of range");
  if (maternal_phenotype.max() > number_of_alleles) Rcpp::stop("maternal allele out of range");
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
  for (auto father : fathers)
  {
    arma::vec paternal_genotype_posterior (number_of_alleles);
    arma::uvec possible_paternal_genotypes (number_of_alleles);
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
      possible_paternal_genotypes[u-1] = u;
    }
    paternal_genotype_posterior -= paternal_genotype_posterior.max();
    paternal_genotypes.at(father) = possible_paternal_genotypes.at(sample(arma::exp(paternal_genotype_posterior)));
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
  //TODO simplify redundant tallying, no need for "counts_of_errors"
  unsigned sampled_phenotypes = 0;
  unsigned sampled_heterozygotes = 0;
  arma::uvec counts_of_errors = {0,0}; //dropouts, mistypes
  arma::uvec dropout_errors (paternity.n_elem+1, arma::fill::zeros);
  arma::uvec mistype_errors (paternity.n_elem+1, arma::fill::zeros);
  if (arma::prod(maternal_phenotype))
  {
    sampled_phenotypes++;
    if (maternal_genotype[0] != maternal_genotype[1]) sampled_heterozygotes++;
    arma::uvec errors =
      simulate_genotyping_errors(maternal_phenotype, maternal_genotype[0], 
          maternal_genotype[1], number_of_alleles, dropout_rate, mistyping_rate);
    dropout_errors.at(0) = errors.at(0);
    mistype_errors.at(0) = errors.at(1);
    counts_of_errors += errors;
  }
  for (unsigned sib=0; sib<paternity.n_elem; ++sib)
  {
    arma::uvec offspring_phenotype = offspring_phenotypes.col(sib);
    if (arma::prod(offspring_phenotype)) { 
      sampled_phenotypes++;
      if (offspring_genotypes.at(0,sib) != offspring_genotypes.at(1,sib)) sampled_heterozygotes++;
      arma::uvec errors =
        simulate_genotyping_errors(offspring_phenotype, offspring_genotypes.at(0,sib), 
          offspring_genotypes.at(1,sib), number_of_alleles, dropout_rate, mistyping_rate);
      dropout_errors.at(sib+1) = errors.at(0);
      mistype_errors.at(sib+1) = errors.at(1);
      counts_of_errors += errors;
    }
  }
  arma::uvec dropout_counts = {counts_of_errors[0], sampled_heterozygotes-counts_of_errors[0]};
  arma::uvec mistype_counts = {counts_of_errors[1], 2*sampled_phenotypes-counts_of_errors[1]};

  // tally maternal/paternal alleles
  arma::uvec allele_counts (number_of_alleles, arma::fill::zeros);
  allele_counts[maternal_genotype[0]-1]++;
  allele_counts[maternal_genotype[1]-1]++;
  for (unsigned i=0; i<paternal_genotypes.n_elem; ++i) allele_counts[paternal_genotypes[i]-1]++;

  // output
  std::vector<arma::uvec> counts;
  counts.push_back(dropout_counts);
  counts.push_back(mistype_counts);
  counts.push_back(allele_counts);
  counts.push_back(dropout_errors);
  counts.push_back(mistype_errors);
  return counts;
}

// [[Rcpp::export]]
Rcpp::List sample_error_rates_given_paternity
 (arma::ucube phenotypes, 
  arma::uvec paternity, 
  const unsigned mother = 1,
  const unsigned number_of_mcmc_samples = 1000,
  const unsigned global_genotyping_error_rates = false,
  const bool random_allele_frequencies = true,
  const bool add_unsampled_allele = true)
{
  if (mother > phenotypes.n_cols || mother < 1) Rcpp::stop("1-based index of mother out of range");

  const unsigned max_iter = number_of_mcmc_samples;
  const unsigned number_of_loci = phenotypes.n_slices;
  const unsigned num_offspring = phenotypes.n_cols - 1;

  // priors (hardcoded for now)
  const arma::vec dropout_rate_prior = {{1.,1.}}; //beta(number of dropout homozygotes, number of heterozygotes)
  const arma::vec mistyping_rate_prior = {{1.,1.}}; //beta(number of mistypes, number of correct calls)
  std::vector<arma::vec> allele_frequencies = 
    collapse_alleles_and_generate_genotype_prior(phenotypes, add_unsampled_allele); //creates uniform frequency prior

  // split maternal, offspring phenotypes
  arma::umat maternal_phenotype = phenotypes.tube(arma::span::all, arma::span(mother-1));
  arma::ucube offspring_phenotypes = phenotypes; offspring_phenotypes.shed_col(mother-1);

  // initialize (could draw from prior instead)
  arma::vec dropout_rate (number_of_loci); dropout_rate.fill(0.05);
  arma::vec mistyping_rate (number_of_loci); mistyping_rate.fill(0.05);
  
  arma::mat dropout_rate_samples (number_of_loci, max_iter);
  arma::mat dropout_errors (paternity.n_elem+1, number_of_loci, arma::fill::zeros);
  arma::mat mistyping_rate_samples (number_of_loci, max_iter);
  arma::mat mistyping_errors (paternity.n_elem+1, number_of_loci, arma::fill::zeros);
  paternity = recode_to_contiguous_integers(paternity);
  for (unsigned iter=0; iter<max_iter; ++iter)
  {
    arma::uvec global_dropout_counts (2, arma::fill::zeros);
    arma::uvec global_mistype_counts (2, arma::fill::zeros);
    for (unsigned locus=0; locus<number_of_loci; ++locus)
    {
      std::vector<arma::uvec> error_counts =
        sample_genotyping_errors_and_allele_counts_given_paternity(paternity, offspring_phenotypes.slice(locus), 
            maternal_phenotype.col(locus), allele_frequencies[locus], dropout_rate[locus], mistyping_rate[locus]);

      // track number of errors
      for (unsigned i=0; i<paternity.n_elem+1; ++i)
      {
        dropout_errors.at(i,locus) += double(error_counts[3][i]);
        mistyping_errors.at(i,locus) += double(error_counts[4][i]);
      }

      // update error rates
      global_dropout_counts += error_counts[0]; global_mistype_counts += error_counts[1];
      dropout_rate[locus] = 0.5 * R::rbeta(1. + error_counts[0][0], 1. + error_counts[0][1]);
      mistyping_rate[locus] = R::rbeta(1. + error_counts[1][0], 1. + error_counts[1][1]);

      // update allele frequencies
      if (random_allele_frequencies)
      {
        for (unsigned allele=0; allele<allele_frequencies[locus].n_elem; ++allele)
        {
          allele_frequencies[locus][allele] = R::rgamma(1. + error_counts[2][allele], 1.);
        }
        allele_frequencies[locus] /= arma::accu(allele_frequencies[locus]);
      }
    }
    if (global_genotyping_error_rates)
    {
      // overwrite per-locus rates with global rate
      dropout_rate.fill(0.5 * R::rbeta(1. + global_dropout_counts[0], 1. + global_dropout_counts[1]));
      mistyping_rate.fill(R::rbeta(1. + global_mistype_counts[0], 1. + global_mistype_counts[1]));
    }
    
    dropout_rate_samples.col(iter) = dropout_rate;
    mistyping_rate_samples.col(iter) = mistyping_rate;
    if (iter % 100 == 0) Rcpp::Rcout << "[" << iter << "]" << std::endl;
  }

  // posterior expectation of error counts & reshuffle matrices to mirror input
  dropout_errors /= double(max_iter);
  arma::rowvec maternal_dropout_errors = dropout_errors.row(0); 
  dropout_errors.shed_row(0); dropout_errors.insert_rows(mother-1, maternal_dropout_errors);
  mistyping_errors /= double(max_iter);
  arma::rowvec maternal_mistyping_errors = mistyping_errors.row(0); 
  mistyping_errors.shed_row(0); mistyping_errors.insert_rows(mother-1, maternal_mistyping_errors);

  return Rcpp::List::create(
      Rcpp::_["dropout_rate"] = dropout_rate_samples,
      Rcpp::_["mistyping_rate"] = mistyping_rate_samples,
      Rcpp::_["dropout_errors"] = arma::trans(dropout_errors),
      Rcpp::_["mistyping_errors"] = arma::trans(mistyping_errors)
      );
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
  if (offspring_phenotypes.max() > number_of_alleles) Rcpp::stop("offspring allele out of range");
  if (maternal_phenotype.max() > number_of_alleles) Rcpp::stop("maternal allele out of range");
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
 (arma::ucube phenotypes,
  arma::vec dropout_rate,
  arma::vec mistyping_rate,
  const unsigned mother = 1)
{
  if (mother > phenotypes.n_cols || mother < 1) Rcpp::stop("1-based index of mother out of range");

  const unsigned num_offspring = phenotypes.n_cols - 1;

  std::vector<arma::vec> allele_frequencies = 
    collapse_alleles_and_generate_genotype_prior(phenotypes); //creates uniform frequency prior
  arma::umat maternal_phenotype = phenotypes.tube(arma::span::all, arma::span(mother-1));
  arma::ucube offspring_phenotypes = phenotypes; offspring_phenotypes.shed_col(mother-1);
  arma::uvec paternity = arma::ones<arma::uvec>(num_offspring);

  const unsigned max_iter = 1000;
  const double convergence_tolerance = 1e-8;

  unsigned iter;
  double current_loglik = -arma::datum::inf;
  double old_loglik = -arma::datum::inf;
  for (iter=0; iter<=max_iter; ++iter)
  {
    paternity = recode_to_contiguous_integers(paternity);
    old_loglik = current_loglik;
    for (unsigned sib=0; sib<num_offspring; ++sib)
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
arma::mat loglikelihood_of_error_rates_given_paternity
 (arma::ucube phenotypes,
  arma::uvec paternity,
  arma::mat grid_of_error_rates,
  const unsigned mother = 1)
{
  if (mother > phenotypes.n_cols || mother < 1) Rcpp::stop("1-based index of mother out of range");

  const unsigned num_offspring = phenotypes.n_cols - 1;
  const unsigned number_of_loci = phenotypes.n_slices;

  std::vector<arma::vec> allele_frequencies = 
    collapse_alleles_and_generate_genotype_prior(phenotypes); //creates uniform frequency prior
  arma::umat maternal_phenotype = phenotypes.tube(arma::span::all, arma::span(mother-1));
  arma::ucube offspring_phenotypes = phenotypes; offspring_phenotypes.shed_col(mother-1);

  arma::mat log_likelihood (grid_of_error_rates.n_rows, number_of_loci);
  paternity = recode_to_contiguous_integers(paternity);
  for (unsigned i=0; i<grid_of_error_rates.n_rows; ++i)
  {
    for (unsigned j=0; j<number_of_loci; ++j)
    {
      log_likelihood.at(i,j) = 
        paternity_loglikelihood_by_locus(paternity, offspring_phenotypes.slice(j), maternal_phenotype.col(j), 
            allele_frequencies[j], grid_of_error_rates.at(i,0), grid_of_error_rates.at(i,1));
    }
  }
  return log_likelihood;
}

// [[Rcpp::export]]
Rcpp::List collapse_alleles_and_generate_prior_wrapper 
 (arma::ucube phenotypes, 
  const unsigned mother = 1,
  const bool add_unsampled_allele = false)
{
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
  const bool global_genotyping_error_rates = true,
  const double concentration = 1.)
{
  // samples from posterior distribution of full sib groups with Dirichlet process prior,
  // using algorithm 8 from Neal 2000 JCGS with m = 1
  
  if (mother > phenotypes.n_cols || mother < 1) Rcpp::stop("1-based index of mother out of range");

  const unsigned max_iter = number_of_mcmc_samples;
  const unsigned num_loci = phenotypes.n_slices;
  const unsigned num_offspring = phenotypes.n_cols - 1;

  // priors (hardcoded for now)
  const double alpha = std::fabs(concentration); //dirichlet process concentration parameter
  const double gamma = std::fabs(concentration); //mfm concentration parameter
  const double delta = 1.; //allele frequency concentration parameter
  const arma::vec dropout_rate_prior = {{1.,1.}}; //beta(number of dropout homozygotes, number of heterozygotes)
  const arma::vec mistyping_rate_prior = {{1.,1.}}; //beta(number of mistypes, number of correct calls)
  std::vector<arma::vec> allele_frequencies = collapse_alleles_and_generate_genotype_prior(phenotypes, true);

  // calculate coefficients needed for the MFM prior
  const unsigned max_number_of_fathers = num_offspring;
  arma::vec log_mfm_prior (num_offspring, arma::fill::zeros);
  if (gamma > 0.)
  {
    for (unsigned i=1; i<=num_offspring; ++i)
    {
      log_mfm_prior[i-1] = log_uniform_MFM_prior(num_offspring, i, gamma, max_number_of_fathers);
    }
  }
  if (arma::any(log_mfm_prior > 0.)) Rcpp::stop("problem with prior, this should not have happened");

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
  arma::uvec number_of_fathers_samples (max_iter);
  arma::mat dropout_errors (num_offspring+1, num_loci, arma::fill::zeros);
  arma::mat mistyping_errors (num_offspring+1, num_loci, arma::fill::zeros);

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

        // "restraunt process" prior
        double log_prior = 0.;
        if (concentration > 0.) // Dirichlet process prior
        {
          log_prior = offspring_per_father[father] > 0 ? 
            log(double(offspring_per_father[father])) - log(double(paternity.n_elem)-1.+alpha): 
            log(alpha) - log(double(paternity.n_elem)-1.+alpha);
        } 
        if (concentration < 0.) // Mixture-of-finite-mixture prior
        { 
          unsigned t = current_number_of_fathers - unsigned(!sib_is_not_singleton);
          log_prior = offspring_per_father[father] > 0 ? 
            log(double(offspring_per_father[father]) + gamma) : 
            log(gamma) + log_mfm_prior[t+1-1] - log_mfm_prior[t-1];
        }
        log_likelihood[father] += log_prior;
      }

      // sample new father
      paternity[sib] = sample(arma::exp(log_likelihood - log_likelihood.max()));
      deviance = -2 * log_likelihood[paternity[sib]];
      paternity = recode_to_contiguous_integers(paternity); 
      //why recode? indices will increase, if pre-existing singleton is moved to a father with a higher index
    }

    // update error rates and allele frequencies via data augmentation
    arma::uvec global_dropout_counts (2, arma::fill::zeros);
    arma::uvec global_mistype_counts (2, arma::fill::zeros);
    for (unsigned locus=0; locus<num_loci; ++locus)
    {
      std::vector<arma::uvec> error_counts =
        sample_genotyping_errors_and_allele_counts_given_paternity(paternity, offspring_phenotypes.slice(locus), 
            maternal_phenotype.col(locus), allele_frequencies[locus], dropout_rate[locus], mistyping_rate[locus]);

      // track number of errors
      for (unsigned i=0; i<paternity.n_elem+1; ++i)
      {
        dropout_errors.at(i,locus) += double(error_counts[3][i]);
        mistyping_errors.at(i,locus) += double(error_counts[4][i]);
      }

      // update error rates
      global_dropout_counts += error_counts[0]; global_mistype_counts += error_counts[1];
      dropout_rate[locus] = 0.5 * R::rbeta(1. + error_counts[0][0], 1. + error_counts[0][1]);
      mistyping_rate[locus] = R::rbeta(1. + error_counts[1][0], 1. + error_counts[1][1]);

      // update allele frequencies
      for (unsigned allele=0; allele<allele_frequencies[locus].n_elem; ++allele)
      {
        allele_frequencies[locus][allele] = R::rgamma(1. + error_counts[2][allele], 1.);
      }
      allele_frequencies[locus] /= arma::accu(allele_frequencies[locus]);
    }
    if (global_genotyping_error_rates)
    {
      // overwrite per-locus rates with global rate
      dropout_rate.fill(0.5 * R::rbeta(1. + global_dropout_counts[0], 1. + global_dropout_counts[1]));
      mistyping_rate.fill(R::rbeta(1. + global_mistype_counts[0], 1. + global_mistype_counts[1]));
    }

    // store state
    paternity_samples.col(iter) = paternity;
    dropout_rate_samples.col(iter) = dropout_rate;
    mistyping_rate_samples.col(iter) = mistyping_rate;
    deviance_samples.at(iter) = deviance;
    number_of_fathers_samples.at(iter) = paternity.max() + 1;

    if (iter % 100 == 0) Rcpp::Rcout << "[" << iter << "] " << "deviance: " << deviance << std::endl; //print # fathers too
  }

  // posterior expectation of error counts & reorder matrices to mirror input
  dropout_errors /= double(max_iter);
  arma::rowvec maternal_dropout_errors = dropout_errors.row(0); 
  dropout_errors.shed_row(0); dropout_errors.insert_rows(mother-1, maternal_dropout_errors);
  mistyping_errors /= double(max_iter);
  arma::rowvec maternal_mistyping_errors = mistyping_errors.row(0); 
  mistyping_errors.shed_row(0); mistyping_errors.insert_rows(mother-1, maternal_mistyping_errors);

  return Rcpp::List::create(
    Rcpp::_["paternity"] = paternity_samples,
    Rcpp::_["dropout_rate"] = dropout_rate_samples,
    Rcpp::_["mistyping_rate"] = mistyping_rate_samples,
    Rcpp::_["number_of_fathers"] = number_of_fathers_samples,
    Rcpp::_["dropout_errors"] = dropout_errors,
    Rcpp::_["mistyping_errors"] = mistyping_errors,
    Rcpp::_["deviance"] = deviance_samples);
}

// ------------------------- alternative implementation without DP ------------------------ //

// [[Rcpp::export]]
arma::uvec sample_matrix (arma::mat probabilities)
{
  arma::uvec draw (2);
  probabilities /= arma::accu(probabilities);
  arma::vec row_sums = arma::sum(probabilities, 1);
  draw[0] = sample(row_sums);
  draw[1] = sample(arma::trans(probabilities.row(draw[0])));
  return draw;
}

// [[Rcpp::export]]
arma::ucube select_columns_from_cube (arma::ucube input, arma::uvec which)
{
  if (which.max() >= input.n_cols) Rcpp::stop("Out of column range");
  arma::ucube output = input;
  arma::uvec drop = arma::regspace<arma::uvec>(0, input.n_cols-1);
  drop.shed_rows(which);
  for (unsigned i=drop.n_elem; i>0; --i) output.shed_col(drop[i-1]);
  return output;
}

// [[Rcpp::export]]
double phenotype_error_model
 (const arma::uvec& phenotype,
  const arma::uvec& genotype,
  const unsigned& number_of_alleles,
  const double& dropout_rate, 
  const double& mistyping_rate)
{
  // from Eqs 1 & 2 in Wang 2004 Genetics

  if (phenotype.n_elem != genotype.n_elem) Rcpp::stop("Genotype and phenotype are different lengths");
  if (phenotype.n_elem > 2) Rcpp::stop("Ploidy cannot be greater than two");
  if (number_of_alleles == 1) return 1.; //monomorphic loci

  const double e1 = dropout_rate;
  const double e2 = mistyping_rate/double(number_of_alleles-1);
  const double E2 = mistyping_rate;

  const bool genotype_is_diploid = phenotype.n_elem == 2;

  if (genotype_is_diploid)
  {
    const bool phenotype_is_homozygous = phenotype[0] == phenotype[1];
    const bool genotype_is_homozygous = genotype[0] == genotype[1];
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
  } else { // haploid, only mistyping errors
    if (phenotype[0] == genotype[0])
    {
      return 1.-E2; // 1 - Pr(other allele) - Pr(other other allele) - ... = 1 - e2 - e2 - .... = 1 - prob of error
    } else {
      return e2; // Pr(other allele) = e2
    }
  }
  return 0.;
}

// [[Rcpp::export]]
arma::uvec sample_phenotype_errors
 (const arma::uvec& phenotype,
  const arma::uvec& genotype,
  const unsigned& number_of_alleles,
  const double& dropout_rate, 
  const double& mistyping_rate)
{
  if (phenotype.n_elem != genotype.n_elem) Rcpp::stop("Genotype and phenotype are different lengths");
  if (phenotype.n_elem > 2) Rcpp::stop("Ploidy cannot be greater than two");
  if (number_of_alleles == 1) return arma::zeros<arma::uvec>(phenotype.n_elem); // monomorphic loci

  const double e1 = dropout_rate;
  const double e2 = mistyping_rate/double(number_of_alleles-1);
  const double E2 = mistyping_rate;

  const bool genotype_is_diploid = phenotype.n_elem == 2;

  if (genotype_is_diploid)
  {
    const bool phenotype_is_homozygous = phenotype[0] == phenotype[1];
    const bool genotype_is_homozygous = genotype[0] == genotype[1];
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
  } else {
    if (phenotype[0] == genotype[0])
    {
      return arma::uvec({0,0});
    } else {
      return arma::uvec({0,1});
    }
  }
  return arma::zeros<arma::uvec>(2);
}

// [[Rcpp::export]]
double mendelian_genotype_model
 (const arma::uvec& offspring_phenotype,
  const arma::uvec& maternal_genotype,
  const arma::uvec& paternal_genotype,
  const unsigned& number_of_alleles,
  const double& dropout_rate, 
  const double& mistyping_rate)
{
  const bool offspring_is_haploid = offspring_phenotype.n_elem == 1;
  arma::uvec offspring_genotype (offspring_phenotype.n_elem);
  double likelihood = 0.;
  unsigned number_of_possible_genotypes = 0;
  for (auto maternal_allele : maternal_genotype)
  {
    offspring_genotype[0] = maternal_allele;
    if (offspring_is_haploid)
    {
      likelihood += phenotype_error_model(offspring_phenotype, offspring_genotype, number_of_alleles, dropout_rate, mistyping_rate);
      number_of_possible_genotypes++;
    } else {
      for (auto paternal_allele : paternal_genotype)
      {
        offspring_genotype[1] = paternal_allele;
        likelihood += phenotype_error_model(offspring_phenotype, offspring_genotype, number_of_alleles, dropout_rate, mistyping_rate);
        number_of_possible_genotypes++;
      }
    }
  }
  likelihood /= double(number_of_possible_genotypes);
  return likelihood;
}

// [[Rcpp::export]]
arma::uvec sample_mendelian_genotype
 (const arma::uvec& offspring_phenotype,
  const arma::uvec& maternal_genotype,
  const arma::uvec& paternal_genotype,
  const unsigned& number_of_alleles,
  const double& dropout_rate, 
  const double& mistyping_rate)
{
  // this is poorly written, clean it up
  const bool offspring_is_haploid = offspring_phenotype.n_elem == 1;
  arma::uvec offspring_genotype (offspring_phenotype.n_elem);
  arma::mat likelihood = offspring_is_haploid ? 
    arma::zeros(maternal_genotype.n_elem, 1) : 
    arma::zeros(maternal_genotype.n_elem, paternal_genotype.n_elem);
  for (unsigned maternal_allele=0; maternal_allele<maternal_genotype.n_elem; ++maternal_allele)
  {
    offspring_genotype[0] = maternal_genotype[maternal_allele];
    if (offspring_is_haploid)
    {
      likelihood.at(maternal_allele, 1) = 
        phenotype_error_model(offspring_phenotype, offspring_genotype, number_of_alleles, dropout_rate, mistyping_rate);
    } else {
      for (unsigned paternal_allele=0; paternal_allele<paternal_genotype.n_elem; ++paternal_allele)
      {
        offspring_genotype[1] = paternal_genotype[paternal_allele];
        likelihood.at(maternal_allele, paternal_allele) =
          phenotype_error_model(offspring_phenotype, offspring_genotype, number_of_alleles, dropout_rate, mistyping_rate);
      }
    }
  }
  arma::uvec new_phenotype = sample_matrix(likelihood).head(offspring_genotype.n_elem);
  new_phenotype[0] = maternal_genotype[new_phenotype[0]];
  if (!offspring_is_haploid)
  {
    new_phenotype[1] = paternal_genotype[new_phenotype[1]];
  }
  return new_phenotype;
}

// [[Rcpp::export]]
Rcpp::List sample_parentage_and_error_rates_from_joint_posterior
 (arma::ucube phenotypes, 
  arma::uvec mothers,
  arma::uvec fathers,
  arma::uvec holdouts,
  const unsigned number_of_mcmc_samples = 1000,
  const unsigned burn_in_samples = 100,
  const unsigned thinning_interval = 1,
  const bool global_genotyping_error_rates = true)
{
  if (mothers.n_elem == 0 || fathers.n_elem == 0) Rcpp::stop("Must have at least one potential father and mother");
  if (mothers.max() > phenotypes.n_cols || mothers.min() < 1) Rcpp::stop("1-based index of mothers out of range");
  if (fathers.max() > phenotypes.n_cols || fathers.min() < 1) Rcpp::stop("1-based index of fathers out of range");
  if (holdouts.max() > phenotypes.n_cols || holdouts.min() < 1) Rcpp::stop("1-based index of holdouts out of range");

  //this way of specifying parental phenotypes could allow monoiecious mating systems in the future
  mothers = arma::unique(mothers);
  fathers = arma::unique(fathers);
  holdouts = arma::unique(holdouts);
  arma::uvec parents_and_holdouts = arma::unique(arma::join_vert(arma::join_vert(mothers, fathers), holdouts));
  arma::uvec offspring = arma::regspace<arma::uvec>(1, phenotypes.n_cols);
  offspring.shed_rows(parents_and_holdouts-1);
  if (parents_and_holdouts.n_elem != holdouts.n_elem + mothers.n_elem + fathers.n_elem) Rcpp::stop("Overlap in mother/father/holdouts");
  if (phenotypes.n_cols != offspring.n_elem + parents_and_holdouts.n_elem) Rcpp::stop("Splitting issues?");

  // recode alleles to integers, split parental and offspring phenotypes
  std::vector<arma::vec> allele_frequencies = collapse_alleles_and_generate_genotype_prior(phenotypes, false);
  mothers -= 1; fathers -= 1; holdouts -= 1; parents_and_holdouts -= 1; offspring -= 1; // convert to 0-based indices
  arma::ucube maternal_phenotypes = select_columns_from_cube(phenotypes, mothers); 
  arma::ucube paternal_phenotypes = select_columns_from_cube(phenotypes, fathers); 
  arma::ucube holdout_phenotypes = select_columns_from_cube(phenotypes, holdouts); 
  arma::ucube offspring_phenotypes = select_columns_from_cube(phenotypes, offspring); 
  paternal_phenotypes.shed_row(1); //diploid->haploid

  // priors (hardcoded for now)
  const double delta = 1.; //allele frequency concentration parameter
  const arma::vec dropout_rate_prior = {{1.,1.}}; //beta(number of dropout homozygotes, number of heterozygotes)
  const arma::vec mistyping_rate_prior = {{1.,1.}}; //beta(number of mistypes, number of correct calls)

  // dimensions
  const unsigned max_iter = number_of_mcmc_samples;
  const unsigned num_loci = phenotypes.n_slices;
  const unsigned num_mothers = mothers.n_elem;
  const unsigned num_fathers = fathers.n_elem;
  const unsigned num_holdouts = holdouts.n_elem;
  const unsigned num_offspring = offspring.n_elem;
  const unsigned num_samples = phenotypes.n_cols;
  arma::uvec num_alleles (num_loci); 
  for (unsigned locus=0; locus<num_loci; ++locus) num_alleles[locus] = allele_frequencies[locus].n_elem;

  // initialize (could draw from prior instead)
  arma::uvec maternity = arma::zeros<arma::uvec>(num_offspring);
  arma::uvec paternity = arma::zeros<arma::uvec>(num_offspring);
  arma::vec dropout_rate (num_loci); dropout_rate.fill(0.05);
  arma::vec mistyping_rate (num_loci); mistyping_rate.fill(0.05);
  arma::ucube maternal_genotypes = maternal_phenotypes; maternal_genotypes.replace(0, 1);
  arma::ucube paternal_genotypes = paternal_phenotypes; paternal_genotypes.replace(0, 1);
  arma::ucube offspring_genotypes = offspring_phenotypes; offspring_genotypes.replace(0, 1);
  // sample missing (0) from HWE prior TODO use MAP estimate for reproducible

  // storage
  arma::imat paternity_samples (num_samples, max_iter); paternity_samples.fill(arma::datum::nan);
  arma::imat maternity_samples (num_samples, max_iter); maternity_samples.fill(arma::datum::nan);
  arma::mat dropout_rate_samples (num_loci, max_iter); dropout_rate_samples.fill(arma::datum::nan);
  arma::mat mistyping_rate_samples (num_loci, max_iter); mistyping_rate_samples.fill(arma::datum::nan);
  arma::mat deviance_samples (num_samples, max_iter); deviance_samples.fill(arma::datum::nan);
  arma::mat holdout_deviance_samples (num_samples, max_iter); holdout_deviance_samples.fill(arma::datum::nan);
  std::vector<arma::mat> allele_frequencies_samples;
  for (unsigned locus=0; locus<num_loci; ++locus) allele_frequencies_samples.push_back(arma::zeros<arma::mat>(num_alleles[locus], max_iter));
  arma::mat dropout_error_expectation (num_samples, num_loci, arma::fill::zeros);
  arma::mat mistyping_error_expectation (num_samples, num_loci, arma::fill::zeros);
  std::vector<arma::cube> genotypes_expectation;
  for (unsigned locus=0; locus<num_loci; ++locus) 
  {
    genotypes_expectation.push_back(arma::zeros<arma::cube>(num_alleles[locus], num_alleles[locus], num_samples));
  }

  // Gibbs sampler
  int start = -int(burn_in_samples);
  for (int iter=start; iter<int(max_iter); ++iter)
  {
    for (unsigned thin=0; thin<thinning_interval; ++thin)
    {
      std::vector<arma::uvec> allele_counts;
      for (unsigned locus=0; locus<num_loci; ++locus) allele_counts.push_back(arma::zeros<arma::uvec>(num_alleles[locus]));
      arma::uvec num_phenotyped_alleles (num_loci, arma::fill::zeros); 
      arma::uvec num_heterozygotes (num_loci, arma::fill::zeros); //not counting haploids (fathers)
      arma::umat maternal_dropout_errors (num_mothers, num_loci, arma::fill::zeros);
      arma::umat paternal_dropout_errors (num_fathers, num_loci, arma::fill::zeros);
      arma::umat offspring_dropout_errors (num_offspring, num_loci, arma::fill::zeros);
      arma::umat maternal_mistyping_errors (num_mothers, num_loci, arma::fill::zeros);
      arma::umat paternal_mistyping_errors (num_fathers, num_loci, arma::fill::zeros);
      arma::umat offspring_mistyping_errors (num_offspring, num_loci, arma::fill::zeros);
      arma::vec deviance (num_offspring, arma::fill::zeros);
      arma::vec holdout_deviance (num_holdouts, arma::fill::zeros);

      // ------ sample maternal genotypes ------
      for (unsigned mother=0; mother<num_mothers; ++mother)
      {
        arma::uvec children = arma::find(maternity == mother);
        for (unsigned locus=0; locus<num_loci; ++locus)
        {
          arma::mat maternal_genotype_probabilities(num_alleles[locus], num_alleles[locus]);
          maternal_genotype_probabilities.fill(-arma::datum::inf);
          for (unsigned first_allele=0; first_allele<num_alleles[locus]; first_allele++)
          {
            for (unsigned second_allele=first_allele; second_allele<num_alleles[locus]; second_allele++)
            {
              // hardy-weinberg prior
              arma::uvec maternal_genotype ({first_allele+1, second_allele+1}); // 1-based allele indexing
              maternal_genotype_probabilities.at(first_allele, second_allele) = 
                log(2. - int(first_allele == second_allele)) +
                log(allele_frequencies[locus].at(first_allele)) + 
                log(allele_frequencies[locus].at(second_allele)); 

              // phenotype likelihoods
              if (maternal_phenotypes.at(0, mother, locus))
              {
                maternal_genotype_probabilities.at(first_allele, second_allele) += 
                  log(phenotype_error_model(maternal_phenotypes.slice(locus).col(mother), maternal_genotype,
                        num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]));
              }
              for (auto sib : children)
              {
                if (offspring_phenotypes.at(0, sib, locus))
                {
                  maternal_genotype_probabilities.at(first_allele, second_allele) += 
                    log(mendelian_genotype_model(offspring_phenotypes.slice(locus).col(sib), maternal_genotype,
                          paternal_genotypes.slice(locus).col(paternity[sib]), 
                          num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]));
                }
              }
            }
          }
          arma::uvec new_alleles = sample_matrix(arma::exp(maternal_genotype_probabilities - maternal_genotype_probabilities.max()));
          for (unsigned i=0; i<2; ++i)
          {
            allele_counts[locus].at(new_alleles[i])++;
            maternal_genotypes.at(i, mother, locus) = new_alleles[i] + 1; //1-based allele indexing
          }
        }
      }

      // ------ sample paternal genotypes ------ 
      for (unsigned father=0; father<num_fathers; ++father)
      {
        arma::uvec children = arma::find(paternity == father);
        for (unsigned locus=0; locus<num_loci; ++locus)
        {
          arma::vec paternal_genotype_probabilities(num_alleles[locus]);
          paternal_genotype_probabilities.fill(-arma::datum::inf);
          for (unsigned first_allele=0; first_allele<num_alleles[locus]; first_allele++)
          {
            // hardy-weinberg prior
            arma::uvec paternal_genotype ({first_allele + 1}); // 1-based allele indexing
            paternal_genotype_probabilities.at(first_allele) = log(allele_frequencies[locus].at(first_allele)); 

            // phenotype likelihoods
            if (paternal_phenotypes.at(0, father, locus))
            {
              paternal_genotype_probabilities.at(first_allele) += 
                log(phenotype_error_model(paternal_phenotypes.slice(locus).col(father), paternal_genotype,
                      num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]));
            }
            for (auto sib : children)
            {
              if (offspring_phenotypes.at(0, sib, locus))
              {
                paternal_genotype_probabilities.at(first_allele) += 
                  log(mendelian_genotype_model(offspring_phenotypes.slice(locus).col(sib), 
                        maternal_genotypes.slice(locus).col(maternity[sib]), paternal_genotype,
                        num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]));
              }
            }
          }
          unsigned new_allele = sample(arma::exp(paternal_genotype_probabilities - paternal_genotype_probabilities.max()));
          allele_counts[locus].at(new_allele)++;
          paternal_genotypes.at(0, father, locus) = new_allele + 1; // 1-based allele indexing
        }
      }

      // ------ sample allele frequencies ------
      for (unsigned locus=0; locus<num_loci; ++locus)
      {
        for (unsigned allele=0; allele<num_alleles[locus]; ++allele)
        {
          allele_frequencies[locus][allele] = R::rgamma(1. + allele_counts[locus][allele], 1.);
        }
        allele_frequencies[locus] /= arma::accu(allele_frequencies[locus]);
      }

      // ------ sample offspring genotypes ------
      for (unsigned sib=0; sib<num_offspring; ++sib)
      {
        for (unsigned locus=0; locus<num_loci; ++locus)
        {
          offspring_genotypes.slice(locus).col(sib) = 
            sample_mendelian_genotype(offspring_phenotypes.slice(locus).col(sib), 
                maternal_genotypes.slice(locus).col(maternity[sib]), 
                paternal_genotypes.slice(locus).col(paternity[sib]), 
                num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]);
        }
      }

      // ------ sample errors and error rates ------ 
      for (unsigned locus=0; locus<num_loci; ++locus)
      {
        // maternal phenotyping errors
        for (unsigned mother=0; mother<num_mothers; ++mother)
        {
          if (maternal_phenotypes.at(0,mother,locus))
          {
            num_phenotyped_alleles.at(locus) += 2;
            num_heterozygotes.at(locus) += int(maternal_genotypes.at(0,mother,locus) != maternal_genotypes.at(1,mother,locus));
            arma::uvec phenotyping_errors = 
              sample_phenotype_errors(maternal_phenotypes.slice(locus).col(mother),
                                      maternal_genotypes.slice(locus).col(mother),
                                      num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]);
            maternal_dropout_errors.at(mother, locus) += phenotyping_errors[0];
            maternal_mistyping_errors.at(mother, locus) += phenotyping_errors[1];
          }
        }

        // paternal phenotyping errors
        for (unsigned father=0; father<num_fathers; ++father)
        {
          if (paternal_phenotypes.at(0,father,locus))
          {
            num_phenotyped_alleles.at(locus) += 1;
            arma::uvec phenotyping_errors = 
              sample_phenotype_errors(paternal_phenotypes.slice(locus).col(father),
                                      paternal_genotypes.slice(locus).col(father),
                                      num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]);
            paternal_dropout_errors.at(father, locus) += phenotyping_errors[0];
            paternal_mistyping_errors.at(father, locus) += phenotyping_errors[1];
          }
        }

        // offspring phenotyping errors
        for (unsigned sib=0; sib<num_offspring; ++sib)
        {
          if (offspring_phenotypes.at(0,sib,locus))
          {
            num_phenotyped_alleles.at(locus) += 2;
            num_heterozygotes.at(locus) += int(offspring_genotypes.at(0,sib,locus) != offspring_genotypes.at(1,sib,locus));
            arma::uvec phenotyping_errors = 
              sample_phenotype_errors(offspring_phenotypes.slice(locus).col(sib),
                                      offspring_genotypes.slice(locus).col(sib),
                                      num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]);
            offspring_dropout_errors.at(sib, locus) += phenotyping_errors[0];
            offspring_mistyping_errors.at(sib, locus) += phenotyping_errors[1];
          }
        }

        // update error rates
        unsigned dropout_errors = arma::accu(offspring_dropout_errors.col(locus)) +
          arma::accu(maternal_dropout_errors.col(locus)) + arma::accu(paternal_dropout_errors.col(locus));
        unsigned mistyping_errors = arma::accu(offspring_mistyping_errors.col(locus)) +
          arma::accu(maternal_mistyping_errors.col(locus)) + arma::accu(paternal_mistyping_errors.col(locus));
        dropout_rate[locus] = 0.5 * R::rbeta(1. + dropout_errors, 1. + num_heterozygotes[locus] - dropout_errors);
        mistyping_rate[locus] = R::rbeta(1. + mistyping_errors, 1. + num_phenotyped_alleles[locus] - mistyping_errors);
      }
      if (global_genotyping_error_rates) // overwrite per-locus rates with global rate
      {
        unsigned dropout_errors = arma::accu(offspring_dropout_errors) +
          arma::accu(maternal_dropout_errors) + arma::accu(paternal_dropout_errors);
        unsigned mistyping_errors = arma::accu(offspring_mistyping_errors) +
          arma::accu(maternal_mistyping_errors) + arma::accu(paternal_mistyping_errors);
        dropout_rate.fill(0.5 * R::rbeta(1. + dropout_errors, 1. + arma::accu(num_heterozygotes) - dropout_errors));
        mistyping_rate.fill(R::rbeta(1. + mistyping_errors, 1. + arma::accu(num_phenotyped_alleles) - mistyping_errors));
      }

      // ------ sample parentage ------
      for (unsigned sib=0; sib<num_offspring; ++sib)
      {
        arma::mat log_likelihood (num_fathers, num_mothers, arma::fill::zeros);
        for (unsigned father=0; father<num_fathers; ++father)
        {
          for (unsigned mother=0; mother<num_mothers; ++mother)
          {
            for (unsigned locus=0; locus<num_loci; ++locus)
            {
              if (offspring_phenotypes.at(0, sib, locus))
              {
                log_likelihood.at(father, mother) += 
                      log(mendelian_genotype_model(offspring_phenotypes.slice(locus).col(sib), 
                                                   maternal_genotypes.slice(locus).col(mother),
                                                   paternal_genotypes.slice(locus).col(father),
                                                   num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]));
              }
            }
          }
        }
        arma::uvec new_parentage = sample_matrix(arma::exp(log_likelihood - log_likelihood.max()));
        paternity[sib] = new_parentage[0];
        maternity[sib] = new_parentage[1];
        deviance.at(sib) = log_likelihood.at(paternity[sib], maternity[sib]);
      }

      if (thin == 0 && iter >= 0)
      {
        // calculate posterior predictive on holdout set
        for (unsigned extra=0; extra<num_holdouts; ++extra)
        {
          arma::mat posterior_predictive (num_fathers, num_mothers, arma::fill::zeros);
          for (unsigned father=0; father<num_fathers; ++father)
          {
            for (unsigned mother=0; mother<num_mothers; ++mother)
            {
              for (unsigned locus=0; locus<num_loci; ++locus)
              {
                if (holdout_phenotypes.at(0, extra, locus))
                {
                  posterior_predictive.at(father, mother) += 
                    log(mendelian_genotype_model(holdout_phenotypes.slice(locus).col(extra), 
                          maternal_genotypes.slice(locus).col(mother),
                          paternal_genotypes.slice(locus).col(father),
                          num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]));
                }
              }
            }
          }
          holdout_deviance.at(extra) = 
            log(arma::accu(arma::exp(posterior_predictive - posterior_predictive.max()))) + 
            posterior_predictive.max() - log(posterior_predictive.n_elem);
        }

        // store state, mapping samples back to original indicies
        dropout_rate_samples.col(iter) = dropout_rate;
        mistyping_rate_samples.col(iter) = mistyping_rate;
        for (unsigned locus=0; locus<num_loci; ++locus)
        {
          allele_frequencies_samples[locus].col(iter) = allele_frequencies[locus];
          for (unsigned mother=0; mother<num_mothers; ++mother)
          {
            dropout_error_expectation.at(mothers[mother], locus) += maternal_dropout_errors.at(mother, locus)/double(max_iter);
            mistyping_error_expectation.at(mothers[mother], locus) += maternal_mistyping_errors.at(mother, locus)/double(max_iter);
            genotypes_expectation[locus].at(maternal_genotypes.at(0, mother, locus)-1, 
                maternal_genotypes.at(1, mother, locus)-1, mothers[mother]) += 1./double(max_iter);
          }
          for (unsigned father=0; father<num_fathers; ++father)
          {
            dropout_error_expectation.at(fathers[father], locus) += paternal_dropout_errors.at(father, locus)/double(max_iter);
            mistyping_error_expectation.at(fathers[father], locus) += paternal_mistyping_errors.at(father, locus)/double(max_iter);
            genotypes_expectation[locus].at(paternal_genotypes.at(0, father, locus)-1, 
                paternal_genotypes.at(0, father, locus)-1, fathers[father]) += 1./double(max_iter);
          }
          for (unsigned sib=0; sib<num_offspring; ++sib)
          {
            dropout_error_expectation.at(offspring[sib], locus) += offspring_dropout_errors.at(sib, locus)/double(max_iter);
            mistyping_error_expectation.at(offspring[sib], locus) += offspring_mistyping_errors.at(sib, locus)/double(max_iter);
            genotypes_expectation[locus].at(offspring_genotypes.at(0, sib, locus)-1, 
                offspring_genotypes.at(1, sib, locus)-1, offspring[sib]) += 1./double(max_iter);
          }
        }
        for (unsigned sib=0; sib<num_offspring; ++sib)
        {
          paternity_samples.at(offspring[sib], iter) = paternity[sib];
          maternity_samples.at(offspring[sib], iter) = maternity[sib];
          deviance_samples.at(offspring[sib], iter) = deviance.at(sib);
        }
        for (unsigned extra=0; extra<num_holdouts; ++extra)
        {
          holdout_deviance_samples.at(holdouts[extra], iter) = holdout_deviance.at(extra);
        }

        // progress
        if (iter % 100 == 0) Rcpp::Rcout << "[" << iter << "] " << "deviance: " << arma::accu(deviance) << std::endl;
      }
    }
  }

  // Rcpp collapses dimensions for elements in std::vector, so copy these over to Rcpp::List
  Rcpp::List genotypes_expectation_wrapped = Rcpp::List::create();
  Rcpp::List allele_frequencies_samples_wrapped = Rcpp::List::create();
  for (unsigned locus=0; locus<num_loci; ++locus) {
    genotypes_expectation_wrapped.push_back(genotypes_expectation[locus]);
    allele_frequencies_samples_wrapped.push_back(allele_frequencies_samples[locus]);
  }

  return Rcpp::List::create(
    Rcpp::_["paternity"] = paternity_samples,
    Rcpp::_["maternity"] = maternity_samples,
    Rcpp::_["dropout_rate"] = dropout_rate_samples,
    Rcpp::_["mistyping_rate"] = mistyping_rate_samples,
    Rcpp::_["dropout_errors"] = dropout_error_expectation,
    Rcpp::_["mistyping_errors"] = mistyping_error_expectation,
    Rcpp::_["genotypes"] = genotypes_expectation_wrapped,
    Rcpp::_["allele_frequencies"] = allele_frequencies_samples_wrapped,
    Rcpp::_["holdout_deviance"] = holdout_deviance_samples,
    Rcpp::_["deviance"] = deviance_samples);
}

// ----------- stupid attempt to use DP prior (still gets stuck/doesnt converge) ------------- //

// [[Rcpp::export]]
Rcpp::List sample_parentage_and_error_rates_from_joint_posterior_alt
 (arma::ucube phenotypes, 
  arma::uvec mothers,
  arma::uvec fathers,
  const double concentration = 1.,
  const unsigned number_of_mcmc_samples = 1000,
  const unsigned burn_in_samples = 100,
  const unsigned thinning_interval = 1,
  const bool global_genotyping_error_rates = true,
  const bool sample_from_prior = false,
  const bool random_initialization = false)
{ // uses DP prior on mothers fathers, so pad with lots of missing values in phenotypes
  if (mothers.n_elem == 0 || fathers.n_elem == 0) Rcpp::stop("Must have at least one potential father and mother");
  if (mothers.max() > phenotypes.n_cols || mothers.min() < 1) Rcpp::stop("1-based index of mothers out of range");
  if (fathers.max() > phenotypes.n_cols || fathers.min() < 1) Rcpp::stop("1-based index of fathers out of range");

  //this way of specifying parental phenotypes could allow monoiecious mating systems in the future
  mothers = arma::unique(mothers);
  fathers = arma::unique(fathers);
  arma::uvec parents = arma::unique(arma::join_vert(mothers, fathers));
  arma::uvec offspring = arma::regspace<arma::uvec>(1, phenotypes.n_cols);
  offspring.shed_rows(parents-1);
  if (parents.n_elem != mothers.n_elem + fathers.n_elem) Rcpp::stop("Overlap in mother/father");
  if (phenotypes.n_cols != offspring.n_elem + parents.n_elem) Rcpp::stop("Splitting issues?");

  // recode alleles to integers, split parental and offspring phenotypes
  std::vector<arma::uvec> allele_lengths = unique_alleles(phenotypes, false);
  std::vector<arma::vec> allele_frequencies = collapse_alleles_and_generate_genotype_prior(phenotypes, false);
  mothers -= 1; fathers -= 1; parents -= 1; offspring -= 1; // convert to 0-based indices
  arma::ucube maternal_phenotypes = select_columns_from_cube(phenotypes, mothers); 
  arma::ucube paternal_phenotypes = select_columns_from_cube(phenotypes, fathers); 
  arma::ucube offspring_phenotypes = select_columns_from_cube(phenotypes, offspring); 
  paternal_phenotypes.shed_row(1); //diploid->haploid

  // if sampling from prior, set all phenotype data to missing
  if (sample_from_prior)
  {
    maternal_phenotypes.fill(0);
    paternal_phenotypes.fill(0);
    offspring_phenotypes.fill(0);
  }

  // priors (hardcoded for now)
  //const double concentration = 1.; //dirichlet process concentration parameter
  const double delta = 1.; //allele frequency concentration parameter
  const arma::vec dropout_rate_prior = {{1.,1.}}; //beta(number of dropout homozygotes, number of heterozygotes)
  const arma::vec mistyping_rate_prior = {{1.,1.}}; //beta(number of mistypes, number of correct calls)

  // dimensions
  const unsigned max_iter = number_of_mcmc_samples;
  const unsigned num_loci = phenotypes.n_slices;
  const unsigned num_mothers = mothers.n_elem;
  const unsigned num_fathers = fathers.n_elem;
  const unsigned num_offspring = offspring.n_elem;
  const unsigned num_samples = phenotypes.n_cols;
  arma::uvec num_alleles (num_loci); 
  for (unsigned locus=0; locus<num_loci; ++locus) num_alleles[locus] = allele_frequencies[locus].n_elem;

  // initialize (could draw from prior instead)
  arma::uvec maternity = random_initialization ? 
    arma::randi<arma::uvec>(num_offspring, arma::distr_param(0, num_mothers-1)) :
    arma::zeros<arma::uvec>(num_offspring);
  arma::uvec paternity = random_initialization ? 
    arma::randi<arma::uvec>(num_offspring, arma::distr_param(0, num_fathers-1)) :
    arma::zeros<arma::uvec>(num_offspring);
  arma::vec dropout_rate (num_loci); dropout_rate.fill(0.05);
  arma::vec mistyping_rate (num_loci); mistyping_rate.fill(0.05);
  arma::ucube maternal_genotypes = maternal_phenotypes; maternal_genotypes.replace(0, 1);
  arma::ucube paternal_genotypes = paternal_phenotypes; paternal_genotypes.replace(0, 1);
  arma::ucube offspring_genotypes = offspring_phenotypes; offspring_genotypes.replace(0, 1);
  // sample missing (0) from HWE prior TODO use MAP estimate for reproducible

  // storage
  arma::imat paternity_samples (num_samples, max_iter); paternity_samples.fill(arma::datum::nan);
  arma::imat maternity_samples (num_samples, max_iter); maternity_samples.fill(arma::datum::nan);
  arma::mat dropout_rate_samples (num_loci, max_iter); dropout_rate_samples.fill(arma::datum::nan);
  arma::mat mistyping_rate_samples (num_loci, max_iter); mistyping_rate_samples.fill(arma::datum::nan);
  arma::mat deviance_samples (num_samples, max_iter); deviance_samples.fill(arma::datum::nan);
  arma::mat holdout_deviance_samples (num_samples, max_iter); holdout_deviance_samples.fill(arma::datum::nan);
  std::vector<arma::mat> allele_frequencies_samples;
  for (unsigned locus=0; locus<num_loci; ++locus) allele_frequencies_samples.push_back(arma::zeros<arma::mat>(num_alleles[locus], max_iter));
  arma::mat dropout_error_expectation (num_samples, num_loci, arma::fill::zeros);
  arma::mat mistyping_error_expectation (num_samples, num_loci, arma::fill::zeros);
  std::vector<arma::cube> genotypes_expectation;
  for (unsigned locus=0; locus<num_loci; ++locus) 
  {
    genotypes_expectation.push_back(arma::zeros<arma::cube>(num_alleles[locus], num_alleles[locus], num_samples));
  }

  // Gibbs sampler
  int start = -int(burn_in_samples);
  for (int iter=start; iter<int(max_iter); ++iter)
  {
    for (unsigned thin=0; thin<thinning_interval; ++thin)
    {
      std::vector<arma::uvec> allele_counts;
      for (unsigned locus=0; locus<num_loci; ++locus) allele_counts.push_back(arma::zeros<arma::uvec>(num_alleles[locus]));
      arma::uvec num_phenotyped_alleles (num_loci, arma::fill::zeros); 
      arma::uvec num_heterozygotes (num_loci, arma::fill::zeros); //not counting haploids (fathers)
      arma::umat maternal_dropout_errors (num_mothers, num_loci, arma::fill::zeros);
      arma::umat paternal_dropout_errors (num_fathers, num_loci, arma::fill::zeros);
      arma::umat offspring_dropout_errors (num_offspring, num_loci, arma::fill::zeros);
      arma::umat maternal_mistyping_errors (num_mothers, num_loci, arma::fill::zeros);
      arma::umat paternal_mistyping_errors (num_fathers, num_loci, arma::fill::zeros);
      arma::umat offspring_mistyping_errors (num_offspring, num_loci, arma::fill::zeros);
      arma::vec deviance (num_offspring, arma::fill::zeros);

      // ------ sample maternal genotypes ------
      for (unsigned mother=0; mother<num_mothers; ++mother)
      {
        arma::uvec children = arma::find(maternity == mother);
        for (unsigned locus=0; locus<num_loci; ++locus)
        {
          arma::mat maternal_genotype_probabilities(num_alleles[locus], num_alleles[locus]);
          maternal_genotype_probabilities.fill(-arma::datum::inf);
          for (unsigned first_allele=0; first_allele<num_alleles[locus]; first_allele++)
          {
            for (unsigned second_allele=first_allele; second_allele<num_alleles[locus]; second_allele++)
            {
              // hardy-weinberg prior
              arma::uvec maternal_genotype ({first_allele+1, second_allele+1}); // 1-based allele indexing
              maternal_genotype_probabilities.at(first_allele, second_allele) = 
                log(2. - int(first_allele == second_allele)) +
                log(allele_frequencies[locus].at(first_allele)) + 
                log(allele_frequencies[locus].at(second_allele)); 

              // phenotype likelihoods
              if (maternal_phenotypes.at(0, mother, locus))
              {
                maternal_genotype_probabilities.at(first_allele, second_allele) += 
                  log(phenotype_error_model(maternal_phenotypes.slice(locus).col(mother), maternal_genotype,
                        num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]));
              }
              for (auto sib : children)
              {
                if (offspring_phenotypes.at(0, sib, locus))
                {
                  maternal_genotype_probabilities.at(first_allele, second_allele) += 
                    log(mendelian_genotype_model(offspring_phenotypes.slice(locus).col(sib), maternal_genotype,
                          paternal_genotypes.slice(locus).col(paternity[sib]), 
                          num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]));
                }
              }
            }
          }
          arma::uvec new_alleles = sample_matrix(arma::exp(maternal_genotype_probabilities - maternal_genotype_probabilities.max()));
          for (unsigned i=0; i<2; ++i)
          {
            allele_counts[locus].at(new_alleles[i])++;
            maternal_genotypes.at(i, mother, locus) = new_alleles[i] + 1; //1-based allele indexing
          }
        }
      }

      // ------ sample paternal genotypes ------ 
      for (unsigned father=0; father<num_fathers; ++father)
      {
        arma::uvec children = arma::find(paternity == father);
        for (unsigned locus=0; locus<num_loci; ++locus)
        {
          arma::vec paternal_genotype_probabilities(num_alleles[locus]);
          paternal_genotype_probabilities.fill(-arma::datum::inf);
          for (unsigned first_allele=0; first_allele<num_alleles[locus]; first_allele++)
          {
            // hardy-weinberg prior
            arma::uvec paternal_genotype ({first_allele + 1}); // 1-based allele indexing
            paternal_genotype_probabilities.at(first_allele) = log(allele_frequencies[locus].at(first_allele)); 

            // phenotype likelihoods
            if (paternal_phenotypes.at(0, father, locus))
            {
              paternal_genotype_probabilities.at(first_allele) += 
                log(phenotype_error_model(paternal_phenotypes.slice(locus).col(father), paternal_genotype,
                      num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]));
            }
            for (auto sib : children)
            {
              if (offspring_phenotypes.at(0, sib, locus))
              {
                paternal_genotype_probabilities.at(first_allele) += 
                  log(mendelian_genotype_model(offspring_phenotypes.slice(locus).col(sib), 
                        maternal_genotypes.slice(locus).col(maternity[sib]), paternal_genotype,
                        num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]));
              }
            }
          }
          unsigned new_allele = sample(arma::exp(paternal_genotype_probabilities - paternal_genotype_probabilities.max()));
          allele_counts[locus].at(new_allele)++;
          paternal_genotypes.at(0, father, locus) = new_allele + 1; // 1-based allele indexing
        }
      }

      // ------ sample allele frequencies ------
      for (unsigned locus=0; locus<num_loci; ++locus)
      {
        for (unsigned allele=0; allele<num_alleles[locus]; ++allele)
        {
          allele_frequencies[locus][allele] = R::rgamma(1. + allele_counts[locus][allele], 1.);
        }
        allele_frequencies[locus] /= arma::accu(allele_frequencies[locus]);
      }

      // ------ sample offspring genotypes ------
      for (unsigned sib=0; sib<num_offspring; ++sib)
      {
        for (unsigned locus=0; locus<num_loci; ++locus)
        {
          offspring_genotypes.slice(locus).col(sib) = 
            sample_mendelian_genotype(offspring_phenotypes.slice(locus).col(sib), 
                maternal_genotypes.slice(locus).col(maternity[sib]), 
                paternal_genotypes.slice(locus).col(paternity[sib]), 
                num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]);
        }
      }

      // ------ sample errors and error rates ------ 
      for (unsigned locus=0; locus<num_loci; ++locus)
      {
        // maternal phenotyping errors
        for (unsigned mother=0; mother<num_mothers; ++mother)
        {
          if (maternal_phenotypes.at(0,mother,locus))
          {
            num_phenotyped_alleles.at(locus) += 2;
            num_heterozygotes.at(locus) += int(maternal_genotypes.at(0,mother,locus) != maternal_genotypes.at(1,mother,locus));
            arma::uvec phenotyping_errors = 
              sample_phenotype_errors(maternal_phenotypes.slice(locus).col(mother),
                                      maternal_genotypes.slice(locus).col(mother),
                                      num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]);
            maternal_dropout_errors.at(mother, locus) += phenotyping_errors[0];
            maternal_mistyping_errors.at(mother, locus) += phenotyping_errors[1];
          }
        }

        // paternal phenotyping errors
        for (unsigned father=0; father<num_fathers; ++father)
        {
          if (paternal_phenotypes.at(0,father,locus))
          {
            num_phenotyped_alleles.at(locus) += 1;
            arma::uvec phenotyping_errors = 
              sample_phenotype_errors(paternal_phenotypes.slice(locus).col(father),
                                      paternal_genotypes.slice(locus).col(father),
                                      num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]);
            paternal_dropout_errors.at(father, locus) += phenotyping_errors[0];
            paternal_mistyping_errors.at(father, locus) += phenotyping_errors[1];
          }
        }

        // offspring phenotyping errors
        for (unsigned sib=0; sib<num_offspring; ++sib)
        {
          if (offspring_phenotypes.at(0,sib,locus))
          {
            num_phenotyped_alleles.at(locus) += 2;
            num_heterozygotes.at(locus) += int(offspring_genotypes.at(0,sib,locus) != offspring_genotypes.at(1,sib,locus));
            arma::uvec phenotyping_errors = 
              sample_phenotype_errors(offspring_phenotypes.slice(locus).col(sib),
                                      offspring_genotypes.slice(locus).col(sib),
                                      num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]);
            offspring_dropout_errors.at(sib, locus) += phenotyping_errors[0];
            offspring_mistyping_errors.at(sib, locus) += phenotyping_errors[1];
          }
        }

        // update error rates
        unsigned dropout_errors = arma::accu(offspring_dropout_errors.col(locus)) +
          arma::accu(maternal_dropout_errors.col(locus)) + arma::accu(paternal_dropout_errors.col(locus));
        unsigned mistyping_errors = arma::accu(offspring_mistyping_errors.col(locus)) +
          arma::accu(maternal_mistyping_errors.col(locus)) + arma::accu(paternal_mistyping_errors.col(locus));
        dropout_rate[locus] = 0.5 * R::rbeta(1. + dropout_errors, 1. + num_heterozygotes[locus] - dropout_errors);
        mistyping_rate[locus] = R::rbeta(1. + mistyping_errors, 1. + num_phenotyped_alleles[locus] - mistyping_errors);
      }
      if (global_genotyping_error_rates) // overwrite per-locus rates with global rate
      {
        unsigned dropout_errors = arma::accu(offspring_dropout_errors) +
          arma::accu(maternal_dropout_errors) + arma::accu(paternal_dropout_errors);
        unsigned mistyping_errors = arma::accu(offspring_mistyping_errors) +
          arma::accu(maternal_mistyping_errors) + arma::accu(paternal_mistyping_errors);
        dropout_rate.fill(0.5 * R::rbeta(1. + dropout_errors, 1. + arma::accu(num_heterozygotes) - dropout_errors));
        mistyping_rate.fill(R::rbeta(1. + mistyping_errors, 1. + arma::accu(num_phenotyped_alleles) - mistyping_errors));
      }

      // ------ sample parentage ------
      arma::uvec paternity_counts(num_fathers, arma::fill::zeros);
      arma::uvec maternity_counts(num_mothers, arma::fill::zeros);
      for (unsigned sib=0; sib<num_offspring; ++sib)
      {
        paternity_counts[paternity[sib]]++; maternity_counts[maternity[sib]]++;
      }
      for (unsigned sib=0; sib<num_offspring; ++sib)
      {
        paternity_counts[paternity[sib]]--; maternity_counts[maternity[sib]]--;
        arma::uvec realized_fathers = arma::find(paternity_counts > 0);
        arma::uvec realized_mothers = arma::find(maternity_counts > 0);
        unsigned auxiliary_fathers = num_fathers - realized_fathers.n_elem;
        unsigned auxiliary_mothers = num_mothers - realized_mothers.n_elem;
        arma::mat log_prior (num_fathers, num_mothers, arma::fill::zeros);
        arma::mat log_likelihood (num_fathers, num_mothers, arma::fill::zeros);
        for (unsigned father=0; father<num_fathers; ++father)
        {
          for (unsigned mother=0; mother<num_mothers; ++mother)
          {
            for (unsigned locus=0; locus<num_loci; ++locus)
            {
              if (offspring_phenotypes.at(0, sib, locus))
              {
                log_likelihood.at(father, mother) += 
                      log(mendelian_genotype_model(offspring_phenotypes.slice(locus).col(sib), 
                                                   maternal_genotypes.slice(locus).col(mother),
                                                   paternal_genotypes.slice(locus).col(father),
                                                   num_alleles[locus], dropout_rate[locus], mistyping_rate[locus]));
              }
            }
            log_prior.at(father, mother) += paternity_counts[father] > 0 ? 
              log(double(paternity_counts[father]))-log(double(num_offspring)-1.+concentration): 
              log(concentration)-log(double(auxiliary_fathers))-log(double(num_offspring)-1.+concentration);
            log_prior.at(father, mother) += maternity_counts[mother] > 0 ? 
              log(double(maternity_counts[mother]))-log(double(num_offspring)-1.+concentration): 
              log(concentration)-log(double(auxiliary_mothers))-log(double(num_offspring)-1.+concentration);
          }
        }
        arma::mat log_posterior = log_likelihood + log_prior;
        arma::uvec new_parentage = sample_matrix(arma::exp(log_posterior - log_posterior.max()));
        paternity[sib] = new_parentage[0]; maternity[sib] = new_parentage[1];
        deviance.at(sib) = log_likelihood.at(paternity[sib], maternity[sib]);
        paternity_counts[paternity[sib]]++; maternity_counts[maternity[sib]]++;
      }

      // store state, mapping samples back to original indicies
      if (thin == 0 && iter >= 0)
      {
        dropout_rate_samples.col(iter) = dropout_rate;
        mistyping_rate_samples.col(iter) = mistyping_rate;
        for (unsigned locus=0; locus<num_loci; ++locus)
        {
          allele_frequencies_samples[locus].col(iter) = allele_frequencies[locus];
          for (unsigned mother=0; mother<num_mothers; ++mother)
          {
            dropout_error_expectation.at(mothers[mother], locus) += maternal_dropout_errors.at(mother, locus)/double(max_iter);
            mistyping_error_expectation.at(mothers[mother], locus) += maternal_mistyping_errors.at(mother, locus)/double(max_iter);
            genotypes_expectation[locus].at(maternal_genotypes.at(0, mother, locus)-1, 
                maternal_genotypes.at(1, mother, locus)-1, mothers[mother]) += 1./double(max_iter);
          }
          for (unsigned father=0; father<num_fathers; ++father)
          {
            dropout_error_expectation.at(fathers[father], locus) += paternal_dropout_errors.at(father, locus)/double(max_iter);
            mistyping_error_expectation.at(fathers[father], locus) += paternal_mistyping_errors.at(father, locus)/double(max_iter);
            genotypes_expectation[locus].at(paternal_genotypes.at(0, father, locus)-1, 
                paternal_genotypes.at(0, father, locus)-1, fathers[father]) += 1./double(max_iter);
          }
          for (unsigned sib=0; sib<num_offspring; ++sib)
          {
            dropout_error_expectation.at(offspring[sib], locus) += offspring_dropout_errors.at(sib, locus)/double(max_iter);
            mistyping_error_expectation.at(offspring[sib], locus) += offspring_mistyping_errors.at(sib, locus)/double(max_iter);
            genotypes_expectation[locus].at(offspring_genotypes.at(0, sib, locus)-1, 
                offspring_genotypes.at(1, sib, locus)-1, offspring[sib]) += 1./double(max_iter);
          }
        }
        for (unsigned sib=0; sib<num_offspring; ++sib)
        {
          paternity_samples.at(offspring[sib], iter) = fathers[paternity[sib]]+1;
          maternity_samples.at(offspring[sib], iter) = mothers[maternity[sib]]+1;
          deviance_samples.at(offspring[sib], iter) = deviance.at(sib);
        }

        // progress
        if (iter % 100 == 0) Rcpp::Rcout << "[" << iter << "] " << "deviance: " << arma::accu(deviance) << std::endl;
      }
    }
  }

  // Rcpp collapses dimensions for elements in std::vector, so copy these over to Rcpp::List
  Rcpp::List genotypes_expectation_wrapped = Rcpp::List::create();
  Rcpp::List allele_frequencies_samples_wrapped = Rcpp::List::create();
  for (unsigned locus=0; locus<num_loci; ++locus) {
    genotypes_expectation_wrapped.push_back(genotypes_expectation[locus]);
    allele_frequencies_samples_wrapped.push_back(allele_frequencies_samples[locus]);
  }

  return Rcpp::List::create(
    Rcpp::_["paternity"] = paternity_samples,
    Rcpp::_["maternity"] = maternity_samples,
    Rcpp::_["dropout_rate"] = dropout_rate_samples,
    Rcpp::_["mistyping_rate"] = mistyping_rate_samples,
    Rcpp::_["dropout_errors"] = dropout_error_expectation,
    Rcpp::_["mistyping_errors"] = mistyping_error_expectation,
    Rcpp::_["genotypes"] = genotypes_expectation_wrapped,
    Rcpp::_["allele_frequencies"] = allele_frequencies_samples_wrapped,
    Rcpp::_["holdout_deviance"] = holdout_deviance_samples,
    Rcpp::_["deviance"] = deviance_samples,
    Rcpp::_["allele_lengths"] = allele_lengths);
}

// --------- yet another attempt, now allowing the number of mothers to vary -------- //

std::tuple<arma::umat, arma::uvec, arma::uvec, arma::uvec, arma::uvec, arma::uvec> 
sample_genotypes_given_parentage
 (const arma::uvec& paternity,
  const arma::uvec& maternity,
  const arma::umat& offspring_phenotypes, 
  const arma::uvec& maternal_phenotype, 
  const arma::vec& allele_frequencies, 
  const double& dropout_rate, 
  const double& mistyping_rate)
{
  // assumes a single "known" maternal phenotype; index 0 in "maternity" refers to the associated mother
  
  const unsigned number_of_alleles = allele_frequencies.n_elem;
  const arma::uvec fathers = arma::unique(paternity);
  const arma::uvec mothers = arma::unique(arma::join_vert(arma::uvec({0}), maternity)); //always include 0'th index, corresponding to maternal phenotype
  const arma::vec allele_frequencies_normalized = allele_frequencies / arma::accu(allele_frequencies);

  if (offspring_phenotypes.n_rows != 2) Rcpp::stop("offspring phenotypes must have 2 rows");
  if (offspring_phenotypes.n_cols != paternity.n_elem) Rcpp::stop("offspring phenotypes must have column for each individual");
  if (paternity.n_elem != maternity.n_elem) Rcpp::stop("maternity/paternity vectors must be the same length");
  if (maternal_phenotype.n_elem != 2) Rcpp::stop("maternal phenotype must have 2 elements");
  if (offspring_phenotypes.max() > number_of_alleles) Rcpp::stop("offspring allele out of range");
  if (maternal_phenotype.max() > number_of_alleles) Rcpp::stop("maternal allele out of range");
  if (arma::any(allele_frequencies_normalized < 0.)) Rcpp::stop("negative allele frequencies");
  if (dropout_rate <= 0. || mistyping_rate <= 0.) Rcpp::stop("negative genotyping error rates");

  // alternatively pass in as mutable argument
  arma::umat maternal_genotypes (2, mothers.n_elem);
  arma::umat offspring_genotypes (2, paternity.n_elem);
  arma::uvec paternal_genotypes (fathers.n_elem);
  maternal_genotypes.fill(arma::datum::nan);
  offspring_genotypes.fill(arma::datum::nan);
  paternal_genotypes.fill(arma::datum::nan);

  // tabulate sib groups
  arma::umat offspring_counts (fathers.max()+1, mothers.max()+1, arma::fill::zeros);
  for (unsigned sib=0; sib<paternity.n_elem; ++sib)
  {
    offspring_counts.at(paternity[sib],maternity[sib])++;
  }

  // simulate maternal genotype; there are choose(k,2)+k possible genotypes
  // marginalize over paternal & offspring genotypes to calculate posterior genotype probabilities
  unsigned number_of_genotypes = number_of_alleles*(number_of_alleles + 1)/2;
  for (auto mother : mothers)
  {
    arma::uvec mated_fathers = arma::find(offspring_counts.col(mother) > 0);
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
        if (mother == 0 && arma::prod(maternal_phenotype)) //index 0 is "phenotyped" mother
        {
          double maternal_phenotype_probability = 
            genotyping_error_model(maternal_phenotype, w, v, number_of_alleles, dropout_rate, mistyping_rate);
          log_halfsib_likelihood += log(maternal_phenotype_probability);
        }
        for (auto father : mated_fathers)
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
    maternal_genotypes.col(mother) = 
      possible_maternal_genotypes.col(sample(arma::exp(maternal_genotype_posterior)));
  }

  // simulate paternal genotype; there are k possible genotypes
  // paternal genotypes are conditionally independent with fixed maternal genotype
  // calculate posterior genotype probabilities by marginalizing over offspring genotypes
  for (auto father : fathers)
  {
    arma::vec paternal_genotype_posterior (number_of_alleles);
    arma::uvec possible_paternal_genotypes (number_of_alleles);
    arma::uvec offspring_from_father = arma::find(paternity == father);
    for (unsigned u=1; u<=number_of_alleles; ++u) // paternal allele
    {
      double paternal_genotype_probability = 
          allele_frequencies_normalized[u-1]; //hwe prior
      double log_fullsib_likelihood = log(paternal_genotype_probability);
      for (auto offspring : offspring_from_father)
      {
        arma::uvec maternal_genotype = maternal_genotypes.col(maternity[offspring]);
        arma::uvec offspring_phenotype = offspring_phenotypes.col(offspring);
        if (arma::prod(offspring_phenotype)) { 
          double offspring_phenotype_probability = // Mendelian segregation probs * phenotype probabilities
            0.5 * genotyping_error_model(offspring_phenotype, maternal_genotype[0], u, number_of_alleles, dropout_rate, mistyping_rate) + 
            0.5 * genotyping_error_model(offspring_phenotype, maternal_genotype[1], u, number_of_alleles, dropout_rate, mistyping_rate); 
          log_fullsib_likelihood += log(offspring_phenotype_probability);
        }
      }
      paternal_genotype_posterior[u-1] = log_fullsib_likelihood;
      possible_paternal_genotypes[u-1] = u;
    }
    paternal_genotype_posterior -= paternal_genotype_posterior.max();
    paternal_genotypes.at(father) = possible_paternal_genotypes.at(sample(arma::exp(paternal_genotype_posterior)));
  }

  // simulate offspring genotypes
  for (unsigned sib=0; sib<paternity.n_elem; ++sib)
  {
    arma::uvec maternal_genotype = maternal_genotypes.col(maternity[sib]);
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
  //0-index is the mother, remaining indices are offspring
  arma::uvec dropouts (paternity.n_elem+1, arma::fill::zeros);
  arma::uvec mistypes (paternity.n_elem+1, arma::fill::zeros);
  arma::uvec heterozygous (paternity.n_elem+1, arma::fill::zeros);
  arma::uvec nonmissing (paternity.n_elem+1, arma::fill::zeros);
  if (arma::prod(maternal_phenotype))
  {
    arma::uvec maternal_genotype = maternal_genotypes.col(0); //genotype of phenotyped mother
    arma::uvec errors =
      simulate_genotyping_errors(maternal_phenotype, maternal_genotype[0], 
          maternal_genotype[1], number_of_alleles, dropout_rate, mistyping_rate);
    dropouts.at(0) = errors.at(0);
    mistypes.at(0) = errors.at(1);
    nonmissing.at(0) += 2;
    if (maternal_genotype[0] != maternal_genotype[1]) heterozygous.at(0) += 1;
  }
  for (unsigned sib=0; sib<paternity.n_elem; ++sib)
  {
    arma::uvec offspring_phenotype = offspring_phenotypes.col(sib);
    if (arma::prod(offspring_phenotype)) { 
      arma::uvec errors =
        simulate_genotyping_errors(offspring_phenotype, offspring_genotypes.at(0,sib), 
          offspring_genotypes.at(1,sib), number_of_alleles, dropout_rate, mistyping_rate);
      dropouts.at(sib+1) = errors.at(0);
      mistypes.at(sib+1) = errors.at(1);
      nonmissing.at(sib+1) += 2;
      if (offspring_genotypes.at(0,sib) != offspring_genotypes.at(1,sib)) heterozygous.at(sib+1) += 1;
    }
  }

  // tally maternal/paternal alleles
  arma::uvec allele_counts (number_of_alleles, arma::fill::zeros);
  for (unsigned mother=0; mother<mothers.n_elem; ++mother)
  {
    allele_counts[maternal_genotypes.at(0,mother)-1]++;
    allele_counts[maternal_genotypes.at(1,mother)-1]++;
  }
  for (unsigned father=0; father<fathers.n_elem; ++father)
  {
    allele_counts[paternal_genotypes[father]-1]++;
  }

  // output
  arma::umat genotypes = arma::join_horiz(maternal_genotypes.col(0), offspring_genotypes);
  return std::make_tuple(genotypes, allele_counts, dropouts, heterozygous, mistypes, nonmissing);
}

double parentage_loglikelihood_by_locus 
 (const arma::uvec& paternity,
  const arma::uvec& maternity,
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
  const arma::uvec mothers = arma::unique(arma::join_vert(arma::uvec({0}), maternity)); //always include 0'th index, corresponding to maternal phenotype
  const arma::vec allele_frequencies_normalized = allele_frequencies / arma::accu(allele_frequencies);

  if (paternity.n_elem != maternity.n_elem) Rcpp::stop("maternity/paternity vectors must be the same length");
  if (offspring_phenotypes.n_rows != 2) Rcpp::stop("offspring phenotypes must have 2 rows");
  if (offspring_phenotypes.n_cols != paternity.n_elem) Rcpp::stop("offspring phenotypes must have column for each individual");
  if (maternal_phenotype.n_elem != 2) Rcpp::stop("maternal phenotype must have 2 elements");
  if (offspring_phenotypes.max() > number_of_alleles) Rcpp::stop("offspring allele out of range");
  if (maternal_phenotype.max() > number_of_alleles) Rcpp::stop("maternal allele out of range");
  if (arma::any(allele_frequencies_normalized < 0.)) Rcpp::stop("negative allele frequencies");
  if (dropout_rate <= 0. || mistyping_rate <= 0.) Rcpp::stop("negative genotyping error rates");

  // tabulate sib groups
  arma::umat offspring_counts (fathers.max()+1, mothers.max()+1, arma::fill::zeros);
  for (unsigned sib=0; sib<paternity.n_elem; ++sib)
  {
    offspring_counts.at(paternity[sib],maternity[sib])++;
  }

  double log_likelihood = 0;
  for (auto mother : mothers)
  {
    arma::uvec mated_fathers = arma::find(offspring_counts.col(mother) > 0);
    double halfsib_likelihood = 0.;
    double halfsib_running_maximum = -arma::datum::inf;
    for (unsigned w=1; w<=number_of_alleles; ++w) // first maternal allele 
    { 
      for (unsigned v=w; v<=number_of_alleles; ++v) // second maternal allele
      {
        double maternal_genotype_probability = 
          (2.-int(w==v)) * allele_frequencies_normalized[w-1] * allele_frequencies_normalized[v-1]; //hwe prior
        double log_halfsib_likelihood = log(maternal_genotype_probability); 
        if (mother == 0 && arma::prod(maternal_phenotype)) //0'th mother is phenotyped
        {
          double maternal_phenotype_probability = 
            genotyping_error_model(maternal_phenotype, w, v, number_of_alleles, dropout_rate, mistyping_rate);
          log_halfsib_likelihood += log(maternal_phenotype_probability);
        }
        for (auto father : mated_fathers)
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
    log_likelihood += log(halfsib_likelihood) + halfsib_running_maximum;
  }
  return log_likelihood;
}

double parentage_loglikelihood 
 (arma::uvec paternity, 
  arma::uvec maternity,
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
      parentage_loglikelihood_by_locus(paternity, maternity, offspring_phenotypes.slice(locus), 
          maternal_phenotype.col(locus), allele_frequencies[locus], dropout_rate[locus], mistyping_rate[locus]);
  }
  return log_likelihood;
}

// [[Rcpp::export]]
Rcpp::List sample_parentage_and_error_rates
 (arma::ucube phenotypes, 
  const unsigned mother = 1,
  const unsigned burn_in = 0,
  const unsigned thinning_interval = 1,
  const unsigned number_of_mcmc_samples = 1000,
  const bool global_genotyping_error_rates = true,
  const double concentration = 1.,
  const double starting_error_rate = 0.01)
{
  // samples from posterior distribution of full sib groups with Dirichlet process prior,
  // using algorithm 8 from Neal 2000 JCGS
  
  if (mother > phenotypes.n_cols || mother < 1) Rcpp::stop("1-based index of mother out of range");

  const unsigned max_iter = number_of_mcmc_samples;
  const unsigned num_loci = phenotypes.n_slices;
  const unsigned num_offspring = phenotypes.n_cols - 1;

  // priors (hardcoded for now)
  const double alpha = std::fabs(concentration); //dirichlet process concentration parameter
  const double delta = 1.; //allele frequency concentration parameter
  const arma::vec dropout_rate_prior = {{1.,1.}}; //beta(number of dropout homozygotes, number of heterozygotes)
  const arma::vec mistyping_rate_prior = {{1.,1.}}; //beta(number of mistypes, number of correct calls)
  std::vector<arma::uvec> allele_lengths = unique_alleles(phenotypes, false);
  std::vector<arma::vec> allele_frequencies = collapse_alleles_and_generate_genotype_prior(phenotypes, false);

  // split maternal, offspring phenotypes
  arma::umat maternal_phenotype = phenotypes.tube(arma::span::all, arma::span(mother-1));
  arma::ucube offspring_phenotypes = phenotypes; offspring_phenotypes.shed_col(mother-1);

  // initialize (could draw from prior instead)
  arma::uvec paternity = arma::zeros<arma::uvec>(num_offspring);//this could be randomized
  arma::uvec maternity = arma::zeros<arma::uvec>(num_offspring);//this could be randomized, but must include 0
  arma::vec dropout_rate (num_loci); dropout_rate.fill(0.01);
  arma::vec mistyping_rate (num_loci); mistyping_rate.fill(0.01);

  // storage
  arma::imat paternity_samples (num_offspring, max_iter);
  arma::imat maternity_samples (num_offspring, max_iter);
  arma::mat dropout_rate_samples (num_loci, max_iter);
  arma::mat mistyping_rate_samples (num_loci, max_iter);
  arma::vec deviance_samples (max_iter);
  arma::mat dropout_errors (num_offspring+1, num_loci, arma::fill::zeros);
  arma::mat mistyping_errors (num_offspring+1, num_loci, arma::fill::zeros);
  std::vector<arma::cube> genotype_posterior;
  for(unsigned locus=0; locus<num_loci; ++locus) 
  {
    genotype_posterior.emplace_back(arma::cube(allele_frequencies[locus].n_elem,allele_frequencies[locus].n_elem,num_offspring+1,arma::fill::zeros));
  }

  double deviance = 0.;
  paternity = recode_to_contiguous_integers(paternity);
  for (int iter=-int(burn_in); iter<int(max_iter); ++iter)
  {
    for (unsigned thin=0; thin<thinning_interval; ++thin)
    {
      // update paternity vector
      for (unsigned sib=0; sib<num_offspring; ++sib)
      {
        // tally size of sib groups
        unsigned current_number_of_fathers = paternity.max() + 1;
        unsigned current_number_of_mothers = maternity.max() + 1;
        arma::umat offspring_per_mating (current_number_of_fathers + 1, current_number_of_mothers + 1, arma::fill::zeros);
        for (unsigned i=0; i<paternity.n_elem; ++i) offspring_per_mating.at(paternity[i],maternity[i])++;
        offspring_per_mating.at(paternity[sib],maternity[sib])--; //remove sib from group
        arma::uvec offspring_per_father = arma::sum(offspring_per_mating, 1);
        arma::urowvec offspring_per_mother = arma::sum(offspring_per_mating, 0);

        // find permissible matings under constraint of single mating/father
        arma::umat permissible_matings (current_number_of_fathers + 1, current_number_of_mothers + 1, arma::fill::zeros);
        unsigned new_fathers = 0;
        for (unsigned father=0; father<permissible_matings.n_rows; ++father)
        {
          for (unsigned mother=0; mother<permissible_matings.n_cols; ++mother)
          {
            const bool existing_father = offspring_per_father.at(father) > 0;
            const bool existing_mother = offspring_per_mother.at(mother) > 0;
            const bool existing_mating = offspring_per_mating.at(father,mother) > 0;
            permissible_matings.at(father,mother) = unsigned(
                ( existing_mating                    ) || 
                ( existing_mother && !existing_father) || 
                (!existing_mother && !existing_father) );
            if (permissible_matings.at(father,mother) && !existing_father) new_fathers++;
          }
        }

        // calculate parentage likelihoods across permissible matings
        arma::mat log_likelihood (current_number_of_fathers + 1, current_number_of_mothers + 1, arma::fill::zeros);
        arma::mat log_prior (current_number_of_fathers + 1, current_number_of_mothers + 1, arma::fill::zeros);
        for (unsigned father=0; father<log_likelihood.n_rows; ++father)
        {
          for (unsigned mother=0; mother<log_likelihood.n_cols; ++mother)
          {
            if (permissible_matings.at(father,mother))
            {
              paternity[sib] = father;
              maternity[sib] = mother;
              log_likelihood.at(father,mother) = 
                parentage_loglikelihood(paternity, maternity, offspring_phenotypes, maternal_phenotype, 
                    allele_frequencies, dropout_rate, mistyping_rate);

              // "restraunt process" prior on number of matings
              // TODO would be useful to have a way to sample from the prior
              if (alpha > 0.)
              {
                log_prior.at(father,mother) = offspring_per_mating.at(father,mother) > 0 ?
                  log(double(offspring_per_mating.at(father,mother))) - log(double(num_offspring)-1.+alpha): 
                  log(alpha) - log(double(new_fathers)) - log(double(num_offspring)-1.+alpha);
              }
            } else { log_likelihood.at(father,mother) = -arma::datum::inf; }
          } 
        }

        // sample new father
        arma::mat log_conditional = log_prior + log_likelihood;
        arma::uvec new_parentage = sample_matrix(arma::exp(log_conditional - log_conditional.max()));
        paternity[sib] = new_parentage.at(0); maternity[sib] = new_parentage.at(1);
        deviance = -2 * log_likelihood.at(paternity[sib],maternity[sib]);

        // this song and dance forces the 0-index to correspond to the known mother
        paternity = recode_to_contiguous_integers(paternity); 
        arma::uvec maternity_aux = recode_to_contiguous_integers(arma::join_vert(maternity, arma::uvec({0})));
        maternity = maternity_aux.head(num_offspring);
        //why recode? indices will increase, if pre-existing singleton is moved to a father with a higher index
      }

      // update error rates and allele frequencies via data augmentation
      unsigned global_dropouts = 0, global_heterozygous = 0, global_mistypes = 0, global_nonmissing = 0;
      for (unsigned locus=0; locus<num_loci; ++locus)
      {
        arma::umat genotypes;
        arma::uvec allele_counts, dropouts, heterozygous, mistypes, nonmissing;
        std::tie(genotypes, allele_counts, dropouts, heterozygous, mistypes, nonmissing) =
          sample_genotypes_given_parentage(paternity, maternity, offspring_phenotypes.slice(locus), 
              maternal_phenotype.col(locus), allele_frequencies[locus], dropout_rate[locus], mistyping_rate[locus]);

        // track expected errors, genotypes
        if (iter >= 0 && thin == 0)
        {
          for (unsigned i=0; i<num_offspring+1; ++i)
          {
            dropout_errors.at(i,locus) += double(dropouts.at(i));
            mistyping_errors.at(i,locus) += double(mistypes.at(i));
            genotype_posterior[locus].at(genotypes.at(0,i)-1, genotypes.at(1,i)-1, i) += 1.0;//genotypes are 1-indexed so convert
          }
        }

        // update error rates
        global_dropouts += arma::accu(dropouts); global_mistypes += arma::accu(mistypes);
        global_heterozygous += arma::accu(heterozygous); global_nonmissing += arma::accu(nonmissing);
        dropout_rate[locus] = 0.5 * R::rbeta(1. + arma::accu(dropouts), 1. + arma::accu(heterozygous) - arma::accu(dropouts));
        mistyping_rate[locus] = R::rbeta(1. + arma::accu(mistypes), 1. + arma::accu(nonmissing) - arma::accu(mistypes));

        // update allele frequencies
        for (unsigned allele=0; allele<allele_frequencies[locus].n_elem; ++allele)
        {
          allele_frequencies[locus][allele] = R::rgamma(1. + allele_counts[allele], 1.);
        }
        allele_frequencies[locus] /= arma::accu(allele_frequencies[locus]);
      }
      if (global_genotyping_error_rates)
      {
        // overwrite per-locus rates with global rate
        dropout_rate.fill(0.5 * R::rbeta(1. + global_dropouts, 1. + global_heterozygous - global_dropouts));
        mistyping_rate.fill(R::rbeta(1. + global_mistypes, 1. + global_nonmissing - global_mistypes));
      }

      // store state
      if (iter >= 0 && thin == 0)
      {
        paternity_samples.col(iter) = arma::conv_to<arma::ivec>::from(paternity);
        maternity_samples.col(iter) = arma::conv_to<arma::ivec>::from(maternity);
        dropout_rate_samples.col(iter) = dropout_rate;
        mistyping_rate_samples.col(iter) = mistyping_rate;
        deviance_samples.at(iter) = deviance;
      }

      if (thin == 0 && iter >= 0 && iter % 100 == 0) Rcpp::Rcout << "sampling [" << iter << "] " << "deviance: " << deviance << std::endl;
    }
  }

  // MAP estimates for genotypes, reorder to mirror input
  arma::ucube imputed_genotypes (arma::size(phenotypes));
  for (unsigned locus=0; locus<num_loci; ++locus)
  {
    arma::umat genotypes (2, num_offspring+1);
    for (unsigned i=0; i<num_offspring+1; ++i)
    {
      genotypes.at(0, i) = allele_lengths[locus].at(arma::index_max(arma::sum(genotype_posterior[locus].slice(i), 0)));
      genotypes.at(1, i) = allele_lengths[locus].at(arma::index_max(arma::sum(genotype_posterior[locus].slice(i), 1)));
    }
    arma::uvec maternal_genotype = genotypes.col(0);
    genotypes.shed_col(0); genotypes.insert_cols(mother-1, maternal_genotype);
    imputed_genotypes.slice(locus) = genotypes;
  }

  // posterior expectation of error counts & reorder matrices to mirror input (0'th index is mother, offspring are 1-based indices)
  dropout_errors /= double(max_iter);
  arma::rowvec maternal_dropout_errors = dropout_errors.row(0); 
  dropout_errors.shed_row(0); dropout_errors.insert_rows(mother-1, maternal_dropout_errors);

  mistyping_errors /= double(max_iter);
  arma::rowvec maternal_mistyping_errors = mistyping_errors.row(0); 
  mistyping_errors.shed_row(0); mistyping_errors.insert_rows(mother-1, maternal_mistyping_errors);

  arma::irowvec maternal_parentage (max_iter); maternal_parentage.fill(arma::datum::nan);
  maternity_samples.insert_rows(mother-1, maternal_parentage);
  paternity_samples.insert_rows(mother-1, maternal_parentage);

  return Rcpp::List::create(
    //Rcpp::_["allele_lengths"] = allele_lengths,
    Rcpp::_["paternity"] = paternity_samples,
    Rcpp::_["maternity"] = maternity_samples,
    Rcpp::_["dropout_rate"] = dropout_rate_samples,
    Rcpp::_["mistyping_rate"] = mistyping_rate_samples,
    Rcpp::_["dropout_errors"] = dropout_errors,
    Rcpp::_["mistyping_errors"] = mistyping_errors,
    Rcpp::_["imputed_genotypes"] = imputed_genotypes,
    Rcpp::_["deviance"] = deviance_samples);
}
