// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// log_ascending_factorial
double log_ascending_factorial(const double x, const unsigned r);
RcppExport SEXP _sydneyPaternity_log_ascending_factorial(SEXP xSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type x(xSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(log_ascending_factorial(x, r));
    return rcpp_result_gen;
END_RCPP
}
// log_descending_factorial
double log_descending_factorial(const double x, const unsigned r);
RcppExport SEXP _sydneyPaternity_log_descending_factorial(SEXP xSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type x(xSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(log_descending_factorial(x, r));
    return rcpp_result_gen;
END_RCPP
}
// log_uniform_MFM_prior
double log_uniform_MFM_prior(const unsigned n, const unsigned t, const double gamma, const unsigned max_number_of_components);
RcppExport SEXP _sydneyPaternity_log_uniform_MFM_prior(SEXP nSEXP, SEXP tSEXP, SEXP gammaSEXP, SEXP max_number_of_componentsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned >::type n(nSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type max_number_of_components(max_number_of_componentsSEXP);
    rcpp_result_gen = Rcpp::wrap(log_uniform_MFM_prior(n, t, gamma, max_number_of_components));
    return rcpp_result_gen;
END_RCPP
}
// genotyping_error_model
double genotyping_error_model(const arma::uvec& phenotype, const unsigned& genotype0, const unsigned& genotype1, const unsigned& number_of_alleles, const double& dropout_rate, const double& mistyping_rate);
RcppExport SEXP _sydneyPaternity_genotyping_error_model(SEXP phenotypeSEXP, SEXP genotype0SEXP, SEXP genotype1SEXP, SEXP number_of_allelesSEXP, SEXP dropout_rateSEXP, SEXP mistyping_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type phenotype(phenotypeSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type genotype0(genotype0SEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type genotype1(genotype1SEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type number_of_alleles(number_of_allelesSEXP);
    Rcpp::traits::input_parameter< const double& >::type dropout_rate(dropout_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type mistyping_rate(mistyping_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(genotyping_error_model(phenotype, genotype0, genotype1, number_of_alleles, dropout_rate, mistyping_rate));
    return rcpp_result_gen;
END_RCPP
}
// genotyping_error_model_class
int genotyping_error_model_class(const arma::uvec& phenotype, const unsigned& genotype0, const unsigned& genotype1);
RcppExport SEXP _sydneyPaternity_genotyping_error_model_class(SEXP phenotypeSEXP, SEXP genotype0SEXP, SEXP genotype1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type phenotype(phenotypeSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type genotype0(genotype0SEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type genotype1(genotype1SEXP);
    rcpp_result_gen = Rcpp::wrap(genotyping_error_model_class(phenotype, genotype0, genotype1));
    return rcpp_result_gen;
END_RCPP
}
// simulate_genotyping_errors
arma::uvec simulate_genotyping_errors(const arma::uvec& phenotype, const unsigned& genotype0, const unsigned& genotype1, const unsigned& number_of_alleles, const double& dropout_rate, const double& mistyping_rate);
RcppExport SEXP _sydneyPaternity_simulate_genotyping_errors(SEXP phenotypeSEXP, SEXP genotype0SEXP, SEXP genotype1SEXP, SEXP number_of_allelesSEXP, SEXP dropout_rateSEXP, SEXP mistyping_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type phenotype(phenotypeSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type genotype0(genotype0SEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type genotype1(genotype1SEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type number_of_alleles(number_of_allelesSEXP);
    Rcpp::traits::input_parameter< const double& >::type dropout_rate(dropout_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type mistyping_rate(mistyping_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_genotyping_errors(phenotype, genotype0, genotype1, number_of_alleles, dropout_rate, mistyping_rate));
    return rcpp_result_gen;
END_RCPP
}
// sample_error_rates_given_paternity
Rcpp::List sample_error_rates_given_paternity(arma::ucube phenotypes, arma::uvec paternity, const unsigned mother, const unsigned number_of_mcmc_samples, const unsigned global_genotyping_error_rates, const bool random_allele_frequencies, const bool add_unsampled_allele);
RcppExport SEXP _sydneyPaternity_sample_error_rates_given_paternity(SEXP phenotypesSEXP, SEXP paternitySEXP, SEXP motherSEXP, SEXP number_of_mcmc_samplesSEXP, SEXP global_genotyping_error_ratesSEXP, SEXP random_allele_frequenciesSEXP, SEXP add_unsampled_alleleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ucube >::type phenotypes(phenotypesSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type paternity(paternitySEXP);
    Rcpp::traits::input_parameter< const unsigned >::type mother(motherSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type number_of_mcmc_samples(number_of_mcmc_samplesSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type global_genotyping_error_rates(global_genotyping_error_ratesSEXP);
    Rcpp::traits::input_parameter< const bool >::type random_allele_frequencies(random_allele_frequenciesSEXP);
    Rcpp::traits::input_parameter< const bool >::type add_unsampled_allele(add_unsampled_alleleSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_error_rates_given_paternity(phenotypes, paternity, mother, number_of_mcmc_samples, global_genotyping_error_rates, random_allele_frequencies, add_unsampled_allele));
    return rcpp_result_gen;
END_RCPP
}
// optimize_paternity_given_error_rates
Rcpp::List optimize_paternity_given_error_rates(arma::ucube phenotypes, arma::vec dropout_rate, arma::vec mistyping_rate, const unsigned mother);
RcppExport SEXP _sydneyPaternity_optimize_paternity_given_error_rates(SEXP phenotypesSEXP, SEXP dropout_rateSEXP, SEXP mistyping_rateSEXP, SEXP motherSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ucube >::type phenotypes(phenotypesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dropout_rate(dropout_rateSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mistyping_rate(mistyping_rateSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type mother(motherSEXP);
    rcpp_result_gen = Rcpp::wrap(optimize_paternity_given_error_rates(phenotypes, dropout_rate, mistyping_rate, mother));
    return rcpp_result_gen;
END_RCPP
}
// loglikelihood_of_error_rates_given_paternity
arma::mat loglikelihood_of_error_rates_given_paternity(arma::ucube phenotypes, arma::uvec paternity, arma::mat grid_of_error_rates, const unsigned mother);
RcppExport SEXP _sydneyPaternity_loglikelihood_of_error_rates_given_paternity(SEXP phenotypesSEXP, SEXP paternitySEXP, SEXP grid_of_error_ratesSEXP, SEXP motherSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ucube >::type phenotypes(phenotypesSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type paternity(paternitySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type grid_of_error_rates(grid_of_error_ratesSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type mother(motherSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikelihood_of_error_rates_given_paternity(phenotypes, paternity, grid_of_error_rates, mother));
    return rcpp_result_gen;
END_RCPP
}
// collapse_alleles_and_generate_prior_wrapper
Rcpp::List collapse_alleles_and_generate_prior_wrapper(arma::ucube phenotypes, const unsigned mother, const bool add_unsampled_allele);
RcppExport SEXP _sydneyPaternity_collapse_alleles_and_generate_prior_wrapper(SEXP phenotypesSEXP, SEXP motherSEXP, SEXP add_unsampled_alleleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ucube >::type phenotypes(phenotypesSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type mother(motherSEXP);
    Rcpp::traits::input_parameter< const bool >::type add_unsampled_allele(add_unsampled_alleleSEXP);
    rcpp_result_gen = Rcpp::wrap(collapse_alleles_and_generate_prior_wrapper(phenotypes, mother, add_unsampled_allele));
    return rcpp_result_gen;
END_RCPP
}
// sample_paternity_and_error_rates_from_joint_posterior
Rcpp::List sample_paternity_and_error_rates_from_joint_posterior(arma::ucube phenotypes, const unsigned mother, const unsigned number_of_mcmc_samples, const bool global_genotyping_error_rates, const double concentration);
RcppExport SEXP _sydneyPaternity_sample_paternity_and_error_rates_from_joint_posterior(SEXP phenotypesSEXP, SEXP motherSEXP, SEXP number_of_mcmc_samplesSEXP, SEXP global_genotyping_error_ratesSEXP, SEXP concentrationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ucube >::type phenotypes(phenotypesSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type mother(motherSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type number_of_mcmc_samples(number_of_mcmc_samplesSEXP);
    Rcpp::traits::input_parameter< const bool >::type global_genotyping_error_rates(global_genotyping_error_ratesSEXP);
    Rcpp::traits::input_parameter< const double >::type concentration(concentrationSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_paternity_and_error_rates_from_joint_posterior(phenotypes, mother, number_of_mcmc_samples, global_genotyping_error_rates, concentration));
    return rcpp_result_gen;
END_RCPP
}
// sample_matrix
arma::uvec sample_matrix(arma::mat probabilities);
RcppExport SEXP _sydneyPaternity_sample_matrix(SEXP probabilitiesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type probabilities(probabilitiesSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_matrix(probabilities));
    return rcpp_result_gen;
END_RCPP
}
// select_columns_from_cube
arma::ucube select_columns_from_cube(arma::ucube input, arma::uvec which);
RcppExport SEXP _sydneyPaternity_select_columns_from_cube(SEXP inputSEXP, SEXP whichSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ucube >::type input(inputSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type which(whichSEXP);
    rcpp_result_gen = Rcpp::wrap(select_columns_from_cube(input, which));
    return rcpp_result_gen;
END_RCPP
}
// phenotype_error_model
double phenotype_error_model(const arma::uvec& phenotype, const arma::uvec& genotype, const unsigned& number_of_alleles, const double& dropout_rate, const double& mistyping_rate);
RcppExport SEXP _sydneyPaternity_phenotype_error_model(SEXP phenotypeSEXP, SEXP genotypeSEXP, SEXP number_of_allelesSEXP, SEXP dropout_rateSEXP, SEXP mistyping_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type phenotype(phenotypeSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type genotype(genotypeSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type number_of_alleles(number_of_allelesSEXP);
    Rcpp::traits::input_parameter< const double& >::type dropout_rate(dropout_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type mistyping_rate(mistyping_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(phenotype_error_model(phenotype, genotype, number_of_alleles, dropout_rate, mistyping_rate));
    return rcpp_result_gen;
END_RCPP
}
// sample_phenotype_errors
arma::uvec sample_phenotype_errors(const arma::uvec& phenotype, const arma::uvec& genotype, const unsigned& number_of_alleles, const double& dropout_rate, const double& mistyping_rate);
RcppExport SEXP _sydneyPaternity_sample_phenotype_errors(SEXP phenotypeSEXP, SEXP genotypeSEXP, SEXP number_of_allelesSEXP, SEXP dropout_rateSEXP, SEXP mistyping_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type phenotype(phenotypeSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type genotype(genotypeSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type number_of_alleles(number_of_allelesSEXP);
    Rcpp::traits::input_parameter< const double& >::type dropout_rate(dropout_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type mistyping_rate(mistyping_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_phenotype_errors(phenotype, genotype, number_of_alleles, dropout_rate, mistyping_rate));
    return rcpp_result_gen;
END_RCPP
}
// mendelian_genotype_model
double mendelian_genotype_model(const arma::uvec& offspring_phenotype, const arma::uvec& maternal_genotype, const arma::uvec& paternal_genotype, const unsigned& number_of_alleles, const double& dropout_rate, const double& mistyping_rate);
RcppExport SEXP _sydneyPaternity_mendelian_genotype_model(SEXP offspring_phenotypeSEXP, SEXP maternal_genotypeSEXP, SEXP paternal_genotypeSEXP, SEXP number_of_allelesSEXP, SEXP dropout_rateSEXP, SEXP mistyping_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type offspring_phenotype(offspring_phenotypeSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type maternal_genotype(maternal_genotypeSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type paternal_genotype(paternal_genotypeSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type number_of_alleles(number_of_allelesSEXP);
    Rcpp::traits::input_parameter< const double& >::type dropout_rate(dropout_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type mistyping_rate(mistyping_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(mendelian_genotype_model(offspring_phenotype, maternal_genotype, paternal_genotype, number_of_alleles, dropout_rate, mistyping_rate));
    return rcpp_result_gen;
END_RCPP
}
// sample_mendelian_genotype
arma::uvec sample_mendelian_genotype(const arma::uvec& offspring_phenotype, const arma::uvec& maternal_genotype, const arma::uvec& paternal_genotype, const unsigned& number_of_alleles, const double& dropout_rate, const double& mistyping_rate);
RcppExport SEXP _sydneyPaternity_sample_mendelian_genotype(SEXP offspring_phenotypeSEXP, SEXP maternal_genotypeSEXP, SEXP paternal_genotypeSEXP, SEXP number_of_allelesSEXP, SEXP dropout_rateSEXP, SEXP mistyping_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uvec& >::type offspring_phenotype(offspring_phenotypeSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type maternal_genotype(maternal_genotypeSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type paternal_genotype(paternal_genotypeSEXP);
    Rcpp::traits::input_parameter< const unsigned& >::type number_of_alleles(number_of_allelesSEXP);
    Rcpp::traits::input_parameter< const double& >::type dropout_rate(dropout_rateSEXP);
    Rcpp::traits::input_parameter< const double& >::type mistyping_rate(mistyping_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_mendelian_genotype(offspring_phenotype, maternal_genotype, paternal_genotype, number_of_alleles, dropout_rate, mistyping_rate));
    return rcpp_result_gen;
END_RCPP
}
// sample_parentage_and_error_rates_from_joint_posterior
Rcpp::List sample_parentage_and_error_rates_from_joint_posterior(arma::ucube phenotypes, arma::uvec mothers, arma::uvec fathers, arma::uvec holdouts, const unsigned number_of_mcmc_samples, const unsigned burn_in_samples, const unsigned thinning_interval, const bool global_genotyping_error_rates);
RcppExport SEXP _sydneyPaternity_sample_parentage_and_error_rates_from_joint_posterior(SEXP phenotypesSEXP, SEXP mothersSEXP, SEXP fathersSEXP, SEXP holdoutsSEXP, SEXP number_of_mcmc_samplesSEXP, SEXP burn_in_samplesSEXP, SEXP thinning_intervalSEXP, SEXP global_genotyping_error_ratesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ucube >::type phenotypes(phenotypesSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type mothers(mothersSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fathers(fathersSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type holdouts(holdoutsSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type number_of_mcmc_samples(number_of_mcmc_samplesSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type burn_in_samples(burn_in_samplesSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type thinning_interval(thinning_intervalSEXP);
    Rcpp::traits::input_parameter< const bool >::type global_genotyping_error_rates(global_genotyping_error_ratesSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_parentage_and_error_rates_from_joint_posterior(phenotypes, mothers, fathers, holdouts, number_of_mcmc_samples, burn_in_samples, thinning_interval, global_genotyping_error_rates));
    return rcpp_result_gen;
END_RCPP
}
// sample_parentage_and_error_rates_from_joint_posterior_alt
Rcpp::List sample_parentage_and_error_rates_from_joint_posterior_alt(arma::ucube phenotypes, arma::uvec mothers, arma::uvec fathers, const double concentration, const unsigned number_of_mcmc_samples, const unsigned burn_in_samples, const unsigned thinning_interval, const bool global_genotyping_error_rates, const bool sample_from_prior, const bool random_initialization);
RcppExport SEXP _sydneyPaternity_sample_parentage_and_error_rates_from_joint_posterior_alt(SEXP phenotypesSEXP, SEXP mothersSEXP, SEXP fathersSEXP, SEXP concentrationSEXP, SEXP number_of_mcmc_samplesSEXP, SEXP burn_in_samplesSEXP, SEXP thinning_intervalSEXP, SEXP global_genotyping_error_ratesSEXP, SEXP sample_from_priorSEXP, SEXP random_initializationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ucube >::type phenotypes(phenotypesSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type mothers(mothersSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type fathers(fathersSEXP);
    Rcpp::traits::input_parameter< const double >::type concentration(concentrationSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type number_of_mcmc_samples(number_of_mcmc_samplesSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type burn_in_samples(burn_in_samplesSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type thinning_interval(thinning_intervalSEXP);
    Rcpp::traits::input_parameter< const bool >::type global_genotyping_error_rates(global_genotyping_error_ratesSEXP);
    Rcpp::traits::input_parameter< const bool >::type sample_from_prior(sample_from_priorSEXP);
    Rcpp::traits::input_parameter< const bool >::type random_initialization(random_initializationSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_parentage_and_error_rates_from_joint_posterior_alt(phenotypes, mothers, fathers, concentration, number_of_mcmc_samples, burn_in_samples, thinning_interval, global_genotyping_error_rates, sample_from_prior, random_initialization));
    return rcpp_result_gen;
END_RCPP
}
// sample_parentage_and_error_rates
Rcpp::List sample_parentage_and_error_rates(arma::ucube phenotypes, const unsigned mother, const unsigned burn_in, const unsigned thinning_interval, const unsigned number_of_mcmc_samples, const bool global_genotyping_error_rates, const double concentration, const double starting_error_rate);
RcppExport SEXP _sydneyPaternity_sample_parentage_and_error_rates(SEXP phenotypesSEXP, SEXP motherSEXP, SEXP burn_inSEXP, SEXP thinning_intervalSEXP, SEXP number_of_mcmc_samplesSEXP, SEXP global_genotyping_error_ratesSEXP, SEXP concentrationSEXP, SEXP starting_error_rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ucube >::type phenotypes(phenotypesSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type mother(motherSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type thinning_interval(thinning_intervalSEXP);
    Rcpp::traits::input_parameter< const unsigned >::type number_of_mcmc_samples(number_of_mcmc_samplesSEXP);
    Rcpp::traits::input_parameter< const bool >::type global_genotyping_error_rates(global_genotyping_error_ratesSEXP);
    Rcpp::traits::input_parameter< const double >::type concentration(concentrationSEXP);
    Rcpp::traits::input_parameter< const double >::type starting_error_rate(starting_error_rateSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_parentage_and_error_rates(phenotypes, mother, burn_in, thinning_interval, number_of_mcmc_samples, global_genotyping_error_rates, concentration, starting_error_rate));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sydneyPaternity_log_ascending_factorial", (DL_FUNC) &_sydneyPaternity_log_ascending_factorial, 2},
    {"_sydneyPaternity_log_descending_factorial", (DL_FUNC) &_sydneyPaternity_log_descending_factorial, 2},
    {"_sydneyPaternity_log_uniform_MFM_prior", (DL_FUNC) &_sydneyPaternity_log_uniform_MFM_prior, 4},
    {"_sydneyPaternity_genotyping_error_model", (DL_FUNC) &_sydneyPaternity_genotyping_error_model, 6},
    {"_sydneyPaternity_genotyping_error_model_class", (DL_FUNC) &_sydneyPaternity_genotyping_error_model_class, 3},
    {"_sydneyPaternity_simulate_genotyping_errors", (DL_FUNC) &_sydneyPaternity_simulate_genotyping_errors, 6},
    {"_sydneyPaternity_sample_error_rates_given_paternity", (DL_FUNC) &_sydneyPaternity_sample_error_rates_given_paternity, 7},
    {"_sydneyPaternity_optimize_paternity_given_error_rates", (DL_FUNC) &_sydneyPaternity_optimize_paternity_given_error_rates, 4},
    {"_sydneyPaternity_loglikelihood_of_error_rates_given_paternity", (DL_FUNC) &_sydneyPaternity_loglikelihood_of_error_rates_given_paternity, 4},
    {"_sydneyPaternity_collapse_alleles_and_generate_prior_wrapper", (DL_FUNC) &_sydneyPaternity_collapse_alleles_and_generate_prior_wrapper, 3},
    {"_sydneyPaternity_sample_paternity_and_error_rates_from_joint_posterior", (DL_FUNC) &_sydneyPaternity_sample_paternity_and_error_rates_from_joint_posterior, 5},
    {"_sydneyPaternity_sample_matrix", (DL_FUNC) &_sydneyPaternity_sample_matrix, 1},
    {"_sydneyPaternity_select_columns_from_cube", (DL_FUNC) &_sydneyPaternity_select_columns_from_cube, 2},
    {"_sydneyPaternity_phenotype_error_model", (DL_FUNC) &_sydneyPaternity_phenotype_error_model, 5},
    {"_sydneyPaternity_sample_phenotype_errors", (DL_FUNC) &_sydneyPaternity_sample_phenotype_errors, 5},
    {"_sydneyPaternity_mendelian_genotype_model", (DL_FUNC) &_sydneyPaternity_mendelian_genotype_model, 6},
    {"_sydneyPaternity_sample_mendelian_genotype", (DL_FUNC) &_sydneyPaternity_sample_mendelian_genotype, 6},
    {"_sydneyPaternity_sample_parentage_and_error_rates_from_joint_posterior", (DL_FUNC) &_sydneyPaternity_sample_parentage_and_error_rates_from_joint_posterior, 8},
    {"_sydneyPaternity_sample_parentage_and_error_rates_from_joint_posterior_alt", (DL_FUNC) &_sydneyPaternity_sample_parentage_and_error_rates_from_joint_posterior_alt, 10},
    {"_sydneyPaternity_sample_parentage_and_error_rates", (DL_FUNC) &_sydneyPaternity_sample_parentage_and_error_rates, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_sydneyPaternity(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
