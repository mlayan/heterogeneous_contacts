#ifndef DEF_UTILS__H
#define DEF_UTILS__H
#include <random>


extern double mIncub;
extern double sdIncub;

extern double mu_f;
extern double sigma_f; 
extern double alpha_f;
extern double m_inc_f;

extern double cf;


// Cumulative probability density distributions
double plnorm(double x, double mean, double sd);
double pgamma(double x, double shape, double scale);
double ptost(double tost, double inc);

// Probability density distributions
double dlnorm(double x, double mean, double sd);
double dgamma(double x, double shape, double scale);
double dtost(double tost, double inc);

// Log probability density distributions
double logdexp(double x, double rate);
double logdlnorm(double x, double mean, double sd);

// Random generators
double rnorm(std::mt19937_64& gen, double mean, double sd);
double runif(std::mt19937_64& gen, double lower_value = 0.0, double upper_value = 1.0);
double rlnorm(std::mt19937_64& gen, double mean, double sd);
double rgamma(std::mt19937_64& gen, double shape, double scale);
double rexp(std::mt19937_64& gen, double rate);

#endif
